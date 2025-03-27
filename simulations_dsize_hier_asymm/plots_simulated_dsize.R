library(tidyverse)
library(spatstat)
library(cmdstanr)
library(Matrix)
library(posterior)
library(ggdist)
library(survival)
devtools::load_all()

path <- "./dsize_hier_asymm/"

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=5, width=8, units="in")
}
figures_folder <- paste0(path,"figures/simulated_dsize/")

out <- readRDS(paste0(path,"analysis_simulated_dsize.rds"))
rmse_tb <- out$rmse_tb
sic_tb <- out$sic_tb
rmse_tb <- rmse_tb %>%
  mutate(coefficient = factor(coefficient)) %>%
  mutate(coefficient = fct_recode(coefficient,
                                  "cohort level"="beta_global",
                                  "subject level"="beta_indiv",
                                  "image level"="beta_local"))
sic_tb <- sic_tb %>%
  mutate(Curve = ifelse(Curve == "Post. Mean","Posterior Mean",Curve))
### RMSE
  
# paper
rmse_tb %>%
  group_by(coefficient,scale,num_pts_per_group) %>%
  summarise(rmse=rvar_mean(rmse)) %>%
  ungroup() %>%
  filter(scale == "all") %>%
  ggplot(aes(num_pts_per_group,ydist=rmse)) +
  stat_halfeye() +
  facet_wrap(~coefficient,scales="free_y") +
  scale_x_continuous(trans="log10") +
  guides(x="axis_logticks") +
  labs(x="# of patients per patient group",y="Average RMSE")
fsave("avg_rmse")

rmse_tb %>%
  group_by(coefficient,scale,num_pts_per_group,images_per_pt) %>%
  summarise(rmse=rvar_mean(rmse)) %>%
  ungroup() %>%
  filter(scale == "all") %>%
  rename(images=images_per_pt,
         coef=coefficient) %>%
  mutate(images = paste0("images: ",images)) %>%
  ggplot(aes(num_pts_per_group,ydist=rmse)) +
  stat_halfeye() +
  facet_grid(coef~images,scales="free_y") +
  scale_x_continuous(trans="log10") +
  guides(x="axis_logticks") +
  labs(x="# of patients per patient group",y="Average RMSE")
fsave("avg_rmse_over_images")

# paper or supplement
# rmse_tb %>%
#   group_by(coefficient,scale,num_pts_per_group,images_per_pt) %>%
#   summarise(rmse=rvar_mean(rmse)) %>%
#   ungroup() %>%
#   filter(scale != "all") %>%
#   filter(coefficient == "beta_global") %>%
#   mutate(scale=factor(scale,levels=c("all","close","mid","far"))) %>%
#   rename(images=images_per_pt,
#          `patients per group`=num_pts_per_group) %>%
#   ggplot(aes(scale,ydist=rmse)) +
#   stat_halfeye() +
#   facet_grid(`patients per group`~images,labeller = label_both) +
#   labs(x="Spatial scale of coefficient",y="Average RMSE") +
#   scale_y_continuous(trans="log10") +
#   guides(y="axis_logticks")
# fsave("avg_rmse_over_scales.png")
# 
# # supplement
# rmse_tb %>%
#   group_by(coefficient,scale,num_pts_per_group,images_per_pt) %>%
#   summarise(rmse=rvar_mean(rmse)) %>%
#   ungroup() %>%
#   filter(scale != "all") %>%
#   filter(coefficient == "beta_indiv") %>%
#   mutate(scale=factor(scale,levels=c("all","close","mid","far"))) %>%
#   rename(images=images_per_pt,
#          `patients per group`=num_pts_per_group) %>%
#   ggplot(aes(scale,ydist=rmse)) +
#   stat_halfeye() +
#   facet_grid(`patients per group`~images,labeller = label_both) +
#   labs(x="Spatial scale of coefficient",y="Average RMSE") +
#   scale_y_continuous(trans="log10") +
#   guides(y="axis_logticks")
# # scale_x_continuous(trans='log10')
# fsave("avg_rmse_over_scales_indiv.png")

# supplement
rmse_tb %>%
  group_by(coefficient,scale,num_pts_per_group,images_per_pt) %>%
  summarise(rmse=rvar_mean(rmse)) %>%
  ungroup() %>%
  filter(scale != "all") %>%
  # filter(coefficient == "beta_local") %>%
  mutate(scale=factor(scale,levels=c("all","close","mid","far"))) %>%
  rename(images=images_per_pt,
         `patients per group`=num_pts_per_group) %>%
  mutate(images = paste0("images: ",images)) %>%
  rename(Scale=scale) %>%
  ggplot(aes(`patients per group`,ydist=rmse,color=Scale,fill=Scale)) +
  stat_halfeye() +
  facet_grid(coefficient~images) +
  labs(x="Number of patients",y="Average RMSE") +
  scale_y_continuous(trans="log10") +
  guides(y="axis_logticks")
# scale_x_continuous(trans='log10')
fsave("avg_rmse_over_scales_all")

### plot example SICs estimate vs ground truth

sic_tb %>%
  filter(sim == 4) %>%
  rename(images=images_per_pt,
         `patients per group`=num_pts_per_group) %>%  
  ggplot(aes(x)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey70") +
  geom_line(aes(y=value,color=Curve),linewidth=1) +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_grid(images~`patients per group`,labeller=label_both,scales="free_y") +
  labs(x="Distance (microns)",y="Log-intensity") +
  theme(legend.position = "top")
fsave("estimated_SIC")
