library(tidyverse)
library(spatstat)
library(cmdstanr)
library(Matrix)
library(posterior)
library(ggdist)
library(survival)
library(SHADE)

# Load utility functions and constants
source("../utils.R")

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_dummy_points/data/"
  grid <- expand.grid(ratio=c(0.5,1,2,5,10),
                      num_points_per_type=c(20,80,150,300,500),
                      sim=1:5)  
} else {
  path <- "./sim_dummy_points/data/"
  grid <- expand.grid(ratio=c(0.5,1,2,5,10),
                      num_points_per_type=c(20,80,150,300,500),
                      sim=1:5)
}

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              )
)

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device = cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./manuscript/images/sim_dummy_points_figures/"

# Load standardized analysis summary
analysis_summary <- readRDS(paste0(path, "analysis_summary.rds"))
rmse_tb <- analysis_summary$rmse_tb
sic_tb <- analysis_summary$sic_tb
rmse_tb <- rmse_tb %>%
  mutate(coefficient = factor(coefficient)) %>%
  mutate(coefficient = fct_recode(coefficient,
                                  "cohort level"="beta_global",
                                  "subject level"="beta_indiv",
                                  "image level"="beta_local")) %>%
  mutate(num_points_per_type = paste0("Count: ",num_points_per_type)) %>%
  mutate(num_points_per_type = factor(num_points_per_type,levels = paste0("Count: ",sort(unique(grid$num_points_per_type)))))

sic_tb <- sic_tb %>%
  mutate(Curve = ifelse(Curve == "Post. Mean","Posterior Mean",Curve)) %>%
  mutate(num_points_per_type = factor(num_points_per_type),
         ratio = factor(ratio)) %>%
  mutate(
    num_points_per_type = fct_relabel(num_points_per_type, ~ paste0("Count: ", .x)),
    ratio = fct_relabel(ratio, ~ paste0("Ratio: ", .x)),
  )

### RMSE

# paper
rmse_tb %>%
  group_by(coefficient,scale,ratio,num_points_per_type) %>%
  summarise(rmse=rvar_mean(rmse)) %>%
  ungroup() %>%
  filter(scale == "all") %>%
  rename(`# of points per type`=num_points_per_type) %>%
  ggplot(aes(ratio,ydist=rmse)) +
  stat_halfeye() +
  facet_grid(coefficient~`# of points per type`,scales="free_y",labeller=label_wrap_gen(width = 20)) +
  labs(x="Ratio of dummy points to real points",y="Average RMSE") +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  guides(x="axis_logticks",y="axis_logticks")
fsave("avg_rmse")

# paper or supplement
rmse_tb %>%
  group_by(coefficient,scale,ratio,num_points_per_type) %>%
  summarise(rmse=rvar_mean(rmse)) %>%
  ungroup() %>%
  filter(scale != "all") %>%
  # filter(coefficient == "beta_global") %>%
  mutate(scale=factor(scale,levels=c("all","close","mid","far"))) %>%
  # rename(`# of points per type`=num_points_per_type) %>%
  ggplot(aes(ratio,ydist=rmse,color=scale,fill=scale)) +
  stat_halfeye() +
  facet_grid(coefficient~num_points_per_type,labeller = label_wrap_gen(width = 20)) +
  labs(x="Ratio of dummy points to real points",y="Average RMSE",color="Scale",fill="Scale") +
  scale_y_continuous(trans="log10",  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(trans="log10") +
  guides(y="axis_logticks",x="axis_logticks") +
  theme(legend.position = "top")
fsave("avg_rmse_over_scales",height=6)

# supplement
# rmse_tb %>%
#   group_by(coefficient,scale,ratio,num_points_per_type) %>%
#   summarise(rmse=rvar_mean(rmse)) %>%
#   ungroup() %>%
#   filter(scale != "all") %>%
#   filter(coefficient == "beta_indiv") %>%
#   mutate(scale=factor(scale,levels=c("all","close","mid","far"))) %>%
#   rename(`# of points per type`=num_points_per_type) %>%
#   ggplot(aes(scale,ydist=rmse)) +
#   stat_halfeye() +
#   facet_grid(ratio~`# of points per type`,labeller = label_both) +
#   labs(x="Spatial scale of coefficient",y="Average RMSE") +
#   scale_y_continuous(trans="log10") +
#   guides(y="axis_logticks")
# # scale_x_continuous(trans='log10')
# fsave("avg_rmse_over_scales_indiv.png")
# 
# # supplement
# rmse_tb %>%
#   group_by(coefficient,scale,ratio,num_points_per_type) %>%
#   summarise(rmse=rvar_mean(rmse)) %>%
#   ungroup() %>%
#   filter(scale != "all") %>%
#   filter(coefficient == "beta_local") %>%
#   mutate(scale=factor(scale,levels=c("all","close","mid","far"))) %>%
#   rename(`# of points per type`=num_points_per_type) %>%
#   ggplot(aes(scale,ydist=rmse)) +
#   stat_halfeye() +
#   facet_grid(ratio~`# of points per type`,labeller = label_both) +
#   labs(x="Spatial scale of coefficient",y="Average RMSE") +
#   scale_y_continuous(trans="log10") +
#   guides(y="axis_logticks")
# # scale_x_continuous(trans='log10')
# fsave("avg_rmse_over_scales_local.png")

### plot example SICs estimate vs ground truth

sic_tb %>%
  filter(x >= MIN_INTERACTION_RADIUS) %>%
  ggplot(aes(x)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey70") +
  geom_line(aes(y=value,color=Curve),linewidth=1) +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_grid(ratio~num_points_per_type,labeller=label_wrap_gen(width=12),scales="free_y") +
  labs(x="Distance (microns)",y="Log-intensity") +
  theme(legend.position = "top")
  # guides(color = guide_legend(override.aes = list( = 2)))

fsave("estimated_SIC",height=6)
