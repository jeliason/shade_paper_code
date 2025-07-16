library(tidyverse)

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./sim_shade_gcross/gx_figures/"

grid <- expand.grid(
  t_density = c("high", "low"), 
  tumor_density = c("high", "low"),
  num_images = c(1,2,3),
  sim_rep = 1:30
) %>%
  mutate(sim_idx = 1:nrow(.))

df <- pmap(grid,\(t_density,tumor_density,num_images,sim_rep,sim_idx) {
  file_out <- paste0("./sim_shade_gcross/data/sim_",sim_idx,".rds")
  
  out <- readRDS(file_out)
  
  tibble(shade_image_power=out$shade_image_power,
         simple_shade_image_power=out$simple_shade_image_power,
         gcross_image_power=out$gcross_image_power,
         sim_idx=sim_idx)
}) %>%
  bind_rows() %>%
  left_join(grid)

df %>%
  mutate(num_images = factor(num_images)) %>%
  rename(Gcross=gcross_image_power,SHADE=shade_image_power,Flat=simple_shade_image_power) %>%
  pivot_longer(c(SHADE,Flat,Gcross),names_to="Method") %>%
  mutate(Method=factor(Method,levels=c("Gcross","Flat","SHADE"))) %>%
  rename(`Tumor density`=tumor_density,`T cell density`=t_density) %>%
  mutate(value = value / 3) %>%
  ggplot(aes(x=num_images,y=value,fill=Method,color=Method)) +
  geom_boxplot() +
  facet_grid(`T cell density` ~ `Tumor density`,labeller=label_both) +
  labs(x="Number of images per patient",y="Average proportion of correct detections")
fsave("prop_detections")
