library(tidyverse)
library(spatstat)
library(Matrix)
library(cmdstanr)
library(SHADE)

source("crc_analysis/utils.R")
theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./crc_analysis/CRC_analysis_paper/"

# CRC dataset needs to be loaded here
# Example:
df_raw <- read_csv("crc_analysis/data/CRC_cleaned.csv")
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

df_raw <- df_raw %>%
  dplyr::mutate(type = as.factor(type)) %>%
  mutate(type = fct_recode(type,"CAFs"="smooth muscle","hybrid E/M"="stroma","TAMs"="CD163+ macros","CTLs"="CD8+ T cells"))

targets <- c(
  "CTLs",
  "memory CD4+ T",
  "granulocytes"
)

sources <- c(
  "TAMs",
  "CAFs",
  "vasculature",
  "hybrid E/M",
  "tumor cells"
)

spots <- pt_data %>%
  pull(Spot)

path <- "./crc_analysis/data/"
n_dummy <- 1e3

types_keep <- c(targets,sources)

num_types <- length(types_keep)

num_pot <- 3
potentials <- make_rbfs(n_basis_functions = num_pot, max_dist = 75, basis_function_sigma = 15)
plot_rbfs(seq(0,80,1),potentials)

dats <- df_raw %>%
  dplyr::filter(Spot %in% spots) %>%
  group_by(Spot) %>%
  group_map(~{
    .x <- .x %>%
      filter(type %in% types_keep) %>%
      droplevels() %>%
      mutate(Spot = .y$Spot)
  })

pats <- lapply(1:length(dats),\(i) {
  df <- dats[[i]]
  tryCatch({
    pat <- make_pat(df$X,df$Y,factor(df$type,levels=types_keep))
    sq_W <- owin(xrange = c(min(df$X),max(df$X)),yrange=c(min(df$Y),max(df$Y)))
    Window(pat) <- sq_W
    pat
  },error=\(e) {
    print(i)
    stop(e)
  })
})

# Distance grid
dists = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)

# Run Gcross estimation across type pairs
Gx <- lapply(1:length(pats),\(i) {
  pp <- pats[[i]]
  res <- expand_grid(t1 = sources, t2 = targets) %>%
    pmap_dfr(function(t1, t2) {
      tryCatch({
        f <- Gcross(pp, i = t1, j = t2)
  
        ix2 <- unlist(lapply(dists, function(x) which.max(f$r[f$r <= x])))
        
        f %>%
          as.data.frame() %>%
          mutate(type1 = t1, type2 = t2) %>%
          slice(ix2) %>%
          mutate(dist = dists)
      }, error = function(e) {
        print(e)
        # print(paste0(t1," ",t2))
        NULL
      })
    }) %>%
    mutate(Spot = unique(dats[[i]]$Spot))
}) %>%
  bind_rows()

df <- Gx %>%
  left_join(pt_data) %>%
  filter(dist %in% c(20,40,60)) %>%
  mutate(Group = factor(Group)) %>%
  group_by(type1,type2,dist) %>%
  group_modify(~{
    glm(Group ~ km, family="binomial",data=.x) %>%
    # glm(Group ~ iso, family="binomial",data=.x) %>%
      broom::tidy() %>%
      filter(term != "(Intercept)")
  })

df %>%
  rename(target=type2,source=type1) %>%
  group_by(target) %>%
  mutate(p.adj=p.adjust(p.value)) %>%
  ungroup() %>%
  mutate(significant = p.adj < 0.05) %>%
  mutate(dist=paste0(dist," microns")) %>%
  rename(Distance=dist) %>%
  ggplot(aes(source,target,fill=estimate)) +
  geom_tile() +
  geom_text(aes(label = ifelse(significant, "*", "")), 
            color = "black", size = 5) +
  scico::scale_fill_scico(palette="bam",midpoint=0) +
  facet_wrap(~Distance,nrow=3,labeller=label_both) +
  labs(x="Source cell type",y="Target cell type",fill="log-OR")
fsave("gx_comparison")

df %>%
  rename(target=type2,source=type1) %>%
  group_by(target) %>%
  mutate(p.adj=p.adjust(p.value)) %>%
  ungroup() %>%
  mutate(significant = p.adj < 0.05) %>%
  mutate(dist=paste0(dist," microns")) %>%
  rename(Distance=dist) %>%
  print(n=nrow(.))
