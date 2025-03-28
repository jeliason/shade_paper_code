library(spatstat)
library(tidyverse)
library(posterior)
library(splines2)
library(survival)
library(survminer)
library(cmdstanr)
library(patchwork)
library(SHADE)

# CRC dataset needs to be loaded here
# Example:
# df_raw <- read_csv("path/to/CRC_data.csv")
cat("NOTE: CRC dataset needs to be loaded here. Replace this with the actual data loading code.\n")

df_raw %>%
  count(Spot,type) %>%
  group_by(type) %>%
  summarise(mn = median(n)) %>%
  filter(mn > 30) %>%
  pull(type) %>%
  as.character() -> types_keep
spots <- unique(df_raw$Spot)

# local_sic_df <- local_sic_df %>%
#   mutate(type1 = fct_recode(type1,"TAMs"="CD163+ macros"),
#          type2 = fct_recode(type2,"TAMs"="CD163+ macros"))

# TODO: Include posterior_plots.R functions directly here or load from SHADE package

path <- "./crc_analysis/data/"

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./crc_analysis/figures/"

file_metadata <- paste0(path,"metadata_CRC_type_",make.names(types_keep[3]),".rds")
metadata <- readRDS(file_metadata)
rbfs <- metadata$potentials
x <- seq(0,100,1)
d <- data.frame(x = x)

# Calculate the values for each basis function at each point
for (i in 1:length(rbfs)) {
  d[, paste0("RBF ", i)] <- rbfs[[i]](x)
}

# Reshape data for ggplot
d <- tidyr::pivot_longer(d, cols = starts_with("RBF"), names_to = "rbf", values_to = "y")

# Plot
ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, color = rbf)) +
  ggplot2::geom_line(linewidth=1) +
  labs(x="Distance (microns)",y="Log-intensity",color="Basis")
fsave("pots")


sic_out <- lapply(types_keep,\(type) {

  file_fit <- paste0(path,"fit_CRC_type_",make.names(type),".rds")
  file_metadata <- paste0(path,"metadata_CRC_type_",make.names(type),".rds")
  metadata <- readRDS(file_metadata)


  spots <- metadata$spots
  types <- metadata$types
  potentials <- metadata$potentials
  coef_names <- metadata$coef_names
  ix_b0 <- grep("beta0",coef_names,value = FALSE)
  combo_names <- coef_names[-ix_b0]

  pt_data %>%
    filter(Spot %in% spots) %>%
    distinct(Patient,Group,OS,OS_Censor,Age,Sex,p_TNM) %>%
    mutate(p_TNM = factor(p_TNM,levels=c("3A","3B","3C","4"),ordered=TRUE)) %>%
    mutate(across(c(Age),~as.vector(scale(.x,center=TRUE)))) -> pt_df

  num_pot <- length(potentials)
  num_types <- length(types)
  num_samples <- length(spots)
  num_combos <- num_types - 1
  num_indiv <- nrow(pt_df)


  fit <- readRDS(file_fit)
  draws <- as_draws_rvars(fit$draws())

  beta_global <- draws$beta_global
  beta_indiv <- draws$beta_indiv
  beta_local <- t(draws$beta_local)
  rownames(beta_global) <- coef_names
  rownames(beta_indiv) <- coef_names
  rownames(beta_local) <- coef_names

  # global-indiv SICs
  get_SIC_df <- function(type1,type2,exponentiate=FALSE) {
    x_seq <- seq(0,100,1)
    x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)
    ix <- grep(paste0("_",type1,"_",type2),rownames(beta_global),fixed = TRUE)
    b_g <- as.matrix(beta_global)[ix,]

    lp_g <- x_des %*% b_g
    lp <- lapply(1:num_indiv,\(s) {
      as.vector(x_des %*% as.vector(beta_indiv[ix,s]))
    }) %>%
      do.call(cbind,.)

    lp <- cbind(lp,lp_g)
    if(exponentiate) {
      lp <- exp(lp)
    }

    lev_g <- levels(factor(pt_df$Group))
    colnames(lp) <- c(paste0("Patient: ",pt_df$Patient),paste0("Global: ",lev_g))

    lp %>%
      as.data.frame() %>%
      mutate(x=x_seq) %>%
      mutate(across(-x,list(
        mn = ~as.vector(E(.)),
        lo = ~as.vector(quantile(.,0.1)),
        hi = ~as.vector(quantile(.,0.9))
      )),.keep="unused") %>%
      pivot_longer(
        cols = ends_with(c("mn", "lo", "hi")),  # Select the new columns
        names_to = c("lp", "transformation"),      # Split into two parts: original name and transformation type
        names_pattern = "(.*)_(.*)"                    # Regex to separate the column name from the transformation
      ) %>%
      pivot_wider(names_from = transformation,values_from=value) %>%
      mutate(is_global = ifelse(str_detect(lp,"Global"),"Global","Individual")) %>%
      separate(lp,c("level","Patient"),sep = ": ",remove = FALSE) %>%
      mutate(global_level = ifelse(Patient %in% unique(pt_df$Group),Patient,NA)) %>%
      mutate(Patient = ifelse(Patient %in% unique(pt_df$Group),NA,Patient)) %>%
      mutate(Patient = as.numeric(Patient)) %>%
      left_join(pt_df %>% distinct(Patient,Group),by="Patient") %>%
      mutate(Group = ifelse(is.na(Group),global_level,Group)) %>%
      mutate(type1=type1,type2=type2)
  }

  sic_df <- expand_grid(type1=type,type2=types[-which(types == type)]) %>%
    pmap(\(type1,type2) {
      cat(type1,", ", type2, "\n")
      # print(type1)
      get_SIC_df(type1,type2,exponentiate = TRUE)
    }) %>% bind_rows()

  # indiv-local SICs
  get_local_SIC_df <- function(type1,type2,exponentiate=FALSE) {
    x_seq <- seq(0,100,1)
    x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)
    ix <- grep(paste0("_",type1,"_",type2),rownames(beta_global),fixed = TRUE)

    lp_indiv <- lapply(1:num_indiv,\(s) {
      as.vector(x_des %*% as.vector(beta_indiv[ix,s]))
    }) %>%
      do.call(cbind,.)

    lp_local <- lapply(1:num_samples,\(s) {
      as.vector(x_des %*% as.vector(beta_local[ix,s]))
    }) %>%
      do.call(cbind,.)

    lp <- cbind(lp_indiv,lp_local)

    if(exponentiate) {
      lp <- exp(lp)
    }
    lp <- as.data.frame(lp)

    colnames(lp) <- c(paste0("Patient: ",pt_df$Patient),paste0("Spot: ",pt_data$Spot))

    lp_df <- lp %>%
      as.data.frame() %>%
      mutate(x=x_seq) %>%
      mutate(across(-x,list(
        mn = ~as.vector(E(.)),
        lo = ~as.vector(quantile(.,0.1)),
        hi = ~as.vector(quantile(.,0.9))
      )),.keep="unused") %>%
      pivot_longer(
        cols = ends_with(c("mn", "lo", "hi")),  # Select the new columns
        names_to = c("lp", "transformation"),      # Split into two parts: original name and transformation type
        names_pattern = "(.*)_(.*)"                    # Regex to separate the column name from the transformation
      ) %>%
      pivot_wider(names_from = transformation,values_from=value) %>%
      mutate(is_indiv = ifelse(str_detect(lp,"Patient"),"Individual","Image")) %>%
      separate(lp,c("level","ID"),sep = ": ",remove = FALSE)

    lp_df %>%
      mutate(Patient = ifelse(level == "Patient",ID,NA)) %>%
      mutate(indiv_level = ifelse(level == "Patient",NA,ID)) %>%
      left_join(pt_data %>% distinct(Patient,Spot) %>% rename(pt = Patient),by=c("indiv_level"="Spot")) %>%
      mutate(Patient = ifelse(is.na(Patient),pt,Patient)) %>%
      mutate(type1=type1,type2=type2) %>%
      select(-c(indiv_level,pt))
  }

  local_sic_df <- expand_grid(type1=type,type2=types[-which(types == type)]) %>%
    pmap(\(type1,type2) {
      cat(type1,", ", type2, "\n")
      # print(type1)
      get_local_SIC_df(type1,type2,exponentiate = TRUE)
    }) %>% bind_rows()

  list(sic_df=sic_df,local_sic_df=local_sic_df)
})

sic_df <- lapply(sic_out,\(o) o$sic_df) %>% bind_rows()
local_sic_df <- lapply(sic_out,\(o) o$local_sic_df) %>% bind_rows()

saveRDS(sic_df,paste0(path,"sic_df.rds"))
saveRDS(local_sic_df,paste0(path,"local_sic_df.rds"))

local_sic_df <- readRDS(paste0(path,"local_sic_df.rds"))
sic_df <- readRDS(paste0(path,"sic_df.rds"))

local_sic_df$Patient <- factor(local_sic_df$Patient,levels=1:35)


sic_df <- sic_df %>%
  mutate(type1 = ifelse(type1 == "CD163+ macros","TAMs",type1),
         type2 = ifelse(type2 == "CD163+ macros","TAMs",type2))
  
# types_filt <- types[c(1,2,5)]
sic_df %>%
  filter(type1 == "tumor cells" & type2 == "TAMs") %>%
  # filter(type1 %in% types_filt & type2 %in% types_filt) %>%
  # filter(type1 != "B cells" & type2 != "B cells") %>%
  mutate(is_global = ifelse(is_global == "Global","Cohort",is_global)) %>%
ggplot(aes(x)) +
  geom_line(aes(y=log(mn),color=Group,group=lp,linetype=is_global,linewidth=is_global,alpha=is_global)) +
  # geom_ribbon(aes(ymin=lo,ymax=hi,fill=is_global,alpha=is_global)) +
  geom_hline(yintercept=0,linetype="dotted") +
  labs(x="Distance (microns)",y="Log-intensity",linetype="") +
  scale_linetype_manual(values=c("solid","dashed")) +
  scale_alpha_manual(values=c(1,0.3)) +
  scale_linewidth_manual(values=c(2,1)) +
  # facet_grid(type1~type2,switch="y") +
  # scale_color_manual(values = as.vector(pals::glasbey())) +
  # scale_fill_manual(values = as.vector(pals::glasbey())) +
  guides(alpha="none",group="none",fill="none",linewidth="none")
fsave("sic_CRC")
  
sic_df %>%
  filter(is_global == "Global") %>%
  ggplot(aes(x)) +
  geom_line(aes(y=log(mn),color=Group,group=lp,linetype=is_global,linewidth=is_global,alpha=is_global)) +
  # geom_ribbon(aes(ymin=lo,ymax=hi,fill=is_global,alpha=is_global)) +
  geom_hline(yintercept=0,linetype="dotted") +
  labs(x="Distance (microns)",y="Log-intensity",linetype="") +
  scale_linetype_manual(values=c("solid","dashed")) +
  scale_alpha_manual(values=c(1,0.3)) +
  scale_linewidth_manual(values=c(1,1)) +
  facet_grid(type1~type2,switch="y",labeller=label_wrap_gen(width = 8)) +
  # scale_color_manual(values = as.vector(pals::glasbey())) +
  # scale_fill_manual(values = as.vector(pals::glasbey())) +
  guides(alpha="none",group="none",fill="none",linewidth="none",linetype="none") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
fsave("sic_CRC_all",height = 10,width=11)

# variance decomposition
# Step 1: compute cohort mean per x / cell type pair
cohort_means <- sic_df %>%
  filter(level == "Global") %>%
  group_by(x, Group, type1, type2) %>%
  summarise(mn_cohort = median(mn, na.rm = TRUE), .groups = "drop")

# Step 2: merge with patient-level data
patient_merged <- sic_df %>%
  filter(level != "Global") %>%
  left_join(cohort_means, by = c("x", "Group", "type1", "type2"))

# Step 3: compute variance of patient deviation from cohort mean, per x
var_within_cohort_x <- patient_merged %>%
  group_by(x, type1, type2) %>%
  summarise(
    var_within_cohort = mad(mn - mn_cohort, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: average across x
var_within_cohort <- var_within_cohort_x %>%
  group_by(type1, type2) %>%
  summarise(var_within_cohort = median(var_within_cohort, na.rm = TRUE), .groups = "drop")

# Step 1: compute patient mean per x / cell type pair
patient_means <- local_sic_df %>%
  filter(level == "Patient") %>%
  group_by(x, Patient, type1, type2) %>%
  summarise(mn_patient = median(mn, na.rm = TRUE), .groups = "drop")

# Step 2: merge with image-level data
image_merged <- local_sic_df %>%
  filter(level != "Patient") %>%
  left_join(patient_means, by = c("x", "Patient", "type1", "type2"))

# Step 3: compute variance of image deviation from patient mean, per x
var_within_patient_x <- image_merged %>%
  group_by(x, type1, type2) %>%
  summarise(
    var_within_patient = mad(mn - mn_patient, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: average across x
var_within_patient <- var_within_patient_x %>%
  group_by(type1, type2) %>%
  summarise(var_within_patient = median(var_within_patient, na.rm = TRUE), .groups = "drop")

# Combine both for plotting
var_combined <- var_within_cohort %>%
  left_join(var_within_patient, by = c("type1", "type2")) %>%
  pivot_longer(cols = starts_with("var_"), names_to = "level", values_to = "variance")

# Clean up labels for plotting
# var_combined <- var_combined %>%
#   mutate(level = recode(level,
#                         var_within_cohort = "Between-patient (within cohort)",
#                         var_within_patient = "Between-image (within patient)"))

# Plot heatmaps
p1 <- var_combined %>%
  filter(level == "var_within_cohort") %>%
ggplot(aes(x = type1, y = type2, fill = variance)) +
  geom_tile() +
  # facet_wrap(~level, ncol = 1) +
  scale_fill_viridis_c(option = "plasma", name = "Within-cohort average MAD",n.breaks=4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Target cell type",
    y = "Source cell type"
  ) +
  theme(legend.position = "top")

p2 <- var_combined %>%
  filter(level == "var_within_patient") %>%
  ggplot(aes(x = type1, y = type2, fill = variance)) +
  geom_tile() +
  # facet_wrap(~level, ncol = 1) +
  scale_fill_viridis_c(option = "plasma", name = "Within-patient average MAD",n.breaks=4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Target cell type",
    y = "Source cell type",
  ) +
  theme(legend.position="top",
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y = element_blank())

p1 + p2
fsave("sic_mad",width=10)

d1 <- local_sic_df %>%
  mutate(mn = log(mn)) %>%
  filter(Patient %in% c(2,19,30)) %>%
  filter(type1 == "tumor cells" & type2 == "hybrid E/M") %>%
  mutate(t1t2 = paste(type1,"-->",type2))

d2 <- local_sic_df %>%
  mutate(mn = log(mn)) %>%
  filter(Patient %in% c(11,28,31)) %>%
  filter(type1 == "CAFs" & type2 == "memory CD4+ T") %>%
  mutate(t1t2 = paste(type1,"-->",type2))

bind_rows(d1,d2) %>%
  ggplot(aes(x)) +
  geom_line(aes(y=mn,color=Patient,group=lp,linetype=is_indiv,linewidth=is_indiv,alpha=is_indiv)) +
  # geom_ribbon(aes(ymin=lo,ymax=hi,fill=is_global,alpha=is_global)) +
  geom_hline(yintercept=0,linetype="dotted") +
  labs(x="Distance (microns)",y="Log-intensity",linetype="") +
  scale_linetype_manual(values=c("dashed","solid")) +
  scale_alpha_manual(values=c(0.3,1)) +
  scale_linewidth_manual(values=c(1,2)) +
  facet_wrap(~t1t2) +
  # facet_grid(type1~type2,switch="y") +
  # scale_color_manual(values = as.vector(pals::trubetskoy())) +
  # scale_fill_manual(values = as.vector(pals::glasbey())) +
  guides(alpha="none",group="none",fill="none",linewidth="none")
fsave("sic_local_example")

# local_sic_df %>%
#   mutate(mn = log(mn)) %>%
#   filter(Patient %in% c(2,19,30)) %>%
#   filter(type1 == "tumor cells" & type2 == "hybrid E/M") %>%
#   
#   # filter(type1 %in% types_filt & type2 %in% types_filt) %>%
#   # filter(type1 != "B cells" & type2 != "B cells") %>%
#   ggplot(aes(x)) +
#   geom_line(aes(y=mn,color=Patient,group=lp,linetype=is_indiv,linewidth=is_indiv,alpha=is_indiv)) +
#   # geom_ribbon(aes(ymin=lo,ymax=hi,fill=is_global,alpha=is_global)) +
#   geom_hline(yintercept=0,linetype="dotted") +
#   labs(x="Distance (microns)",y="Log-intensity",linetype="") +
#   scale_linetype_manual(values=c("dashed","solid")) +
#   scale_alpha_manual(values=c(0.3,1)) +
#   scale_linewidth_manual(values=c(1,2)) +
#   # facet_wrap(~Patient) +
#   # facet_grid(type1~type2,switch="y") +
#   # scale_color_manual(values = as.vector(pals::trubetskoy())) +
#   # scale_fill_manual(values = as.vector(pals::glasbey())) +
#   guides(alpha="none",group="none",fill="none",linewidth="none")
# fsave("sic_local_tumor_EM")


# 1. Create a background dataframe with a vertical gradient
gradient_df <- expand.grid(
  x = seq(0, 100, length.out = 200),
  y = seq(-0.3, 0.3, length.out = 200)
)

max_y <- max(abs(gradient_df$y))

# Use soft colors: blue top, rose bottom
gradient_df <- gradient_df %>%
  mutate(
    fill = case_when(
      y >= 0 ~ alpha("#a6cbe3", y / max_y),  # soft blue
      y <  0 ~ alpha("#f5b8b4", abs(y) / max_y)  # soft rose
    )
  )

# 2. Prep your SIC difference data
diff_df <- sic_df %>%
  filter(is_global == "Global") %>%
  mutate(mn = log(mn)) %>%
  mutate(
    type1 = recode(type1, "granulocytes" = "gran.", "vasculature" = "vasc."),
    type2 = recode(type2, "granulocytes" = "gran.", "vasculature" = "vasc.")
  ) %>%
  group_by(x, type1, type2) %>%
  summarise(
    diff = mn[Group == "CLR"] - mn[Group == "DII"],
    diff_lo = lo[Group == "CLR"] - lo[Group == "DII"],
    diff_hi = hi[Group == "CLR"] - hi[Group == "DII"],
    .groups = "drop"
  ) %>%
  mutate(
    type1 = paste0("Target: ", type1),
    type2 = paste0("Source: ", type2)
  )

# 3. Plot
ggplot() +
  # Gradient background
  geom_raster(data = gradient_df, aes(x = x, y = y, fill = fill)) +
  scale_fill_identity() +
  
  # Foreground plot layers
  geom_line(data = diff_df, aes(x = x, y = diff), color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "#d88c3d") +
  
  # Labels
  labs(
    x = "Distance (microns)",
    y = "Difference in SIC between Groups",
    linetype = ""
  ) +
  
  # Facet grid
  facet_grid(type2 ~ type1, labeller = label_wrap_gen(width = 8)) +
  
  # Axis breaks
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 4) +
  
  # Clean up guides
  guides(alpha = "none", group = "none", fill = "none", linewidth = "none", linetype = "none") +
  
  # Theme tweaks
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  )
fsave("sic_CRC_diff",height = 8,width=9)
  
# compare cell type combos of SICs, in groups
sic_df %>%
  mutate(combo = paste0(pmin(type1,type2),"\n",pmax(type1,type2))) %>%
  # mutate(direction=ifelse(pmin(type1,type2) == type1,"Type 1 is Receiver","Type 1 is Effector")) %>%
  mutate(direction=ifelse(pmin(type1,type2) == type1,TRUE,FALSE)) %>%
  
  # filter(combo == "tumor cells <--> vasculature") %>%
  filter(is_global == "Global") %>%
  # filter(Group == "CLR") %>%
  filter(type1 != type2) %>%
  group_by(x,combo,Group) %>%
  summarise(diff=mn[direction] - mn[!direction],
            diff_lo=lo[direction] - lo[!direction],
            diff_hi=hi[direction] - hi[!direction]) %>%
  ungroup() %>%
  ggplot(aes(x)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_ribbon(aes(ymin=diff_lo,ymax=diff_hi,fill=Group),alpha=0.2) +
  geom_line(aes(y=diff,color=Group)) + 
  facet_wrap(~combo,nrow=4) +#,labeller=label_wrap_gen(width=15)) +
  labs(x="Distance (microns)",y="Difference in SICs between Asymmetric Pairs") +
  theme(legend.position = "top")
fsave("combo_compare",height=6,width=10)

# sic_df %>%
#   mutate(combo = paste0(pmin(type1,type2)," <--> ",pmax(type1,type2))) %>%
#   mutate(direction=ifelse(pmin(type1,type2) == type1,TRUE,FALSE)) %>%
#   # filter(combo == "tumor cells <--> vasculature") %>%
#   filter(is_global == "Global") %>%
#   filter(Group == "DII") %>%
#   filter(type1 != type2) %>%
#   group_by(x,combo) %>%
#   summarise(diff=mn[direction] - mn[!direction],
#             diff_lo=lo[direction] - lo[!direction],
#             diff_hi=hi[direction] - hi[!direction]) %>%
#   ungroup() %>%
#   ggplot(aes(x)) +
#   geom_hline(yintercept = 0,linetype="dotted") +
#   geom_ribbon(aes(ymin=diff_lo,ymax=diff_hi),alpha=0.2) +
#   geom_line(aes(y=diff)) +
#   facet_wrap(~combo,ncol=4) +
#   labs(x="Distance",y="Difference in Relative Intensity between Asymmetric Pairs")
# fsave("combo_compare_DII.png")


# heatmaps of significant beta coefficients

# AUCs for different cell types, at image levels. Can then aggregate into patients and groups.
spots_ix <- 1:length(spots)
dats <- df_raw %>%
  dplyr::filter(Spot %in% spots[spots_ix]) %>%
  group_by(Spot) %>%
  group_map(~{
    .x <- .x %>%
      filter(type %in% types_keep) %>%
      droplevels()
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

n_dummy <- 30
Qs <- lapply(pats,\(pat) make_quadrature(pat,n_dummy = n_dummy,dist = "grid"))

lapply(types_keep,\(type) {

  lapply(Qs,\(Q) {
    log(intensity(Q$dummy)) %>%
      enframe() %>%
      right_join(tibble(name=marks(Q)),by="name") %>%
      filter(name == type) %>%
      pull(value)
  }) %>% unlist() -> offset
  
  # data_lists <- lapply(1:length(Qs),\(i) {
  #   data_list <- make_acyclic_data(Qs[[i]],potentials,type,verbose=FALSE)
  #   # print(i)
  #   data_list
  # })
  # saveRDS(data_lists,paste0(path,"data_lists_auc_grid_type",make.names(type),".rds"))
  
  data_lists <- readRDS(paste0(path,"data_lists_auc_grid_type",make.names(type),".rds"))
  
  x_cells <- do.call(rbind,lapply(data_lists,\(data_list) data_list$data))
  nms <- colnames(x_cells)
  # Define columns to set to zero

  # Create a diagonal matrix with ones, except zeros at specified indices
  D <- Matrix::Diagonal(length(nms))
  x_cells <- x_cells %*% D
  is_cell <- unlist(lapply(data_lists,\(data_list) data_list$response))
  keep_ix <- which(is_cell == 0)
  
  # we need to restrict this to cell types not equal to the target cell type
  # otherwise we have data leakage
  
  x_cells <- x_cells[keep_ix,]
  
  sample_id <- lapply(1:length(data_lists),\(i) {
    rep(i,length(data_lists[[i]]$response))
  }) %>% unlist()
  
  sample_id <- sample_id[keep_ix]
  
  y_start_stop<- lapply(1:length(data_lists),\(i) {
    idx <- sort(which(sample_id == i))
    c(idx[1],idx[length(idx)])
  }) %>% do.call(rbind,.)
  
  offset <- offset[keep_ix]
  
  grid_size <- n_dummy^2 + 4
  
  num_indiv <- 35
  
  n_samples <- length(spots_ix)
  
  
  
  data_gq <- list(
    num_indiv = num_indiv,
    num_types = num_types,
    num_pot = num_pot,
    num_pt_groups = 2,
    n_cells = nrow(x_cells),
    d_cells = ncol(x_cells),
    grid_size = grid_size,
    oset = -offset,
    n_samples = n_samples,
    y_start_stop = y_start_stop,
    x_cells = as.matrix(x_cells)
  )

  file_fit <- paste0(path,"fit_CRC_type_",make.names(type),".rds")
  fit <- readRDS(file_fit)
  mod_gq <- cmdstan_model("scripts/BISTRO_gq.stan")
  
  draws_fit <- fit$draws()
  # draws_fit <- draws_fit[1:100,]
  draws_fit <- as_draws_matrix(t(as.matrix(apply(draws_fit,2,mean))))
  fit_gq <- mod_gq$generate_quantities(draws_fit, data = data_gq, seed = 123)
  
  
  draws_gq <- as_draws_rvars(fit_gq$draws())
  p_pred <- E(draws_gq$p_pred)
  
  dim(p_pred)

  print(type)
  # spot_ix <- which(spots == "22_A")
  spdf <- lapply(1:length(Qs),\(i) {
    Qs[[i]]$dummy %>% 
      as.data.frame() %>%
      filter(marks == type) %>%
      mutate(spot=spots[i]) %>%
      mutate(log_pred=p_pred[i,])
    }) %>%
    bind_rows()

  spdf %>%
    group_by(spot,x) %>%
    mutate(count_x=n()) %>%
    ungroup() %>%
    filter(count_x == n_dummy) %>%
    group_by(spot,y) %>%
    mutate(count_y=n()) %>%
    ungroup() %>%
    filter(count_y == n_dummy) %>%
    select(-c(count_x,count_y)) -> spdf
  
  spdf
  
}) %>%
  bind_rows() -> spdf

saveRDS(spdf,paste0(path,"spdf.rds"))
spdf <- readRDS(paste0(path,"spdf.rds"))

spot_ix <- which(spots == "22_A")
spdf %>%
  filter(spot == spots[[spot_ix]]) %>%
  ggplot(aes(x,y)) +
  geom_raster(aes(fill=log_pred)) +
  geom_point(data=as.data.frame(pats[[spot_ix]]),size=0.5) +
  facet_wrap(~marks,ncol=4) +
  scico::scale_fill_scico(palette = "imola") +
  labs(fill="Log-intensity",x="X (microns)",y="Y (microns)") +
  theme(legend.position = "top")
fsave("example_pred")

types <- types_keep
aucs_df <- lapply(1:length(spots),\(spot_ix) {
  print(spot_ix)
  
  lapply(types_keep,\(type) {
    # print(type)
    covar <- spdf %>%
      filter(marks == type) %>%
      filter(spot == spots[spot_ix]) %>%
      select(-marks,-spot) %>%
      as.data.frame() %>%
      as.im()
    
    pp <- pats[[spot_ix]]
    pp <- subset(pp,marks==type)
    pp$marks <- NULL
    
    if(pp$n > 0) {
      auc(pp,covar)
    } else {
      NA
    }
  }) %>%
    setNames(.,types) %>%
    as_tibble() %>%
    mutate(Spot = spots[spot_ix])

}) %>%
  bind_rows()

aucs_df %>%
  pivot_longer(-Spot) %>%
  group_by(Spot) %>%
  summarise(mn=mean(value)) %>%
  arrange(desc(mn))

aucs_df %>%
  left_join(pt_data %>% select(Spot,Patient,Group)) %>%
  pivot_longer(-c(Spot,Patient,Group)) %>%
  ggplot(aes(x=Group,y=value)) +
  geom_boxplot() +
  geom_jitter(color="orange") +
  facet_wrap(~name) +
  labs(y="AUC")
fsave("auc_boxplot_groups")



