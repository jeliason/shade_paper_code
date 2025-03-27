library(tidyverse)
library(spatstat)
library(cmdstanr)
library(Matrix)
library(posterior)
library(ggdist)
library(survival)
library(kableExtra)
devtools::load_all()

path <- "./no_hier/"

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=5, width=8, units="in")
}
figures_folder <- paste0(path,"figures/simulated_no_hier/")

out <- readRDS(paste0(path,"analysis_simulated_no_hier.rds"))
rmse_tb <- out$rmse_tb
sic_tb <- out$sic_tb
sic_tb <- sic_tb %>%
  mutate(Curve = ifelse(Curve == "Post. Mean (Hier)","Posterior Mean (Hier)",Curve),
         Curve = ifelse(Curve == "Post. Mean (No Hier)","Posterior Mean (No Hier)",Curve))
### RMSE
  
# paper
rmse_table <- rmse_tb %>%
  group_by(scale, hier) %>%
  reframe(
    Mean_RMSE = mean(rvar_mean(rmse, na.rm = TRUE)),
    Lower = mean(rvar_quantile(rmse, probs=0.025, na.rm = TRUE)),
    Upper = mean(rvar_quantile(rmse, probs=0.975, na.rm = TRUE))
  ) %>%
  mutate(
    RMSE_Text = sprintf("%.3f (%.3f, %.3f)", Mean_RMSE, Lower, Upper)  # Format as "mean (lower, upper)"
  ) %>%
  select(hier, scale, RMSE_Text) %>%
  pivot_wider(names_from = scale, values_from = RMSE_Text)

# Create a LaTeX table
kable(rmse_table, format = "latex", booktabs = TRUE, escape = FALSE,  # Prevent LaTeX formatting issues
      col.names = c("Hierarchical", "Scale: All", "Scale: Small", "Scale: Medium", "Scale: Large")) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  column_spec(2:5, width = "4cm") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1, background = "#F0F0F0") %>%  # Manually set gray for first data row

  footnote(general = "RMSE values with 95% credible intervals in parentheses. NA indicates missing data.")

### plot example SICs estimate vs ground truth

ggplot(sic_tb, aes(x = x, y = Value,linewidth=Curve)) +
  geom_ribbon(
    data = sic_tb %>% filter(Curve != "True"),  # Only apply ribbon to posterior means
    aes(ymin = Lower, ymax = Upper, fill = Curve),
    alpha = 0.2
  ) +
  geom_line(aes( color = Curve,linetype = Curve)) +  # Thicker lines for clarity
  geom_hline(yintercept=0,linetype="dotted") +
  facet_wrap(~sim) +
  # Set the same colors for posterior means and a distinct color for the true curve
  scale_linetype_manual(values = c(
    "Posterior Mean (Hier)" = "dashed",
    "Posterior Mean (No Hier)" = "dashed",
    "True" = "solid"
  )) +
  scale_color_manual(values = c(
    "Posterior Mean (Hier)" = "#F8766D",
    "Posterior Mean (No Hier)" = "#7CAE00",
    "True" = "#00BFC4"
  )) +
  # Make the fills match the color of the posterior means, and exclude fill for the true curve
  scale_fill_manual(values = c(
    "Posterior Mean (Hier)" = "#F8766D",
    "Posterior Mean (No Hier)" = "#7CAE00",
    "True" = NA  # No fill for the true curve
  )) +
  scale_linewidth_manual(values = c(
    "Posterior Mean (Hier)" = 0.8,
    "Posterior Mean (No Hier)" = 0.8,
    "True" = 1.5  # Larger linewidth for the true curve
  )) +
  labs(x = "Distance (microns)", y = "Log-intensity", color = "Curve", linetype = "Curve", fill = "Curve") +
  guides(fill="none",linetype="none",linewidth="none") +
  theme(strip.text = element_blank(),
        legend.position = "top")

fsave("estimated_SIC")
