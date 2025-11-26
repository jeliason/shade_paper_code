# Summary Figures

Scripts for generating overview figures that illustrate the SHADE method and summarize results across simulations.

## Scripts

### summary_figure_subplots.R
Generates **Figure 1** (method overview) showing:
- Example point pattern with multiple cell types
- Spatial interaction curves (SICs) at different scales
- Hierarchical model structure illustration

### other_summary_figures.R
Generates **Figure 2** and **Supplementary Figures S1, S2** showing:
- Comparative performance summaries across simulations
- Parameter recovery visualizations
- Uncertainty quantification examples

## Usage

```r
source("summary_figures/summary_figure_subplots.R")  # Figure 1
source("summary_figures/other_summary_figures.R")    # Figures 2, S1, S2
```

## Output

- `summary_figure.pdf` - Main method overview figure
- `group_comp.pdf`, `multi_example.pdf`, `sics_ind.pdf`, etc. - Component plots
- `demo_plots/` - Additional demonstration visualizations
