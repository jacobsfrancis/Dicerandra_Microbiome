# Reorganized R scripts

Run the scripts in numeric order.

1. `01_Data_Cleaning.R` reads raw data, assigns compartments, restricts the
   primary analysis to sites represented in all seven compartments, calculates
   Bray–Curtis distances, and saves cleaned objects.
2. `02_Statistical_Analyses.R` runs the coordination mixed models, site
   PERMANOVAs and dispersion tests, and the compartment-by-site-status model.
3. Scripts `03`–`10` each generate one final figure.

The first two scripts save reusable R objects in:

`../Output/Poster_Common_5_Sites/Intermediate/`

The figure scripts read those objects, so they do not repeat data preparation
or model fitting.

All section headings use RStudio Outline-compatible formatting:

`# 1) Major section ####`

`## 1a) Subsection ####`
