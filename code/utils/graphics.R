# graphics utilities and shared color scales

# condition color scale
condition_colorscale <- c(
  WT_Std = "#1f77b4",
  WT_CDAA = "#d62728"
)

# expression heatmap color scale
expression_colorscale <- colorspace::diverging_hcl(
  palette = "Green-Orange",
  n = 100
)
