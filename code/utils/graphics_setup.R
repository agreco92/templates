library(tidyverse)
library(ggthemes)
library(patchwork)

theme_set(theme_bw(base_size = 12))

theme_update(
  plot.title = element_text(size = rel(1), hjust = 0.5),
  plot.subtitle = element_text(size = rel(0.8), hjust = 0.5),
  legend.text = element_text(size = rel(0.8)),
  legend.title = element_text(size = rel(1)),
  axis.title = element_text(size = rel(0.8)),
  axis.text = element_text(size = rel(.6)),
  #plot.margin = margin(1,1,1,1),
  line = element_line(linewidth = .5 / .pt / .75),
  axis.ticks.length = unit(1, 'pt'),
  # legend.box.margin = margin(t=0, b=0, l=0,r=0),
  legend.spacing = unit(5, 'pt'),
  legend.key.height = unit(5, 'pt'),
  legend.key.width = unit(5, 'pt')
)
ptsize = 1.
