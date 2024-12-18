library(readr)
library(dplyr)
library(ggplot2)

theme_set(theme_bw())

full_df <- read_csv("output/relatives_relatedness.csv")

labs <- full_df |>
  distinct(proband1, proband2, proband_region1, proband_region2, prm, relationship)
labeldat <- full_df |>
  group_by(proband1, proband2, prm) |>
  summarise(ypos = max(grm) + .25) |>
  inner_join(labs)


full_df |>
  filter(type == "branch_recap") |>
  group_by('proband1', 'proband2') |>
  ggplot(aes(x=prm,
             y=grm, 
             group=interaction(prm, proband1,proband2), 
             colour=relationship)) +
  geom_boxplot(width=0.01, position = position_dodge(width=0)) +
#  geom_text(data = labeldat, aes(label = relationship, y = ypos), 
#            position = position_dodge(width = .75), 
#            show.legend = FALSE) +
  facet_grid(rows = vars(proband_region1))
             
ggsave("plots/branch_recap_sim_boxplot.png", height = 8, width = 7, dpi=500)

full_df |>
  filter(type == "site_recap") |>
  group_by('proband1', 'proband2') |>
  ggplot(aes(x=prm,
             y=grm, 
             group=interaction(prm, proband1,proband2), 
             colour=proband_region1)) +
  geom_boxplot(width=0.01, position = position_dodge(width=0)) +
  facet_grid(cols = vars(generation), scales = "free_x")
ggsave("plots/site_recap_sim_boxplot.png")
