# 2) Fraction of probands with at least one founder at each depth
# First, for each depth, find how many distinct probands have a founder
fraction_data <- depths_long %>%
  distinct(proband, depth) %>%
  group_by(depth) %>%
  summarise(count = n_distinct(proband), .groups = "drop")

total_probands <- length(list_of_probands)

# Ensure depth = 1 is included, even if not present
fraction_data <- fraction_data %>%
  mutate(fraction = count / total_probands)

fraction_data_stacked <- fraction_data %>%
  mutate(no_founder = 1 - fraction) %>%
  pivot_longer(cols = c("fraction", "no_founder"),
               names_to = "status", values_to = "val") %>%
  mutate(status = factor(status, levels = c("no_founder", "fraction")))

fraction_plot <- ggplot(fraction_data_stacked, aes(x = factor(depth), y = val, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c( "tomato", "royalblue"), 
                    labels = c("No Founder", "Has Founder")) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(
    title = "Fraction of Probands with Founders at Given Depth",
    x = "Depth",
    y = "Fraction of Probands",
    fill = ""
)

# Compute genetic contribution for each founder
# contribution = 2^-depth
genetic_contributions <- depths_long %>%
  mutate(contribution = 2^(-depth))

# Determine the maximum depth
max_depth <- max(genetic_contributions$depth, na.rm = TRUE)

# For each depth d, sum contributions of founders with depth <= d, by proband
cumulative_contribution <- lapply(1:max_depth, function(d) {
  genetic_contributions %>%
    filter(depth <= d) %>%
    group_by(proband) %>%
    summarise(cumulative_contrib = sum(contribution), .groups = "drop") %>%
    mutate(depth = d)
}) %>% bind_rows()

# Average across probands at each depth
avg_contribution <- cumulative_contribution %>%
  group_by(depth) %>%
  summarise(avg_contrib = mean(cumulative_contrib, na.rm = TRUE), .groups = "drop")

# Plot the average cumulative genetic coverage
cumulative_plot <- ggplot(avg_contribution, aes(x = factor(depth), y = avg_contrib)) +
  geom_line(aes(group=1)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Average Cumulative Genetic Coverage Up to Given Depth",
    x = "Depth",
    y = "Average Cumulative Contribution"
  )

# Combine plots using cowplot (three plots: density_plot, fraction_plot, cumulative_plot)
combined_plot <- plot_grid(
  density_plot,
  fraction_plot,
  cumulative_plot,
  nrow = 3, rel_heights = c(1,1,1),
  align = "v", 
  axis = "lr"
)
