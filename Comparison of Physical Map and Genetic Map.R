

###### Comparison of physical map and genetic map ########

# Read the merged data frame from a text file
merged_df <- read.delim("phys-gen-compare-Provena_A_GS7_B.txt", head = T)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create a flag for mismatched chromosomes in new column
merged_df <- merged_df %>%
  mutate(chrom_mismatch = ifelse(phys_chrom != gen_chrom, "Mismatch", "Match"))


# Plot correlation between physical map and genetic map

ggplot(merged_df, aes(x = gen_pos, y = phys_pos, color = chrom_mismatch)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  facet_wrap(~ phys_chrom, scales = "free_y") +  # Adjust scales to "free_y" for swapped axis
  labs(
    x = expression(bold("Genetic Map Position (cM)")),
    y = expression(bold("Physical Map Position (Mb)")),
    title = expression(bold("Correlation Between Genetic Map and Physical Map by Chromosome"))
  ) +
  scale_y_continuous(labels = function(y) y / 1e6, expand = c(0, 0)) +  # Scale for y-axis
  scale_color_manual(values = c("Match" = "blue", "Mismatch" = "red")) +
  theme_minimal(base_size = 15) + # Larger base font size
  theme(
    panel.background = element_rect(fill = "lightgrey", color = NA),  # Lighter grey background
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14),  # Larger chromosome labels
    legend.position = c(1, 0),  # Legend at the bottom right corner
    legend.justification = c("right", "bottom"),  # Justify to the right and bottom
    legend.box.margin = margin(10, 10, 10, 10),  # Margin around the legend box
    legend.margin = margin(0, 0, 0, 0),  # Margin for the legend (top, right, bottom, left)
    legend.title = element_text(face = "bold"),  # Use 'face' argument for bold text
    legend.text = element_text(size = 12),  # Increase legend text size
    axis.ticks = element_line(size = 0.5),  # Control the size of axis ticks
    axis.text = element_text(size = 12)  # Adjust axis text size
  )


# Create a new column to highlight mismatched markers between 4A and 4D

filtered_df$highlight <- ifelse(
  (filtered_df$phys_chrom == "4A" & filtered_df$gen_chrom == "4D") |
    (filtered_df$phys_chrom == "4D" & filtered_df$gen_chrom == "4A"), 
  "highlight", 
  "normal"
)

# Create a new column to classify markers based on their physical and genetic chromosome

filtered_df$chrom_class <- with(filtered_df, ifelse(phys_chrom == "4A" & gen_chrom == "4A", "4A_to_4A",
                                                    ifelse(phys_chrom == "4A" & gen_chrom == "4D", "4A_to_4D",
                                                           ifelse(phys_chrom == "4D" & gen_chrom == "4A", "4D_to_4A", "4D_to_4D"))))


# Plotting correlation of physical and genetic maps for chromosome 4A and 4D code with filtered data

ggplot(filtered_df, aes(x = gen_pos, y = phys_pos, color = chrom_class)) +
  geom_point(alpha = 0.6, size = 3) +  # Increase point size to 3 (adjust as needed)
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  facet_wrap(~ phys_chrom, scales = "free_y") +  # Adjust scales to "free_y" for swapped axis
  labs(
    x = expression(bold("Genetic Map Position (cM)")),
    y = expression(bold("Physical Map Position (Mb)")),
    title = expression(bold("Correlation Between Genetic Map and Physical Map by Chromosome"))
  ) +
  scale_y_continuous(labels = function(y) y / 1e6, expand = c(0, 0)) +  # Scale for y-axis
  scale_color_manual(values = c("4A_to_4A" = "cyan3", "4A_to_4D" = "purple", 
                                "4D_to_4A" = "cyan3", "4D_to_4D" = "purple")) +  # Different colors for each category
  theme_minimal(base_size = 15) + # Larger base font size
  theme(
    panel.background = element_rect(fill = "lightgrey", color = NA),  # Lighter grey background
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14),  # Larger chromosome labels
    legend.position = "right",  # Legend to the right
    legend.justification = "top",  # Justify legend to the top
    legend.box.margin = margin(10, 10, 10, 10),  # Margin around the legend box
    legend.margin = margin(0, 0, 0, 0),  # Margin for the legend (top, right, bottom, left)
    legend.title = element_text(face = "bold"),  # Use 'face' argument for bold text
    legend.text = element_text(size = 12)  # Increase legend text size
  )
