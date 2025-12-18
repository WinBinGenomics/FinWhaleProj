# author: Kaden Winspear @ CSUSM -> Eastern Pacific Fin Whale Project.
# Date; April 2025
# Usage: Used to create box plots for genomewide het across individuals in each population. 

# libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Color palette.
colors <- brewer.pal(3, "Dark2") 

# Define the file path to your CSV file (make sure the file path is correct / adjust for either xchrom or autosomes.)
file_path <- "/Users/kadenwinspear/Documents/finwhale_proj/genomewide_diversity/Xchromwide_estimates/all70_genomewide_XChrom_heterozygosity_20250404.csv"

# Reads the CSV - noice
data <- read.csv(file_path, stringsAsFactors = FALSE)

# Converts PopId to a factor w/ levels: ENP, GOC, ESP.
data$PopId <- factor(data$PopId, levels = c("ENP", "GOC", "ESP"))

# Sumstats
summary_stats <- data %>%
  group_by(PopId) %>%
  summarize(
    n = n(),
    mean = mean(GenomeHet, na.rm = TRUE),
    sd = sd(GenomeHet, na.rm = TRUE),
    median = median(GenomeHet, na.rm = TRUE)
  )
print(summary_stats)

# Creates the plot -> lets make this thing look swag -> added plot instructions. again, swag.
p <- ggplot(data, aes(x = PopId, y = GenomeHet, fill = PopId)) +
  # Creates the boxplot without plotting outlier points.
  # alpha = 0.7 sets a slight transparency for the fill.
  # color = "black" makes the box outline and median line black.
  # fatten = 2.5 increases the thickness of the median line.
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", fatten = 2.5) +
  
  # Set axis labels.
  labs(
    x = "Population",
    #y = "Autosomal Heterozygosity" # adjust if autosomes or xchrom
    y = "X Chromosome Heterozygosity"
  ) +
  
  theme_minimal() +
  theme(
    # Removes them ugly grid lines.
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Adds x and y axis lines.
    axis.line = element_line(color = "black", size = 0.5),
    # tick marks
    axis.ticks = element_line(color = "black", size = 0.9),
    # Bolds legend text and title.
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    # Customize axis text and title appearance.
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(size = 13, face = "bold")
  ) +
  
  # color for pops
  scale_fill_manual(values = c("ENP" = colors[1], "GOC" = colors[2], "ESP" = colors[3])) +
  
  # y-axis +labels
  scale_y_continuous(
    breaks = seq(0.00, 0.002, by = 0.0002),
    labels = scales::number_format(accuracy = 0.00005)
  )

# Display the plot.
print(p)


ggsave(filename = "/Users/kadenwinspear/Documents/finwhale_proj/genomewide_diversity/Xchromwide_estimates/Xchrom_Het_boxplot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Performs pairwise Wilcoxon (Mann–Whitney U) tests for GenomeHet between populations.
pairwise_results <- pairwise.wilcox.test(data$GenomeHet, data$PopId, p.adjust.method = "BH")
print(pairwise_results)
