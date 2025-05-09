library(tidyverse)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
rm(list=ls())

setwd("C:/Users/bmahardja/OneDrive - DOI/Documents/GitHub/DeltaSmelt_SummerFallX2_VOI/Consequence_Table")

cons_table <- read.csv(file.path("Data","ConsequenceTable_2025-05-09.csv"))


# Standardize score based on local scale

dsm <- cons_table %>% filter(Objective=="DeltaSmelt") %>%
  mutate(std_score=(Score-min(Score))/(max(Score)-min(Score)))

water <- cons_table %>% filter(Objective=="WaterCost") %>%
  mutate(std_score=(max(Score)-Score)/(max(Score)-min(Score)))

# Hypothesis weights
cons_table_std <- bind_rows(dsm, water) %>% mutate(Hypothesis=as.factor(Hypothesis))

# Create dataset necessary for the line plot
line_plot_data <- cons_table_std %>% mutate(fish_weight = case_when(Objective == "DeltaSmelt" ~ 1.0,
                                                                            Objective == "WaterCost" ~ 0))

# Define custom colors
custom_colors <- c("Alt F80" = "black","Alt F74" = "#E69F00", "Alt S74" = "#0072B2", "Alt S74F80"= "#D55E00","Alt NoX2" = "black","Alt NoFlow" = "black")
custom_line <- c("Alt F80" = "solid","Alt F74" = "solid", "Alt S74" = "longdash", "Alt S74F80"= "longdash","Alt NoX2" = "dotted","Alt NoFlow" = "dotted")

# Convert the 'Category' column to a factor with a specific order
line_plot_data$Alternatives <- factor(line_plot_data$Alternatives, levels = c("Alt F80", "Alt F74", "Alt S74",
                                                                              "Alt S74F80","Alt NoX2"))

# Remove the letter 'H' from the 'Hypothesis' column
line_plot_data$Hypothesis <- as.character(line_plot_data$Hypothesis)  # Convert factor to character
line_plot_data$Hypothesis <- gsub("H", "", line_plot_data$Hypothesis)  # Remove letter 'A'
line_plot_data$Hypothesis <- factor(line_plot_data$Hypothesis)          # Convert back to factor


# Split the data frame by the 'cyl' factor
split_data <- split(line_plot_data, line_plot_data$Hypothesis)

# Create a list to store plots
plots <- list()

# Loop through each subset of the data and create a ggplot
for (Hypothesis in names(split_data)) {
  p <- ggplot(split_data[[Hypothesis]], aes(x = fish_weight, y = std_score, color=Alternatives,linetype=Alternatives)) +
    geom_line(linewidth= 1.2) +
    theme_minimal()+
    labs(title = paste("Hypothesis:", Hypothesis),
         x = "Delta Smelt objective weight",
         y = "Utility score (objective-weighted linear value function)") +
    theme(axis.text = element_text(size = 14),  # Increase tick mark font size
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.key.size = unit(1.2, "cm"),
          legend.position="none")+
    scale_color_manual(values = custom_colors) +
    scale_linetype_manual(values = custom_line)
  
  # Store the plot in the list
  plots[[Hypothesis]] <- p
}


# Print the plots
for (Hypothesis in names(plots)) {
  print(plots[[Hypothesis]])
}

# Combine plots using ggarrange
combined_plot <- ggarrange(plotlist=plots,
                           ncol = 2, nrow = 4,
                           common.legend = TRUE, legend = "top")

# Add common x and y labels
combined_plot <- annotate_figure(combined_plot,
                                 bottom = text_grob("Delta Smelt objective weight", size = 14),
                                 left = text_grob("Composite alternative performance score", rot = 90, vjust = 1, size = 14))

# Create a 2x4 panel plot
tiff(filename=file.path("Output","Figure_LinePlot_ByHypothesis.tiff"),
     type="cairo",
     units="in", 
     width=9, #10*1, 
     height=14, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
#grid.arrange(grobs = plots, ncol = 2,
#             left = textGrob("Composite alternative performance score", rot = 90, gp = gpar(fontsize = 14)),
#             bottom = textGrob("Delta Smelt objective weight", gp = gpar(fontsize = 14)))
combined_plot
dev.off()


