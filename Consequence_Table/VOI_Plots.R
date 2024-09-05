library(tidyverse)

setwd("C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_SummerFallX2_VOI/Consequence_Table")

cons_table <- read.csv(file.path("Data","ConsequenceTable_2024-08-21.csv"))

# Standardize score based on local scale

dsm <- cons_table %>% filter(Objective=="DeltaSmelt") %>%
  mutate(std_score=(Score-min(Score))/(max(Score)-min(Score)))

water <- cons_table %>% filter(Objective=="WaterCost") %>%
  mutate(std_score=(max(Score)-Score)/(max(Score)-min(Score)))


# Hypothesis weights
cons_table_std <- bind_rows(dsm, water) %>% mutate(hypo_weight= case_when(Hypothesis == "H1" ~ 0.575,
                                                                          Hypothesis == "H2" ~ 0.175,
                                                                          Hypothesis == "H3" ~ 0.175,
                                                                          Hypothesis == "H4" ~ 0.075)) %>%
  mutate(score_hyp = hypo_weight*std_score) 

# Create dataset necessary for the line plot
line_plot_data <- cons_table_std %>% group_by(Alternatives,Objective) %>%
  summarise(comp_score = sum(score_hyp)) %>% mutate(fish_weight = case_when(Objective == "DeltaSmelt" ~ 1.0,
                                                                            Objective == "WaterCost" ~ 0))
# Define custom colors
custom_colors <- c("Alt F80" = "black","Alt F74" = "#E69F00", "Alt S74" = "#0072B2", "Alt S74F80"= "#D55E00","Alt NoX2" = "black","Alt NoFlow" = "black")
custom_line <- c("Alt F80" = "solid","Alt F74" = "solid", "Alt S74" = "longdash", "Alt S74F80"= "longdash","Alt NoX2" = "dotted","Alt NoFlow" = "dotted")

# Convert the 'Category' column to a factor with a specific order
line_plot_data$Alternatives <- factor(line_plot_data$Alternatives, levels = c("Alt F80", "Alt F74", "Alt S74",
                                                                                "Alt S74F80","Alt NoX2"))
# Change name for AFS
#levels(line_plot_data$Alternatives)[levels(line_plot_data$Alternatives) == "Alt NoX2"] <- "Alt NoFlow"


# Line plot
plot_scoreline <- ggplot(data=line_plot_data, aes(x=fish_weight, y=comp_score, color=Alternatives,linetype=Alternatives)) +
  geom_line(linewidth= 1.2) +
  theme_minimal()+
  labs(title = NULL,
       x = "Delta Smelt objective weight",
       y = "Composite score (objective-weighted linear value function)") +
  theme(axis.text = element_text(size = 14),  # Increase tick mark font size
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  legend.text=element_text(size=14),
  legend.title=element_text(size=14),
  axis.title.x = element_text(size=14),
  axis.title.y = element_text(size=14),
  legend.key.size = unit(1.2, "cm"))+
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = custom_line)
plot_scoreline

# Export plot
tiff(filename=file.path("Output","Figure_LinePlot_2Objectives.tiff"),
     type="cairo",
     units="in", 
     width=10, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_scoreline
dev.off()

######
# Expected Value of Perfect Information calculation


#####
# Initialize the results data frame
EVPI_results <- data.frame(fish_weight = numeric(), EVPI = numeric(), EV_uncertainty = numeric())

# Assuming cons_table_std is already defined and contains the necessary columns

# Define a function to calculate EVPI for a given fish weight
calculate_EVPI <- function(i) {
  cons_table_reconfig <- cons_table_std %>% 
    dplyr::select(Alternatives, Hypothesis, Objective, std_score, hypo_weight) %>% 
    spread(Objective, std_score) %>% 
    mutate(composite_score = (DeltaSmelt * i) + (WaterCost * (1 - i)))
  
  certainty_calc <- cons_table_reconfig %>% 
    group_by(Hypothesis) %>% 
    summarise(composite_score = max(composite_score), hypo_weight = mean(hypo_weight)) %>%
    mutate(hypothesis_score = composite_score * hypo_weight)
  
  EV_under_certainty <- sum(certainty_calc$hypothesis_score)
  
  uncertainty_calc <- cons_table_reconfig %>% 
    mutate(composite_score_hypo = hypo_weight * composite_score) %>% 
    group_by(Alternatives) %>% 
    summarise(composite_score = sum(composite_score_hypo))
  
  EV_under_uncertainty <- max(uncertainty_calc$composite_score)
  
  EV_PI <- EV_under_certainty - EV_under_uncertainty
  
  # Create a new row with the fish weight and EVPI results
  new_row <- data.frame(fish_weight = i,
                        EVPI = EV_PI, 
                        EV_uncertainty = EV_under_uncertainty)
  
  return(new_row)
}

# Use lapply to calculate EVPI for each fish weight
EVPI_results_list <- lapply(seq(0, 1, 0.01), calculate_EVPI)

# Combine the results into a single data frame
EVPI_results <- do.call(rbind, EVPI_results_list)

# Add percentage of value per Brian Healy's figure
EVPI_results <- EVPI_results %>% mutate(percent_EVPI = EVPI/EV_uncertainty*100)

#####
# Create plot for EVPI
plot_EVPI <- ggplot(EVPI_results, aes(x = fish_weight, y = percent_EVPI)) +
  geom_line(color = "navy", linewidth = 1.2) +  # Line color and thickness
  labs(title = NULL,
       x = "Delta Smelt objective weight",
       y = "Expected value of perfect information (%)") +
  theme_classic() +                          # Classic theme
  theme(axis.text = element_text(size = 14),  # Increase tick mark font size
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
plot_EVPI

# Export plot
tiff(filename=file.path("Output","Figure_EVPI_2Objectives.tiff"),
     type="cairo",
     units="in", 
     width=8, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_EVPI
dev.off()

        
######
# Partial Expected Value of Perfect Information calculation

# Partial EVPI for resolving the food hypothesis

# Set up the data frame
EVPI_results_partial_food <- data.frame(fish_weight = numeric(),
                                        EVPI = numeric(), 
                                        Hypothesis = factor())

# Define a function to calculate partial EVPI for a given fish weight
calculate_partial_EVPI <- function(i) {
  cons_table_reconfig <- cons_table_std %>% 
    dplyr::select(Alternatives, Hypothesis, Objective, std_score, hypo_weight) %>% 
    spread(Objective, std_score) %>% 
    mutate(composite_score = (DeltaSmelt * i) + (WaterCost * (1 - i))) %>%
    # Create a new column that defines the food hypothesis
    mutate(Hypothesis_food = case_when(Hypothesis == "H1" ~ "Food",
                                       Hypothesis == "H2" ~ "Food",
                                       Hypothesis == "H3" ~ "No Food",
                                       Hypothesis == "H4" ~ "No Food"))
  
  partial_certainty_calc <- cons_table_reconfig %>% 
    group_by(Hypothesis_food, Alternatives) %>% 
    mutate(composite_score_hypo = hypo_weight * composite_score) %>%
    summarise(composite_score = sum(composite_score_hypo), 
              hypo_weight = sum(hypo_weight), .groups = 'drop') %>%
    # Remember to divide everything by the collapsed/combined hypothesis weights
    mutate(composite_score = composite_score / hypo_weight) %>% 
    ungroup() %>%
    group_by(Hypothesis_food) %>% 
    summarise(hypo_weight = mean(hypo_weight), 
              best_score = max(composite_score), .groups = 'drop') %>%
    mutate(score_for_EVPI = hypo_weight * best_score)
  
  EV_under_partial_certainty <- sum(partial_certainty_calc$score_for_EVPI)
  
  uncertainty_calc <- cons_table_reconfig %>% 
    mutate(composite_score_hypo = hypo_weight * composite_score) %>% 
    group_by(Alternatives) %>% 
    summarise(composite_score = sum(composite_score_hypo), .groups = 'drop')
  
  EV_under_uncertainty <- max(uncertainty_calc$composite_score)
  
  EVPI_partial_food <- EV_under_partial_certainty - EV_under_uncertainty
  
  # Create a new row with the fish weight and EVPI results
  new_row <- data.frame(fish_weight = i,
                        EVPI = EVPI_partial_food,
                        Hypothesis = "Food")
  
  return(new_row)
}

# Use lapply to calculate the partial EVPI for each fish weight
EVPI_results_list <- lapply(seq(0, 1, 0.01), calculate_partial_EVPI)

# Combine the results into a single data frame
EVPI_results_partial_food <- do.call(rbind, EVPI_results_list)


# Set up the data frame
EVPI_results_partial_movement <- data.frame(fish_weight = numeric(),
                                            EVPI = numeric(), 
                                            Hypothesis = factor())

# Define a function to calculate partial EVPI for a given fish weight
calculate_partial_EVPI_movement <- function(i) {
  cons_table_reconfig <- cons_table_std %>% 
    dplyr::select(Alternatives, Hypothesis, Objective, std_score, hypo_weight) %>% 
    spread(Objective, std_score) %>% 
    mutate(composite_score = (DeltaSmelt * i) + (WaterCost * (1 - i))) %>%
    # Create a new column that defines the movement hypothesis
    mutate(Hypothesis_movement = case_when(Hypothesis == "H1" ~ "Movement",
                                           Hypothesis == "H2" ~ "No Movement",
                                           Hypothesis == "H3" ~ "Movement",
                                           Hypothesis == "H4" ~ "No Movement"))
  
  partial_certainty_calc <- cons_table_reconfig %>% 
    group_by(Hypothesis_movement, Alternatives) %>% 
    mutate(composite_score_hypo = hypo_weight * composite_score) %>%
    summarise(composite_score = sum(composite_score_hypo), 
              hypo_weight = sum(hypo_weight), .groups = 'drop') %>%
    # Remember to divide everything by the collapsed/combined hypothesis weights
    mutate(composite_score = composite_score / hypo_weight) %>% 
    ungroup() %>%
    group_by(Hypothesis_movement) %>% 
    summarise(hypo_weight = mean(hypo_weight), 
              best_score = max(composite_score), .groups = 'drop') %>%
    mutate(score_for_EVPI = hypo_weight * best_score)
  
  EV_under_partial_certainty <- sum(partial_certainty_calc$score_for_EVPI)
  
  uncertainty_calc <- cons_table_reconfig %>% 
    mutate(composite_score_hypo = hypo_weight * composite_score) %>% 
    group_by(Alternatives) %>% 
    summarise(composite_score = sum(composite_score_hypo), .groups = 'drop')
  
  EV_under_uncertainty <- max(uncertainty_calc$composite_score)
  
  EVPI_partial_movement <- EV_under_partial_certainty - EV_under_uncertainty
  
  # Create a new row with the fish weight and EVPI results
  new_row <- data.frame(fish_weight = i,
                        EVPI = EVPI_partial_movement,
                        Hypothesis = "Movement")
  
  return(new_row)
}

# Use lapply to calculate the partial EVPI for each fish weight
EVPI_results_list_movement <- lapply(seq(0, 1, 0.01), calculate_partial_EVPI_movement)

# Combine the results into a single data frame
EVPI_results_partial_movement <- do.call(rbind, EVPI_results_list_movement)

# Combine the two datasets together
EVPI_results_combined <- bind_rows(EVPI_results_partial_food,EVPI_results_partial_movement)

EVPI_results_combined <- EVPI_results_combined %>% left_join(EVPI_results %>% select(fish_weight,EV_uncertainty)) %>%
  mutate(percent_EVPI = EVPI/EV_uncertainty*100)

# Define custom colors
custom_colors <- c("Food" = "blue2", "Movement" = "yellow4")

# Stacked area chart
plot_partial_EVPI <- ggplot(EVPI_results_combined, aes(x=fish_weight, y=percent_EVPI, fill=Hypothesis)) + 
  geom_area(alpha=0.6, position='stack') +  
  labs(title = NULL,
       x = "Delta Smelt objective weight",
       y = "Expected value of partial information (%)") +
  theme_classic() +                          # Classic theme 
  theme(axis.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))+
  theme(legend.title = element_blank(),      # Remove legend title
        legend.position = c(0.2, 0.5),            # Position the legend
  legend.background = element_rect(color = "black"), # Legend background
  legend.key = element_rect(fill = "lightgray")) + # Legend key color
  scale_fill_manual(values = custom_colors)     # Use a color palette

plot_partial_EVPI

# Export plot
tiff(filename=file.path("Output","Figure_EVPI_partial_2Objectives.tiff"),
     type="cairo",
     units="in", 
     width=8, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_partial_EVPI
dev.off()

#######
# Alternate way of showing data instead of consequence table
# Using bar chart

# Define custom colors
custom_colors_alt <- c("Alt F80" = "#000000", "Alt F74" = "#E69F00",
                       "Alt S74" = "yellow4" , "Alt S74F80" = "#56B4E9","Alt NoX2"= "#999999")

# Grab data for the water cost plot
data_plot_water <- cons_table %>% filter(Objective=="WaterCost"&Hypothesis=="H1") %>%
  select(-Hypothesis)

# Convert the 'Category' column to a factor with a specific order
data_plot_water$Alternatives <- factor(data_plot_water$Alternatives, levels = c("Alt F80", "Alt F74", "Alt S74",
                                                                                "Alt S74F80","Alt NoX2"))

# Change name for AFS
levels(data_plot_water$Alternatives)[levels(data_plot_water$Alternatives) == "Alt NoX2"] <- "Alt NoFlow"

# Water cost plot
plot_bar_water <- ggplot(data_plot_water, aes(x=Alternatives, y=Score, fill=Alternatives)) +
  geom_bar(stat = "identity") +
  labs(title = "Water Cost Objective",
       x = "Alternative",
       y = "Thousand Acre Feet") + 
  scale_fill_manual(values = custom_colors_alt,guide="none")  +   # Use a color palette
  theme_bw() +                          # Classic theme 
  theme(axis.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
#plot_bar_water

# Export plot
tiff(filename=file.path("Output","Figure_BarPlot_WaterObjective.tiff"),
     type="cairo",
     units="in", 
     width=8, #10*1, 
     height=5, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_bar_water
dev.off()

# Grab data for the delta smelt plot
data_plot_dsm <- dsm

# Convert the 'Category' column to a factor with a specific order
data_plot_dsm$Alternatives <- factor(data_plot_dsm$Alternatives, levels = c("Alt F80", "Alt F74", "Alt S74",
                                                                                "Alt S74F80","Alt NoX2"))
# Change name for AFS
levels(data_plot_dsm$Alternatives)[levels(data_plot_dsm$Alternatives) == "Alt NoX2"] <- "Alt NoFlow"

# Delta Smelt obj plot
plot_bar_dsm <- ggplot(data_plot_dsm, aes(x=Alternatives, y=Score, fill=Alternatives)) +
  geom_bar(stat = "identity") +
  labs(title = "Delta Smelt Objective",
       x = "Alternative",
       y = "Lambda") + 
  facet_grid(cols = vars(Hypothesis)) +
  scale_fill_manual(values = custom_colors_alt,guide="none")  +   # Use a color palette
  theme_bw() +                          # Classic theme 
  coord_cartesian(ylim=c(0.75,1)) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size = 14))
plot_bar_dsm


# Export plot
tiff(filename=file.path("Output","Figure_BarPlot_DeltaSmeltObjective.tiff"),
     type="cairo",
     units="in", 
     width=11, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_bar_dsm
dev.off()













#################################
#Excess code DO NOT USE
# VOI
fish_weight <- 0.87

cons_table_reconfig <- cons_table_std %>% dplyr::select(Alternatives,Hypothesis,Objective,std_score,hypo_weight) %>% 
  spread(Objective,std_score) %>% mutate(composite_score = (DeltaSmelt*fish_weight)+(WaterCost*(1-fish_weight)) )

certainty_calc <- cons_table_reconfig %>% group_by(Hypothesis) %>% summarise(composite_score=max(composite_score),hypo_weight = mean(hypo_weight)) %>%
  mutate(hypothesis_score = composite_score*hypo_weight)

EV_under_certainty <- sum(certainty_calc$hypothesis_score)

uncertainty_calc <- cons_table_reconfig %>% mutate(composite_score_hypo = hypo_weight*composite_score) %>% group_by(Alternatives) %>% 
  summarise(composite_score = sum(composite_score_hypo))

EV_under_uncertainty <- max(uncertainty_calc$composite_score)

EVPI <- EV_under_certainty - EV_under_uncertainty

# Old code using for loop before switching to lapply
# Set up the data frame
EVPI_results <- data.frame(fish_weight = numeric(),
                           EVPI = numeric(), EV_uncertainty=numeric())
# Loop function across objective weights
for (i in seq(0, 1, 0.01)) {
  cons_table_reconfig <- cons_table_std %>% dplyr::select(Alternatives,Hypothesis,Objective,std_score,hypo_weight) %>% 
    spread(Objective,std_score) %>% mutate(composite_score = (DeltaSmelt*i)+(WaterCost*(1-i)) )
  
  certainty_calc <- cons_table_reconfig %>% group_by(Hypothesis) %>% summarise(composite_score=max(composite_score),hypo_weight = mean(hypo_weight)) %>%
    mutate(hypothesis_score = composite_score*hypo_weight)
  
  EV_under_certainty <- sum(certainty_calc$hypothesis_score)
  
  uncertainty_calc <- cons_table_reconfig %>% mutate(composite_score_hypo = hypo_weight*composite_score) %>% group_by(Alternatives) %>% 
    summarise(composite_score = sum(composite_score_hypo))
  
  EV_under_uncertainty <- max(uncertainty_calc$composite_score)
  
  EV_PI <- EV_under_certainty - EV_under_uncertainty
  
  # Create a new row with the iteration number, value, and calculation description
  new_row <- data.frame(fish_weight = i,
                        EVPI = EV_PI, EV_uncertainty=EV_under_uncertainty)
  
  # Append the new row to the results data frame
  EVPI_results <- rbind(EVPI_results, new_row)
}

# Add percentage of value per Brian Healy's figure
EVPI_results <- EVPI_results %>% mutate(percent_EVPI = EVPI/EV_uncertainty*100)

#####
# Food only
fish_weight <- 0.2

cons_table_reconfig <- cons_table_std %>% dplyr::select(Alternatives,Hypothesis,Objective,std_score,hypo_weight) %>% 
  spread(Objective,std_score) %>% mutate(composite_score = (DeltaSmelt*fish_weight)+(WaterCost*(1-fish_weight)) ) %>%
  # Create a new column that defines the food hypothesis
  mutate(Hypothesis_food = case_when(Hypothesis == "H1" ~ "Food",
                                     Hypothesis == "H2" ~ "Food",
                                     Hypothesis == "H3" ~ "No Food",
                                     Hypothesis == "H4" ~ "No Food"))

partial_certainty_calc <- cons_table_reconfig %>% group_by(Hypothesis_food,Alternatives) %>% mutate(composite_score_hypo = hypo_weight*composite_score) %>%
  summarise(composite_score = sum(composite_score_hypo),hypo_weight=sum(hypo_weight)) %>% mutate(composite_score=composite_score/hypo_weight) %>% ungroup() %>%
  group_by(Hypothesis_food) %>% summarise(hypo_weight=mean(hypo_weight),best_score=max(composite_score)) %>% 
  mutate(score_for_EVPI = hypo_weight*best_score)

EV_under_partial_certainty <- sum(partial_certainty_calc$score_for_EVPI)

uncertainty_calc <- cons_table_reconfig %>% mutate(composite_score_hypo = hypo_weight*composite_score) %>% group_by(Alternatives) %>% 
  summarise(composite_score = sum(composite_score_hypo))

EV_under_uncertainty <- max(uncertainty_calc$composite_score)

EVPI_partial_food <- EV_under_partial_certainty - EV_under_uncertainty
EV_under_partial_certainty - EV_under_uncertainty


# Partial EVPI for resolving the movement hypothesis
# Set up the data frame
EVPI_results_partial_movement <- data.frame(fish_weight = numeric(),
                                            EVPI = numeric(), Hypothesis = factor())

# Loop function across objective weights
for (i in seq(0, 1, 0.01)) {
  cons_table_reconfig <- cons_table_std %>% dplyr::select(Alternatives,Hypothesis,Objective,std_score,hypo_weight) %>% 
    spread(Objective,std_score) %>% mutate(composite_score = (DeltaSmelt*i)+(WaterCost*(1-i)) ) %>%
    # Create a new column that defines the food hypothesis
    mutate(Hypothesis_movement = case_when(Hypothesis == "H1" ~ "Movement",
                                           Hypothesis == "H2" ~ "No Movement",
                                           Hypothesis == "H3" ~ "Movement",
                                           Hypothesis == "H4" ~ "No Movement"))
  
  partial_certainty_calc <- cons_table_reconfig %>% group_by(Hypothesis_movement,Alternatives) %>% mutate(composite_score_hypo = hypo_weight*composite_score) %>%
    summarise(composite_score = sum(composite_score_hypo),hypo_weight=sum(hypo_weight)) %>% 
    # Remember to divide everything by the collapsed/combined hypothesis weights
    mutate(composite_score=composite_score/hypo_weight) %>% ungroup() %>%
    group_by(Hypothesis_movement) %>% summarise(hypo_weight=mean(hypo_weight),best_score=max(composite_score)) %>% 
    mutate(score_for_EVPI = hypo_weight*best_score)
  
  EV_under_partial_certainty <- sum(partial_certainty_calc$score_for_EVPI)
  
  uncertainty_calc <- cons_table_reconfig %>% mutate(composite_score_hypo = hypo_weight*composite_score) %>% group_by(Alternatives) %>% 
    summarise(composite_score = sum(composite_score_hypo))
  
  EV_under_uncertainty <- max(uncertainty_calc$composite_score)
  
  EVPI_partial_movement <- EV_under_partial_certainty - EV_under_uncertainty
  
  # Create a new row with the iteration number, value, and calculation description
  new_row <- data.frame(fish_weight = i,
                        EVPI = EVPI_partial_movement,
                        Hypothesis = "Movement")
  
  # Append the new row to the results data frame
  EVPI_results_partial_movement <- rbind(EVPI_results_partial_movement, new_row)
}



