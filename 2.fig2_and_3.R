library(tidyverse)
data_agg <- read_rds("data_agg_depth.rds")

## This script contains code for fig. 2 and 3, as well as some descriptive statistics/simple tests
    # descriptive statistics
range(data_agg$individualCount)
mean(data_agg$individualCount)
mean(data_agg$individualCount[data_agg$phylum == "Rotifera"])
mean(data_agg$individualCount[data_agg$phylum == "Arthropoda"])
mean(log(data_agg$individualCount[data_agg$phylum == "Arthropoda"]))
mean(log(data_agg$individualCount[data_agg$phylum == "Rotifera"]))


## FIGURE 2 - species accumulation curve
### Species accumulation curve 
years <- sort(unique(data_agg$year)) 
spsummed <- c()
spaccum <- as.data.frame(years, spsummed)
for (i in 1994:2019) {
  spaccum$spsummed[spaccum$years == i] <- length(unique(data_agg$species[data_agg$year %in% 1994:i]))
}
sp_year <- data_agg %>% 
  group_by(year) %>% 
  summarise(speciesCount = length(unique(species)),
            mean_logab = mean(log(individualCount)),
            mean_ab = mean(individualCount))

colnames(spaccum) <- c("year", "spsummed")
sp_year <- merge(spaccum, sp_year, by = "year")

plota <- ggplot(data = sp_year, aes(year, spsummed)) +
  geom_point(aes(), col = "blue", size = 3) +  # Increase point size
  geom_line(aes(), col = "blue", linewidth = 1.2) +  # Increase line size
  geom_point(aes(year, speciesCount), size = 3) +  # Increase point size
  geom_line(aes(year, speciesCount), linetype = "dashed", linewidth = 1.2) +  # Increase line size
  ylab("No. species") +
  xlab("Year") +
  scale_y_continuous(limits = c(5, 20), breaks = seq(5, 20, 5)) +
  theme_classic() +  # Remove background guiding lines
  theme(text = element_text(size = 20),  # Increase text size
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),# Increase axis title size
        legend.text = element_text(size = 20),  # Increase legend text size
        legend.title = element_text(size = 20),
        plot.margin = unit(c(1, 2, 1.5, 1), "lines"))  

plota

### Does the cumulative number of species recorded in a year depend on the mean total abundance?
# get model estimates
# Fit a linear regression model
model1 <- lm(speciesCount ~ mean_logab, data = sp_year)

# Print the summary of the model
summary(model1)
## conclusion: it does not

# Does the mean abundance change over time, in a linear way?
model2 <- lm(mean_ab ~ year, data = sp_year)

# Print the summary of the model
summary(model2)


## FIGURE 3 - mean abundance per year, depth, group

# Including depths
# Initialize an empty data frame to store the results
depth_data <- data.frame()
depths <- c(17.5, 2.5, 7.5, 12.5)

# Loop through each depth and calculate means and SDs for each phylum
for (depth in depths) {
  for (i in 1994:2019) {
    rotif_data <- data_agg$individualCount[data_agg$phylum == "Rotifera" & data_agg$year == i & data_agg$depth == depth]
    arthropoda_data <- data_agg$individualCount[data_agg$phylum == "Arthropoda" & data_agg$year == i & data_agg$depth == depth]
    
    depth_data <- rbind(depth_data, data.frame(
      year = i,
      depth = depth,
      mrotif = mean(rotif_data),
      sdrotif = sd(rotif_data),
      mart = mean(arthropoda_data),
      sdart = sd(arthropoda_data)
    ))
  }
}

# Now, depth_data contains the mean and SD values for each phylum at each depth.
# You can use this data to create the plot.

# 8-panel plot
  # years only on lowest row

# Set up the panel layout
par(mfrow = c(4, 2))
par(mar = c(1, 3, 1, 1))  # Adjust the values as needed
par(oma = c(1.5, 3.5, 1.5, 1.5))

guideline_art <- c(10, 20)  # Adjust these heights as needed
guideline_rot <- c(250, 500, 750)


# Create plots for each depth
# Subset the data for the current depth
subset_data_arthropoda <- subset(depth_data, depth_data$depth == 2.5)
subset_data_rotifera <- subset(depth_data, depth_data$depth == 2.5)

# Create the plot for mart in the left column
plot(subset_data_arthropoda$year, subset_data_arthropoda$mart, type = "b",
     xlab = " ", ylab = " ", xlim = range(depth_data$year),
     ylim = c(0,30),
     xaxt='n', cex.axis = 1.2)
abline(h = guideline_art, col = "gray", lty = 2)
# Add tick marks to the x-axis without numbers
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020), labels = FALSE, tick = TRUE)
mtext("0-5 m", side = 2, line = 3, cex = 1)
mtext("Arthropoda", side =3, line = 1)

# Create the plot for mrotif in the right column
plot(subset_data_rotifera$year, subset_data_rotifera$mrotif, type = "b",
     xlab = " ", ylab = " ", xlim = range(depth_data$year),
     ylim = c(0,1000),
     xaxt = 'n', cex.axis = 1.2)
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020), labels = FALSE, tick = TRUE)
abline(h = guideline_rot, col = "gray", lty = 2)
mtext("Rotifera", side =3, line = 1)

# Subset the data for the current depth
subset_data_arthropoda <- subset(depth_data, depth_data$depth == 7.5)
subset_data_rotifera <- subset(depth_data, depth_data$depth == 7.5)

# Create the plot for mart in the left column
plot(subset_data_arthropoda$year, subset_data_arthropoda$mart, type = "b",
     xlab = " ", ylab = " ", xlim = range(depth_data$year),
     ylim = c(0,30),
     xaxt = 'n', cex.axis = 1.2)
abline(h = guideline_art, col = "gray", lty = 2)
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020), labels = FALSE, tick = TRUE)
mtext("5-10 m", side = 2, line = 3, cex = 1)

# Create the plot for mrotif in the right column
plot(subset_data_rotifera$year, subset_data_rotifera$mrotif, type = "b",
     xlab = " ", ylab = " ", xlim = range(depth_data$year),
     ylim = c(0,1000),
     xaxt = 'n', cex.axis = 1.2)
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020), labels = FALSE, tick = TRUE)
abline(h = guideline_rot, col = "gray", lty = 2)


# Subset the data for the current depth
subset_data_arthropoda <- subset(depth_data, depth_data$depth == 12.5)
subset_data_rotifera <- subset(depth_data, depth_data$depth == 12.5)

# Create the plot for mart in the left column
plot(subset_data_arthropoda$year, subset_data_arthropoda$mart, type = "b",
     xlab = " ", ylab = " ", xlim = range(depth_data$year),
     ylim = c(0,30),
     xaxt = 'n', cex.axis = 1.2)
abline(h = guideline_art, col = "gray", lty = 2)
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020), labels = FALSE, tick = TRUE)
mtext("10-15 m", side = 2, line = 3, cex = 1)
mtext("Mean Abundance", side = 2, line = 5, at = 30)

# Create the plot for mrotif in the right column
plot(subset_data_rotifera$year, subset_data_rotifera$mrotif, type = "b",
     xlab = " ", ylab = " ", xlim = range(depth_data$year),
     ylim = c(0,1000),
     xaxt = 'n', cex.axis = 1.2)
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020), labels = FALSE, tick = TRUE)
abline(h = guideline_rot, col = "gray", lty = 2)


  # Subset the data for the current depth
  subset_data_arthropoda <- subset(depth_data, depth_data$depth == 17.5)
  subset_data_rotifera <- subset(depth_data, depth_data$depth == 17.5)
  
  # Create the plot for mart in the left column
  plot(subset_data_arthropoda$year, subset_data_arthropoda$mart, type = "b",
       xlab = "Arthropoda", ylab = " ", xlim = range(depth_data$year),
       ylim = c(0,30), cex.axis = 1.2)
  abline(h = guideline_art, col = "gray", lty = 2)
  mtext("15-20 m", side = 2, line = 3, cex = 1)
  
  # Create the plot for mrotif in the right column
  plot(subset_data_rotifera$year, subset_data_rotifera$mrotif, type = "b",
       xlab = "Rotifera", ylab = " ", xlim = range(depth_data$year),
       ylim = c(0,1000), cex.axis = 1.2)
  abline(h = guideline_rot, col = "gray", lty = 2)
 
 
# Reset the panel layout
par(mfrow = c(1, 1))

# Relationship between mean abundance and depth for the different groups
model3 <- lm(mrotif ~ depth, data = depth_data)

# Print the summary of the model
summary(model3)

model4 <- lm(mart ~ depth, data = depth_data)

# Print the summary of the model
summary(model4)

