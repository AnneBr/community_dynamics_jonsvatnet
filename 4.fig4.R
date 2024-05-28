# Figure 4

# Required packages
library(tidyverse)
library(patchwork)

fig4_data <- tibble(case = c("Kilvatn", "Lille Jonsvatn", "Store Jonsvatn", "Total"),
                     mu = c(2.63, 3.01, 2.48, 2.72),                # Mean log abundance - from data analysis table
                     sd = sqrt(c(2.5, 3.2, 2.56, 2.95))) %>%        # SD = sqrt(variance) - from data analysis table
  mutate(x = list(seq(-5, 11, length.out = 1000))) %>% 
  unnest(x) %>% 
  mutate(ds = dnorm(x, mean = mu, sd = sd))

fig4 <- fig4_data %>% 
  ggplot(aes(x = x, y = ds)) + 
  geom_line(aes(colour = case), linewidth = 1) + 
  geom_line(aes(x = x, y = ds, colour = case), 
            data = fig4_data[fig4_data$case == "Total",], 
            linewidth = 2.5) +
  theme_bw() +
  ylab("Probability density") +
  xlab("Log abundance of individuals per species") +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
  scale_x_continuous(limits = c(-2, 8), breaks = seq(-2, 8, 1)) +
  scale_color_discrete(name = "Total variance:",
                       breaks = c("Kilvatn", "Lille Jonsvatn", "Store Jonsvatn", "Total"),
                       labels = c("Kilvatn", "Lille Jonsvatn", "Store Jonsvatn", "Total")) +
  theme(axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.position = c(0.80, 0.75),
        plot.margin = margin(-5.5, 5.5, 5.5, 5.5, "pt"))

fig4

