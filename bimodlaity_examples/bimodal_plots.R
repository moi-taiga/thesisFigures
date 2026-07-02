# Load necessary libraries
library(ggplot2)
library(stats)
library(gridExtra)

# Set a seed for reproducibility so you get the exact same plot
set.seed(27)

# Define the number of samples for each distribution
n_samples1 <- 1000  # You can change this to any number
n_samples2 <- 400   # You can change this to any number

# Generate the data
dist1 <- rnorm(n_samples1, mean = 0, sd = 0.3)
dist2 <- rnorm(n_samples2, mean = 2, sd = 0.85)

# Combine both distributions into a single long-format data frame
bimodal_data <- data.frame(
  value = c(dist1, dist2),
  source = c(rep("Not Expressed", n_samples1), 
             rep("Expressed", n_samples2))
)

# Calculate means and SDs for each distribution
stats <- data.frame(
  source = c("Not Expressed", "Expressed"),
  mean = c(mean(dist1), mean(dist2)),
  sd = c(sd(dist1), sd(dist2))
)

# Function to calculate the intersection of two normal distributions
find_intersection <- function(mean1, sd1, mean2, sd2, n1 = 1, n2 = 1) {
  a <- 1 / (2 * sd1^2) - 1 / (2 * sd2^2)
  b <- mean2 / (sd2^2) - mean1 / (sd1^2)
  c <- mean1^2 / (2 * sd1^2) - mean2^2 / (2 * sd2^2) - log((n1 * sd2) / (n2 * sd1))
  
  roots <- polyroot(c(c, b, a))
  Re(roots[abs(Im(roots)) < 1e-6])  # Return only real roots
}

# Calculate the intersection points
intersection_points <- find_intersection(mean(dist1), sd(dist1), mean(dist2), sd(dist2), n1 = n_samples1, n2 = n_samples2)

# Generate theoretical normal curves so the plot visualization matches the calculated root
x_range <- seq(min(bimodal_data$value), max(bimodal_data$value), length.out = 1000)
theo_data <- data.frame(
  value = rep(x_range, 2),
  count = c(dnorm(x_range, mean(dist1), sd(dist1)) * n_samples1,
            dnorm(x_range, mean(dist2), sd(dist2)) * n_samples2),
  source = rep(c("Not Expressed", "Expressed"), each = 1000)
)

# 1. Plot the overall bimodal distribution
p1 <- ggplot(bimodal_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(count)), binwidth = 0.2, 
                 fill = "lightgrey", colour = "black", alpha = 0.7) +
  stat_function(fun = function(x) (n_samples1 * dnorm(x, mean(dist1), sd(dist1)) + n_samples2 * dnorm(x, mean(dist2), sd(dist2))) * 0.2,
                colour = "darkblue", linewidth = 1, n = 1000) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(0, NA))

# 2. Plot the two overlapping distributions separately to see how they combine
p2 <- ggplot(theo_data, aes(x = value, y = count, fill = source)) +
  geom_area(alpha = 0.5, colour = "black", position = "identity") +
  geom_vline(data = stats, aes(xintercept = mean, colour = source), show.legend = FALSE) +
  geom_vline(data = stats, aes(xintercept = mean + sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(data = stats, aes(xintercept = mean - sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(xintercept = intersection_points, colour = "black", alpha = 0.7, show.legend = FALSE) +
  annotate("text", x = intersection_points, y = -50, label = "Root", colour = "black", angle = 90, vjust = -0.5, size = 4.8) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  scale_colour_manual(values = c("skyblue", "salmon")) +
  labs(
    fill = "Distribution",
    colour = "Distribution"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(xlim = c(0, NA))

png("A.png", width = 1000, height = 800)
grid.arrange(p1, p2, ncol = 1, 
             top = grid::textGrob("Gene Expression from Two Overlapping Distributions (A)", gp = grid::gpar(fontsize = 20, fontface = "bold")),
             left = grid::textGrob("Frequency", rot = 90, gp = grid::gpar(fontsize = 20), just = "center", x = unit(0.5, "npc")),
             bottom = grid::textGrob("Expression", gp = grid::gpar(fontsize = 20)),
             padding = unit(2, "line"))
dev.off()

# ------------ Example B --------------------- 
# Where the mass of the second distribution is 10% of the first
B_n_samples1 <- 1000
B_n_samples2 <- 100

# Generate the data
B_dist1 <- rnorm(B_n_samples1, mean = 0, sd = 0.3)
B_dist2 <- rnorm(B_n_samples2, mean = 1.4, sd = 0.35)

# Combine both distributions into a single long-format data frame
B_bimodal_data <- data.frame(
  value = c(B_dist1, B_dist2),
  source = c(rep("Not Expressed", B_n_samples1), 
             rep("Expressed", B_n_samples2))
)

# Calculate means and SDs for each distribution
B_stats <- data.frame(
  source = c("Not Expressed", "Expressed"),
  mean = c(mean(B_dist1), mean(B_dist2)),
  sd = c(sd(B_dist1), sd(B_dist2))
)

# Calculate the intersection points
B_intersection_points <- find_intersection(mean(B_dist1), sd(B_dist1), mean(B_dist2), sd(B_dist2), n1 = B_n_samples1, n2 = B_n_samples2)

# Generate theoretical normal curves for B
B_x_range <- seq(min(B_bimodal_data$value), max(B_bimodal_data$value), length.out = 1000)
B_theo_data <- data.frame(
  value = rep(B_x_range, 2),
  count = c(dnorm(B_x_range, mean(B_dist1), sd(B_dist1)) * B_n_samples1,
            dnorm(B_x_range, mean(B_dist2), sd(B_dist2)) * B_n_samples2),
  source = rep(c("Not Expressed", "Expressed"), each = 1000)
)

# 3. Plot the overall bimodal distribution
B_p1 <- ggplot(B_bimodal_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(count)), binwidth = 0.2, 
                 fill = "lightgrey", colour = "black", alpha = 0.7) +
  stat_function(fun = function(x) (B_n_samples1 * dnorm(x, mean(B_dist1), sd(B_dist1)) + B_n_samples2 * dnorm(x, mean(B_dist2), sd(B_dist2))) * 0.2,
                colour = "darkblue", linewidth = 1, n = 1000) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(0, NA))

# 4. Plot the two overlapping distributions separately to see how they combine
B_p2 <- ggplot(B_theo_data, aes(x = value, y = count, fill = source)) +
  geom_area(alpha = 0.5, colour = "black", position = "identity") +
  geom_vline(data = B_stats, aes(xintercept = mean, colour = source), show.legend = FALSE) +
  geom_vline(data = B_stats, aes(xintercept = mean + sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(data = B_stats, aes(xintercept = mean - sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(xintercept = B_intersection_points, colour = "black", alpha = 0.7, show.legend = FALSE) +
  annotate("text", x = B_intersection_points, y = -70, label = "Root", colour = "black", angle = 90, vjust = -0.5, size = 4.8) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  scale_colour_manual(values = c("skyblue", "salmon")) +
  labs(
    fill = "Distribution",
    colour = "Distribution"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(xlim = c(0, NA))

png("B_10.png", width = 1000, height = 800)
grid.arrange(B_p1, B_p2, ncol = 1, 
             top = grid::textGrob("Gene Expression from Two Overlapping Distributions (B - 10%)", gp = grid::gpar(fontsize = 20, fontface = "bold")),
             left = grid::textGrob("Frequency", rot = 90, gp = grid::gpar(fontsize = 20), just = "center", x = unit(0.5, "npc")),
             bottom = grid::textGrob("Expression", gp = grid::gpar(fontsize = 20)),
             padding = unit(2, "line"))
dev.off()


## -- B3 --
# repeat the B analysis with the second distribution having 3% of the mass of the first distribution,

# Where the mass of the second distribution is 10% of the first
B3_n_samples1 <- 1000
B3_n_samples2 <- 30

# Generate the data
B3_dist1 <- rnorm(B3_n_samples1, mean = 0, sd = 0.3)
B3_dist2 <- rnorm(B3_n_samples2, mean = 1.4, sd = 0.35)

# ComB3ine B3oth distriB3utions into a single long-format data frame
B3_B3imodal_data <- data.frame(
  value = c(B3_dist1, B3_dist2),
  source = c(rep("Not Expressed", B3_n_samples1), 
             rep("Expressed", B3_n_samples2))
)

# Calculate means and SDs for each distriB3ution
B3_stats <- data.frame(
  source = c("Not Expressed", "Expressed"),
  mean = c(mean(B3_dist1), mean(B3_dist2)),
  sd = c(sd(B3_dist1), sd(B3_dist2))
)

# Calculate the intersection points
B3_intersection_points <- find_intersection(mean(B3_dist1), sd(B3_dist1), mean(B3_dist2), sd(B3_dist2), n1 = B3_n_samples1, n2 = B3_n_samples2)

# Generate theoretical normal curves for B3
B3_x_range <- seq(min(B3_B3imodal_data$value), max(B3_B3imodal_data$value), length.out = 1000)
B3_theo_data <- data.frame(
  value = rep(B3_x_range, 2),
  count = c(dnorm(B3_x_range, mean(B3_dist1), sd(B3_dist1)) * B3_n_samples1,
            dnorm(B3_x_range, mean(B3_dist2), sd(B3_dist2)) * B3_n_samples2),
  source = rep(c("Not Expressed", "Expressed"), each = 1000)
)

# 3. Plot the overall B3imodal distriB3ution
B3_p1 <- ggplot(B3_B3imodal_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(count)), binwidth = 0.2, 
                 fill = "lightgrey", colour = "black", alpha = 0.7) +
  stat_function(fun = function(x) (B3_n_samples1 * dnorm(x, mean(B3_dist1), sd(B3_dist1)) + B3_n_samples2 * dnorm(x, mean(B3_dist2), sd(B3_dist2))) * 0.2,
                colour = "darkblue", linewidth = 1, n = 1000) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(0, NA))

# 4. Plot the two overlapping distriB3utions separately to see how they comB3ine
B3_p2 <- ggplot(B3_theo_data, aes(x = value, y = count, fill = source)) +
  geom_area(alpha = 0.5, colour = "black", position = "identity") +
  geom_vline(data = B3_stats, aes(xintercept = mean, colour = source), show.legend = FALSE) +
  geom_vline(data = B3_stats, aes(xintercept = mean + sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(data = B3_stats, aes(xintercept = mean - sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(xintercept = B3_intersection_points, colour = "black", alpha = 0.7, show.legend = FALSE) +
  annotate("text", x = B3_intersection_points, y = -70, label = "Root", colour = "black", angle = 90, vjust = -0.5, size = 4.8) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  scale_colour_manual(values = c("skyblue", "salmon")) +
  labs(
    fill = "Distribution",
    colour = "Distribution"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(xlim = c(0, NA))

png("B_3.png", width = 1000, height = 800)
grid.arrange(B3_p1, B3_p2, ncol = 1, 
             top = grid::textGrob("Gene Expression from Two Overlapping Distributions (B - 3%)", gp = grid::gpar(fontsize = 20, fontface = "bold")),
             left = grid::textGrob("Frequency", rot = 90, gp = grid::gpar(fontsize = 20), just = "center", x = unit(0.5, "npc")),
             bottom = grid::textGrob("Expression", gp = grid::gpar(fontsize = 20)),
             padding = unit(2, "line"))
dev.off()


## ----------- example C ---------------------
# Where the mass of the fisrt distribution is 10% of the first
C_n_samples1 <- 100
C_n_samples2 <- 1000

# Generate the data
C_dist1 <- rnorm(C_n_samples1, mean = 0, sd = 0.3)
C_dist2 <- rnorm(C_n_samples2, mean = 2, sd = 0.85)

# ComCine Coth distriCutions into a single long-format data frame
C_Cimodal_data <- data.frame(
  value = c(C_dist1, C_dist2),
  source = c(rep("Not Expressed", C_n_samples1), 
             rep("Expressed", C_n_samples2))
)

# Calculate means and SDs for each distriCution
C_stats <- data.frame(
  source = c("Not Expressed", "Expressed"),
  mean = c(mean(C_dist1), mean(C_dist2)),
  sd = c(sd(C_dist1), sd(C_dist2))
)

# Calculate the intersection points
C_intersection_points <- find_intersection(mean(C_dist1), sd(C_dist1), mean(C_dist2), sd(C_dist2), n1 = C_n_samples1, n2 = C_n_samples2)

# Generate theoretical normal curves for C
C_x_range <- seq(min(C_Cimodal_data$value), max(C_Cimodal_data$value), length.out = 1000)
C_theo_data <- data.frame(
  value = rep(C_x_range, 2),
  count = c(dnorm(C_x_range, mean(C_dist1), sd(C_dist1)) * C_n_samples1,
            dnorm(C_x_range, mean(C_dist2), sd(C_dist2)) * C_n_samples2),
  source = rep(c("Not Expressed", "Expressed"), each = 1000)
)

# 3. Plot the overall Cimodal distriCution
C_p1 <- ggplot(C_Cimodal_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(count)), binwidth = 0.2, 
                 fill = "lightgrey", colour = "black", alpha = 0.7) +
  stat_function(fun = function(x) (C_n_samples1 * dnorm(x, mean(C_dist1), sd(C_dist1)) + C_n_samples2 * dnorm(x, mean(C_dist2), sd(C_dist2))) * 0.2,
                colour = "darkblue", linewidth = 1, n = 1000) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(0, NA))

# 4. Plot the two overlapping distriCutions separately to see how they comCine
C_p2 <- ggplot(C_theo_data, aes(x = value, y = count, fill = source)) +
  geom_area(alpha = 0.5, colour = "black", position = "identity") +
  geom_vline(data = C_stats, aes(xintercept = mean, colour = source), show.legend = FALSE) +
  geom_vline(data = C_stats, aes(xintercept = mean + sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(data = C_stats, aes(xintercept = mean - sd, colour = source), linetype = "dotted", show.legend = FALSE) +
  geom_vline(xintercept = C_intersection_points, colour = "black", alpha = 0.7, show.legend = FALSE) +
  annotate("text", x = C_intersection_points, y = -6, label = "Root", colour = "black", angle = 90, vjust = -0.5, size = 4.8) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  scale_colour_manual(values = c("skyblue", "salmon")) +
  labs(
    fill = "Distribution",
    colour = "Distribution"
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(xlim = c(0, NA))

png("C.png", width = 1000, height = 800)
grid.arrange(C_p1, C_p2, ncol = 1, 
             top = grid::textGrob("Gene Expression from Two Overlapping Distributions (C)", gp = grid::gpar(fontsize = 20, fontface = "bold")),
             left = grid::textGrob("Frequency", rot = 90, gp = grid::gpar(fontsize = 20), just = "center", x = unit(0.5, "npc")),
             bottom = grid::textGrob("Expression", gp = grid::gpar(fontsize = 20)),
             padding = unit(2, "line"))
dev.off()

## ------------ Example D ---------------------
# where the SD of the two distributions overlap
D_n_samples1 <- 1000
D_n_samples2 <- 800

# Generate the data
D_dist1 <- rnorm(D_n_samples1, mean = 0, sd = 0.8)
D_dist2 <- rnorm(D_n_samples2, mean = 1.6, sd = 1.3)

# DomDine Doth distriDutions into a single long-format data frame
D_Dimodal_data <- data.frame(
  value = c(D_dist1, D_dist2),
  source = c(rep("Not Expressed", D_n_samples1), 
             rep("Expressed", D_n_samples2))
)

# Dalculate means and SDs for each distriDution
D_stats <- data.frame(
  source = c("Not Expressed", "Expressed"),
  mean = c(mean(D_dist1), mean(D_dist2)),
  sd = c(sd(D_dist1), sd(D_dist2))
)

# Dalculate the intersection points
D_intersection_points <- find_intersection(mean(D_dist1), sd(D_dist1), mean(D_dist2), sd(D_dist2), n1 = D_n_samples1, n2 = D_n_samples2)

# Generate theoretical normal curves for D
D_x_range <- seq(min(D_Dimodal_data$value), max(D_Dimodal_data$value), length.out = 1000)
D_theo_data <- data.frame(
  value = rep(D_x_range, 2),
  count = c(dnorm(D_x_range, mean(D_dist1), sd(D_dist1)) * D_n_samples1,
            dnorm(D_x_range, mean(D_dist2), sd(D_dist2)) * D_n_samples2),
  source = rep(c("Not Expressed", "Expressed"), each = 1000)
)

# 3. Plot the overall Dimodal distriDution
D_p1 <- ggplot(D_Dimodal_data, aes(x = value)) +
    geom_histogram(aes(y = after_stat(count)), binwidth = 0.2, 
                                 fill = "lightgrey", colour = "black", alpha = 0.7) +
    stat_function(fun = function(x) (D_n_samples1 * dnorm(x, mean(D_dist1), sd(D_dist1)) + D_n_samples2 * dnorm(x, mean(D_dist2), sd(D_dist2))) * 0.2,
                                colour = "darkblue", linewidth = 1, n = 1000) +
    theme_minimal(base_size = 16) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank()) +
    coord_cartesian(xlim = c(0, NA))

# 4. Plot the two overlapping distriDutions separately to see how they comDine
D_p2 <- ggplot(D_theo_data, aes(x = value, y = count, fill = source)) +
    geom_area(alpha = 0.5, colour = "black", position = "identity") +
    geom_vline(data = D_stats, aes(xintercept = mean, colour = source), show.legend = FALSE) +
    geom_vline(data = D_stats, aes(xintercept = mean + sd, colour = source), linetype = "dotted", show.legend = FALSE) +
    geom_vline(data = D_stats, aes(xintercept = mean - sd, colour = source), linetype = "dotted", show.legend = FALSE) +
    geom_vline(xintercept = D_intersection_points, colour = "black", alpha = 0.7, show.legend = FALSE) +
    annotate("text", x = D_intersection_points, y = -6, label = "Root", colour = "black", angle = 90, vjust = -0.5, size = 4.8) +
    theme_minimal(base_size = 16) +
    scale_fill_manual(values = c("skyblue", "salmon")) +
    scale_colour_manual(values = c("skyblue", "salmon")) +
    labs(
        fill = "Distribution",
        colour = "Distribution"
    ) +
    theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                legend.position = c(0.95, 0.95),
                legend.justification = c("right", "top"),
                legend.background = element_rect(fill = alpha("white", 0.5)),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                plot.margin = unit(c(0, 0, 0, 0), "lines")) +
    coord_cartesian(xlim = c(0, NA))

# png("D_Seperate.png", width = 800, height = 600)
# print(D_p2)
# dev.off()


png("D.png", width = 1000, height = 800)
grid.arrange(D_p1, D_p2, ncol = 1, 
                         top = grid::textGrob("Gene Expression from Two Overlapping Distributions (D)", gp = grid::gpar(fontsize = 20, fontface = "bold")),
                         left = grid::textGrob("Frequency", rot = 90, gp = grid::gpar(fontsize = 20), just = "center", x = unit(0.5, "npc")),
                         bottom = grid::textGrob("Expression", gp = grid::gpar(fontsize = 20)),
                         padding = unit(2, "line"))
dev.off()
