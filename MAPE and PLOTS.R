library(forecast)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
a <- c(1,2,16,17,18,21,55,56,57)
b <- a + 1
names <- colnames(train_set[,b])
f1 <- read.csv("SW_FORECAST.csv", header=F)
t1 <- read.csv("TESTSETNOSTATIONARY.csv", header=F)
f21 <- f1[,a]
t21 <- t1[,a]
f1 <- as.matrix(f21)
t1 <- as.matrix(t21)
diff1 <- t1-f1
mape1 <- abs(diff1) / abs(t1) * 100
mape1[is.infinite(mape1)] <- NA
n_horizons1 <- nrow(mape1)
n_variables1 <- ncol(mape1)
horizons1 <- 1:n_horizons1
TEST2 <- read.csv("TESTSETSTATIONARY.csv", header=F)
FHLR <- read.csv("FORECAST_FHLR.csv",header=F)
T2 <- TEST2[,a]
F2 <- FHLR[,a]
train2 <- train_set[258,b]
options(digits = 15)
colnames(F2) <- colnames(train2)
new2 <- rbind(train2,F2)
TS2 <- read.csv("TESTSETNOSTATIONARY.csv",header=F)
TS2 <- TS2[,a]
LAST2 <- as.data.frame(apply(new2, 2, cumsum))
LAST2 <- LAST2[-1,]
f2 <- as.matrix(LAST2)
t2 <- as.matrix(TS2)
diff2 <- t2-f2
mape2 <- abs(diff2) / abs(t2) * 100
mape2[is.infinite(mape2)] <- NA
n_horizons2 <- nrow(mape2)
n_variables2 <- ncol(mape2)
horizons2 <- 1:n_horizons2
TEST3 <- read.csv("TESTSETSTATIONARY.csv", header=F)
FHLZ<- read.csv("FORECAST_FHLZ.csv",header=F)
T3 <- TEST3[,a]
F3 <- FHLZ[,a]
train3 <- train_set[258,b]
options(digits = 15)
colnames(F3) <- colnames(train3)
new3 <- rbind(train3,F3)
TS3 <- read.csv("TESTNOSTATIONARY.csv",header=F)
TS3 <- TS3[,a]
LAST3 <- as.data.frame(apply(new3, 2, cumsum))
LAST3 <- LAST3[-1,]
f3 <- as.matrix(LAST3)
t3 <- as.matrix(TS3)
diff3 <- tsw3-f3
mape3 <- abs(diff3) / abs(tsw3) * 100
mape3[is.infinite(mape3)] <- NA
n_horizons3 <- nrow(mape3)
n_variables3 <- ncol(mape3)
horizons3 <- 1:n_horizons3
plot_list_enhanced <- list()
for (i in 1:5) {
  var_data <- data.frame(
    Horizon = rep(horizons1, 3),
    Method = rep(c("SW", "FHLR", "FHLZ"), each = n_horizons1),
    APE = c(mape1[, i], mape2[, i], mape3[, i])
  )
  p_var_enhanced <- ggplot(var_data, aes(x = Horizon, y = APE, color = Method)) +
    geom_line(size = 1.3, alpha = 0.8, na.rm = TRUE, 
              position = position_dodge(width = 0.05),
              aes(linetype = Method)) +
    geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE,
               
               position = position_dodge(width = 0.05),
               aes(shape = Method)) +
    scale_color_manual(values = c("SW" = "#1B9E77", "FHLR" = "#D95F02", "FHLZ" = "#7570B3")) +
    scale_linetype_manual(values = c("SW" = "solid", "FHLR" = "longdash", "FHLZ" = "dotdash")) +
    scale_shape_manual(values = c("SW" = 16, "FHLR" = 17, "FHLZ" = 15)) +
    scale_x_continuous(breaks = horizons1) +
    labs(title = paste("Absolute percentage error omparison for", var_names[i]),
         x = "Forecast Horizon (Quarters)", 
         y = "Forecast Error (%)",
         color = "Method", linetype = "Method", shape = "Method") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          legend.box = "horizontal") +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
  
  plot_list_enhanced[[i]] <- p_var_enhanced
}
all_data <- data.frame()
for (i in 6:9) {
  var_data <- data.frame(
    Variable = var_names[i],
    Horizon = rep(horizons1, 3),
    Method = rep(c("SW", "FHLR", "FHLZ"), each = n_horizons1),
    APE = c(mape1[, i], mape2[, i], mape3[, i]),
    Enhanced = i <= 7
  )
  all_data <- rbind(all_data, var_data)
}
all_data$Variable <- factor(all_data$Variable, levels = var_names)
plot_hybrid_multiples <- ggplot(all_data, aes(x = Horizon, y = APE, color = Method)) +
  geom_line(data = subset(all_data, Enhanced == TRUE), 
            size = 1.3, alpha = 0.8, na.rm = TRUE, 
            position = position_dodge(width = 0.05),
            aes(linetype = Method)) +
  geom_point(data = subset(all_data, Enhanced == TRUE), 
             size = 2.5, alpha = 0.9, na.rm = TRUE,
             position = position_dodge(width = 0.05),
             aes(shape = Method)) +
  geom_line(data = subset(all_data, Enhanced == FALSE), 
            size = 1.1, na.rm = TRUE) +
  geom_point(data = subset(all_data, Enhanced == FALSE), 
             size = 1.5, na.rm = TRUE) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("SW" = "#1B9E77", "FHLR" = "#D95F02", "FHLZ" = "#7570B3")) +
  scale_linetype_manual(values = c("SW" = "solid", "FHLR" = "longdash", "FHLZ" = "dotdash")) +
  scale_shape_manual(values = c("SW" = 16, "FHLR" = 17, "FHLZ" = 15)) +
  scale_x_continuous(breaks = horizons1) +
  labs(x = "Forecast Horizon (Quarters)", 
       y = "Forecast Error (%)",
       color = "Method", linetype = "Method", shape = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        legend.box = "horizontal") +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
print(plot_hybrid_multiples)