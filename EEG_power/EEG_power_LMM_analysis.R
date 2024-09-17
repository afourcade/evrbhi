# EEG Power LMM analysis
# this script does an LMM analysis of the effects of arousal and movement on the
# ROI EEG spectral power in the classical frequency bands (delta, theta, alpha, beta, gamma)
# 2024 by Antonin Fourcade


# load packages -----------------------------------------------------------

library(ggdist)
library(gghalves)
library(raincloudplots)
library(lme4)
library(lmerTest)
library(hrbrthemes)
library(viridis)
library(car)
library(readxl)
library(performance)
library(rstatix)
library(tidyverse)
library(report)

require(pacman)
pacman::p_load("tidyr", "lattice", "ggplot2", "dplyr", "patchwork", "afex", "emmeans", "ggpol", "psychReport", "knitr",
               "cowplot","readr", "rmarkdown","Rmisc","devtools","gghalves","plotrix", "here","car","lme4","svglite", "stargazer","influence.ME", "MuMIn","gridExtra")

# set paths ---------------------------------------------------------------

basePath <- file.path("E:", "NeVRo", "Data", "Frequency_Bands_Power")
data_dir <- "ROI_EEG_pow"
fig_dir <- "figures"
res_dir <- "results"
dataname <- "nvr_roi_eeg_features"

# plotting parameters ----------------------------------------------------

colors <- c('dodgerblue3', 'firebrick2')
fills <- c('dodgerblue3', 'firebrick2')
line_color = 'gray'
line_alpha = .4
size = 2
alpha = .6
font <- "Avenir"

# get data and format -----------------------------------------------------

data <- read.csv(file.path(basePath, paste0(dataname, ".csv")), sep = ",")

data <- data %>% 
  mutate(arousal = factor(arousal, levels = c("LA", "HA")),
         mov_cond = factor(mov_cond, levels = c("nomov", "mov")),
         SJ = factor(SJ, levels = unique(data$SJ))) %>% 
  pivot_longer(!c(SJ, mov_cond, arousal))

# loop LMM across features -------------------------------------------------

features <- unique(data$name)
features_to_keep <- c("delta", "theta", "alpha", "beta", "gamma") # features to analyze
lmmResults <- vector(mode = "list", length = length(features))

for (i in 1:length(features)){
  # create/set current path
  curDir <- file.path(basePath)
  if (!file.exists(file.path(curDir, fig_dir))) {dir.create(file.path(curDir, fig_dir))}
  if (!file.exists(file.path(curDir, res_dir))) {dir.create(file.path(curDir, res_dir))}
  
  # get data of current features and log-transform
  if(features[i] %in% features_to_keep){
    curDat <- data %>% 
      filter(name == features[i]) %>% 
      mutate(value = log10(value))
    
    # histogram
    p <- ggplot(curDat, aes(x = value)) +
      geom_histogram(color = "black", alpha = .4, bins = 30) +
      xlab(features[i]) +
      theme_bw(base_size = 14)
    ggsave(file.path(curDir, paste0(fig_dir, "/", features[i], "_hist.png")), p, height = 4, width = 5)
    
    # create lmm
    curLMM <- lmer(value ~ arousal * mov_cond + (arousal * mov_cond | SJ), data = curDat)
    summary(curLMM)
    
    # pairwise differences of emmeans for all factors in a linear mixed model
    # In an analysis of covariance model, emmeans (or estimated marginal means) are the group means after having controlled for a covariate
    #emmeans(curLMM, revpairwise ~ arousal, lmerTest.limit = nrow(curDat), pbkrtest.limit = nrow(curDat))
    anova(curLMM,  type = 3)
    emmeans(curLMM, revpairwise ~ arousal, lmerTest.limit = nrow(curDat))
    emmeans(curLMM, revpairwise ~ arousal | mov_cond, lmerTest.limit = nrow(curDat), adjust = "tukey")
    emmeans(curLMM, revpairwise ~ mov_cond, lmerTest.limit = nrow(curDat))
    
    # save results and model checks
    sink(file.path(curDir, paste0(res_dir, "/", features[i], "_lmer.txt")))
    cat(c("Results" , features[i], "\n"))
    cat("\n")
    print(summary(curLMM))
    cat("\n")
    print(anova(curLMM,  type = 3))
    cat("\n")
    print(emmeans(curLMM, revpairwise ~ arousal, lmerTest.limit = nrow(curDat)))
    cat("\n")
    print(emmeans(curLMM, revpairwise ~ arousal | mov_cond, lmerTest.limit = nrow(curDat), adjust = "tukey"))
    cat("\n")
    print(emmeans(curLMM, revpairwise ~ mov_cond, lmerTest.limit = nrow(curDat)))
    cat("\n")
    print(model_performance(curLMM))
    cat("\n")
    cat("Singularity:\n")
    print(check_singularity(curLMM))
    cat("\n")
    cat("Heteroscedasticity:\n")
    print(check_heteroscedasticity(curLMM))
    cat("\n")
    cat("Autocorrelation:\n")
    print(check_autocorrelation(curLMM))
    cat("\n")
    cat("Normality of residuals:\n")
    print(check_normality(curLMM))
    cat("\n")
    cat("Outliers:\n")
    print(check_outliers(curLMM))
    cat("\n")
    cat("Performance accuracy:\n")
    print(performance_accuracy(curLMM))
    sink()
    
    #save results table
    write.csv(report_table(curLMM), file.path(curDir, paste0(res_dir, "/", features[i], "_lmer.csv")), row.names = FALSE)
    
    # save model check figure
    png(file.path(curDir, paste0(fig_dir, "/", features[i], "_checkModel.png")), width = 8, height = 10, unit = "in", res = 300)
    print(check_model(curLMM))
    dev.off()
    
    # raincloud plot
    # create summary dataframe
    dataSummary <- curDat %>% 
      group_by(SJ, arousal, mov_cond) %>% 
      dplyr::summarise(n = n(),
                mean = mean(value, na.rm = T),
                sd = sd(value, na.rm = T),
                se = sd/sqrt(n))
    
    # create 1x1 dataframe for plotting
    df_1x1 <- data_1x1(
      array_1 = as.matrix((dataSummary %>% filter(arousal == "LA"))[,"mean"]),
      array_2 = as.matrix((dataSummary %>% filter(arousal == "HA"))[,"mean"]),
      jit_distance = .09,
      jit_seed = 321)
    df_1x1$mov_cond <- c(as.matrix((dataSummary %>% filter(arousal == "LA"))[,"mov_cond"]), as.matrix((dataSummary %>% filter(arousal == "HA"))[,"mov_cond"]))
    df_1x1 <- as.data.frame(df_1x1)
    
    # plotting
    raincloud <- ggplot(data = df_1x1) + 
      geom_point(data = df_1x1 %>% 
                   dplyr::filter(x_axis == "1"), aes(x = jit, y = y_axis), 
                 color = colors[1], fill = fills[1], size = size, 
                 alpha = alpha, show.legend = FALSE) + 
      geom_point(data = df_1x1 %>% 
                   dplyr::filter(x_axis == "2"), aes(x = jit, y = y_axis), 
                 color = colors[2], fill = fills[2], size = size, 
                 alpha = alpha, show.legend = FALSE) + 
      geom_line(aes(x = jit, y = y_axis, group = id), color = line_color, alpha = line_alpha, 
                show.legend = FALSE) + 
      geom_half_boxplot(data = df_1x1 %>% 
                          dplyr::filter(x_axis == "1"), aes(x = x_axis, y = y_axis), 
                        color = colors[1], fill = fills[1], position = position_nudge(x = 1.18), 
                        side = "r", outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, 
                        width = 0.2, alpha = alpha, show.legend = FALSE) + 
      geom_half_boxplot(data = df_1x1 %>% dplyr::filter(x_axis == "2"), aes(x = x_axis, y = y_axis), color = colors[2], 
                        fill = fills[2], position = position_nudge(x = .32), 
                        side = "r", outlier.shape = NA, center = TRUE, 
                        errorbar.draw = FALSE, width = 0.2, alpha = alpha, 
                        show.legend = FALSE) + 
      geom_half_violin(data = df_1x1 %>% 
                         dplyr::filter(x_axis == "1"), aes(x = x_axis, y = y_axis), 
                       color = colors[1], fill = fills[1], position = position_nudge(x = 1.47), 
                       side = "r", alpha = alpha, show.legend = FALSE) + 
      geom_half_violin(data = df_1x1 %>% dplyr::filter(x_axis == "2"), aes(x = x_axis, y = y_axis), color = colors[2], 
                       fill = fills[2], position = position_nudge(x = 0.47), 
                       side = "r", alpha = alpha, show.legend = FALSE) +
      scale_x_continuous(breaks = c(1, 2), 
                         labels = c("Low", "High"), 
                         limits = c(.75, 3)) +
      xlab("Arousal") + 
      ylab(paste0("Mean ", features[i], " power (log-transform)")) +
      theme_classic(base_size = 16) +
      theme(text = element_text(family = font),
            axis.text.x = element_text(face = "bold", color = colors),
            strip.text.x = element_text(face = "bold")) +
      facet_grid(.~factor(mov_cond, levels = c("nomov", "mov"), labels = c("No movement", "Movement")))
    
    # show plot
    raincloud
    
    # Save plot in .png and .svg
    ggsave(file.path(curDir, paste0(fig_dir, "/", features[i], "_raincloud.png")), raincloud, height = 4, width = 5)
    ggsave(file.path(curDir, paste0(fig_dir, "/", features[i], "_raincloud.svg")), raincloud, height = 4, width = 5)
  }  
}
