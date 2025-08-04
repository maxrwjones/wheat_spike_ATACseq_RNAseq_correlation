suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))

### Read in DF of 'real' putative silencers
real_silencers <- read.table("silencer_SS_profiles.tsv", header=TRUE)

dim(real_silencers)
head(real_silencers)

### Read the 100 DFs of simulated putative silencers into a list
sim_silencers <- grep("simulated_silencer", dir(), value = TRUE) %>% map(~read_delim(.x, delim='\t', show_col_types = FALSE))

length(sim_silencers)
dim(sim_silencers[[1]])
head(sim_silencers[[1]])

### Test that the shuffled datasets all have different first peaks
### (as this would be consistent with peaks having been assigned random parent genes)
first_peaks <- c()
for (i in 1:100){
   first_peaks[i]  <- sim_silencers[[i]][1,2]
}
head(first_peaks)
length(unique(first_peaks))



mean_hist <- ggplot(real_silencers, aes(x=mean)) + 
  geom_histogram(fill="#f8766d", color="#e9ecef") +
    ggtitle("Distribution of SS means") +
    labs(x = "Mean SS across developmental stages",
         y = "Number of genes")


ggsave("real_distribution_SS_means_silencers.svg", height=12, width=12, unit="cm")

median_hist <- ggplot(real_silencers, aes(x=median)) + 
  geom_histogram(fill="#f68060", color="#e9ecef") +
    ggtitle("Distribution of SS medians")

plot_grid(mean_hist, median_hist, labels = "AUTO")


mean_hist <- ggplot(sim_silencers[[1]], aes(x=mean)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef") +
    ggtitle("Distribution of SS means")

median_hist <- ggplot(sim_silencers[[1]], aes(x=median)) + 
  geom_histogram(fill="#f68060", color="#e9ecef") +
    ggtitle("Distribution of SS medians")

plot_grid(mean_hist, median_hist, labels = "AUTO")


mean_hist <- ggplot(sim_silencers[[2]], aes(x=mean)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef") +
    ggtitle("Distribution of SS means")

median_hist <- ggplot(sim_silencers[[2]], aes(x=median)) + 
  geom_histogram(fill="#f68060", color="#e9ecef") +
    ggtitle("Distribution of SS medians")

plot_grid(mean_hist, median_hist, labels = "AUTO")


cat("Cumulative distribution function for real vs simulated SS medians")
plot(ecdf(real_silencers$median))
for (i in 1:10){
    lines(ecdf(sim_silencers[[i]]$median), col="red")
}

cat("Cumulative distribution function for real vs simulated SS means")
plot(ecdf(real_silencers$mean))
for (i in 1:10){
    lines(ecdf(sim_silencers[[i]]$mean), col="red")
}


means_of_simulated_means <- sapply(sim_silencers, function(df) mean(df$mean))
means_of_simulated_medians <- sapply(sim_silencers, function(df) median(df$mean))

### Define some variables for lines I want to overlay
real_data_mean <- mean(real_silencers$mean)
mean_sim <- mean(means_of_simulated_means)
mean_sim_plus_sd <- mean(means_of_simulated_means) + sd(means_of_simulated_means)
mean_sim_minus_sd <- mean(means_of_simulated_means) - sd(means_of_simulated_means)
mean_plus_5sd <- mean(means_of_simulated_means) + 5*sd(means_of_simulated_means)
mean_minus_5sd <- mean(means_of_simulated_means) - 5*sd(means_of_simulated_means)

### Plot the data
ggplot() + aes(means_of_simulated_means) +
    geom_histogram(binwidth=0.001, colour="#000000", fill="#0099F8") +
    geom_segment(aes(x=real_data_mean, y=0, xend=real_data_mean, yend=50), color="red", linetype="solid", size=1) +
    geom_segment(aes(x=mean_sim, y=0, xend=mean_sim, yend=50), color="deeppink", linetype="solid", size=1) +
    geom_segment(aes(x=mean_sim_plus_sd, y=0, xend=mean_sim_plus_sd, yend=50), color="deeppink", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_sim_minus_sd, y=0, xend=mean_sim_minus_sd, yend=50), color="deeppink", linetype="dashed", size=1) +
     geom_segment(aes(x=mean_plus_5sd, y=0, xend=mean_plus_5sd, yend=50), color="purple", linetype="dashed", size=1) +
     geom_segment(aes(x=mean_minus_5sd, y=0, xend=mean_minus_5sd, yend=50), color="purple", linetype="dashed", size=1) +
    labs(title = "Mean simulated SS means from 100 datasets",
         x = "Mean SS mean",
         y = "Count") +
    theme_classic() +
    theme(axis.title = element_text(size = rel(2.5)),
          axis.text = element_text(size = rel(2)),
          axis.ticks = element_line(linewidth = rel(2)),
          plot.title = element_text(size=rel(3))
         ) +
    scale_y_continuous(expand = c(0,0))

ggsave("simulated_SS_means_silencers.svg", height=12, width=50, unit="cm", limitsize = FALSE)


### Calculate some summary statistics
mean_sim
sd(means_of_simulated_means)
mean(means_of_simulated_medians)