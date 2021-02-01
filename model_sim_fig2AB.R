library(deSolve)
library(ggplot2)
library(tidyverse)

# Protein synthesis, Fs_j (External variable / forcing function)
synth_fc_amp <- 3.05
amplitude <- 1 - 2/(synth_fc_amp+1)       # fraction of protein synthesis
frequency <- 1/24                         # oscillations/h
phase <- 0.5*pi
constant <- 1
Fs_j <- function(t, a = amplitude, f = frequency, p = phase, c = constant) {
  a * sin(2*pi*f * t + p) + c
}

# Protein oscillation model
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dFp_j <- Fs_j(time) * (mu + kd_mean) - Fp_j * (mu + kd_j)
    
    # return rate of change
    list(c(dFp_j))
  })
}

# Define state variable and initial value
start_state <- c(Fp_j = 1)

# Define parameters and their values
parameters <- list(c(mu = 0.018, kd_j = 0.01, kd_mean = 0.01),
                   c(mu = 0.018, kd_j = 0.05, kd_mean = 0.05),
                   c(mu = 0.018, kd_j = 0.10, kd_mean = 0.10))

# Set time span and step interval
times <- seq(0, 24*10, by = 0.02)

# Simulate protein oscillations 
p_frac <- list()
for(i in 1:length(parameters)) {
  out <- ode(y = start_state, times = times, func = model, parms = parameters[[i]])
  out <- out %>% as_tibble %>% mutate(time = as.numeric(time),
                                      frac = as.numeric(Fp_j),
                                      frac_of = "Abundance") %>% mutate %>%
    select(-Fp_j)
  
  p_frac <- append(p_frac, list(out))
}
rm(out, i)

# Generate synthesis data for plotting
s_frac <- tibble(time = times,
                 frac = Fs_j(times),
                 frac_of = "Synthesis")

# Combine abundace and synthesis data in each simulation into one tibble
fraction_l <- list()
for(i in 1:length(p_frac)) {
  f <- bind_rows(p_frac[[i]], s_frac)
  fraction_l <- append(fraction_l, list(f))
}
rm(p_frac, s_frac, f, i)

# Create variable representing parameter values in each simulation  
for(i in 1:length(fraction_l)) {
  fraction_l[[i]] <- fraction_l[[i]] %>%
    mutate(muX = "mu = ", mu = parameters[[i]]["mu"],
           kd_jX = "kd_j = ", kd_j = parameters[[i]]["kd_j"],
           kd_meanX = "kd_mean = ", kd_mean = parameters[[i]]["kd_mean"]) %>%
    unite("mu", muX:mu, sep = "", remove = T) %>%
    unite("kd_j", kd_jX:kd_j, sep = "", remove = T) %>%
    unite("kd_mean", kd_meanX:kd_mean, sep = "", remove = T)
  
}

# Combine data from all simulations into one tibble for plotting 
fraction <- bind_rows(fraction_l) %>%
  mutate(mu = as.factor(mu),
         kd_j = as.factor(kd_j),
         kd_mean = as.factor(kd_mean),
         frac_of = factor(frac_of, levels = c("Synthesis", "Abundance"), ordered = T)) 

# Make new time axis to display steady state oscillations only. Convert time unit to days.
fraction <- fraction %>% mutate(time_d = (time - 168)/24 )

# Calculate relative protein amplitude (max/min)
max <- fraction %>% filter(kd_j == "kd_j = 0.1" & time_d > 0 & frac_of == "Abundance") %>%
  summarise(m = max(frac)) %>% as.numeric
min <- fraction %>% filter(kd_j == "kd_j = 0.1" & time_d > 0 & frac_of == "Abundance") %>%
   summarise(m = min(frac)) %>% as.numeric
prot_fc_amp <- max/min
#summary(filter(fraction, time_d > 0 & frac_of == "Abundance" & kd_j == "kd_j = 0.05"))

# Plots (for publication)
p <- ggplot(filter(fraction, kd_j == "kd_j = 0.1"),
            aes(x = time_d, y = frac, group = frac_of, color = frac_of, linetype = frac_of)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1.6),
                     breaks = seq(0, 2, 0.2),
                     expand = expansion(mult = 0)) +
  scale_x_continuous(limits = c(0, 3.6),
                     breaks = seq(0, 3, 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = c("#216b21ff", "#AE2708")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme_bw() +
  theme(legend.position = "none",
        line = element_line(size = 0.4),
        panel.border = element_rect(size = 0.2, fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  annotate(
    "rect", fill="black", alpha=0.1,
    xmin=0.5, xmax=1,
    ymin=0, ymax=1.6) +
  annotate(
    "rect", fill="black", alpha=0.1,
    xmin=1.5, xmax=2,
    ymin=0, ymax=1.6) +
  annotate(
    "rect",fill="black", alpha=0.1,
    xmin=2.5, xmax=3,
    ymin=0, ymax=1.6) +
  annotate(
    "rect", fill="white",
    xmin=3, xmax=3.6,
    ymin=0, ymax=1.6)

# ggsave(file ="path",
#        plot = p, height=5.1, width=5.6, units = "cm")

# ggsave(file ="path",
#        plot = p, height=5.1, width=5.6, units = "cm")