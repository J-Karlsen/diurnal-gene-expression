library(tidyverse)
library(reshape2)
library(ggplot2)

#load data 
cult1 <- read_csv("path_to>cultivation1_data.csv")
cult2 <- read_csv("path_to>cultivation2_data.csv")

# Select relevant measurements
cult1 <- cult1 %>% select("$Date", "$Time", FE1_BATCHTIME, FE1_OD, FE1_FEEDRATEA, FE1_DO, FE1_REDLIGHT,
                        FE1_BLUELIGHT, FE1_CO2FLOW) %>%
  rename("date" = "$Date",
         "clock" = "$Time",
         "time" = "FE1_BATCHTIME",
         "od" = "FE1_OD",
         "feed_rate" = "FE1_FEEDRATEA",
         "do" = "FE1_DO",
         "red_light" = "FE1_REDLIGHT",
         "blue_light" = "FE1_BLUELIGHT",
         "co2_flow" = "FE1_CO2FLOW") %>%
  mutate(clock = as.character(clock))

cult2 <- cult2 %>% select("$Date", "$Time", FE1_BATCHTIME, FE1_OD, FE1_FEEDRATEA, FE1_DO, FE1_REDLIGHT,
                          FE1_BLUELIGHT, FE1_CO2FLOW) %>%
  rename("date" = "$Date",
         "clock" = "$Time",
         "time" = "FE1_BATCHTIME",
         "od" = "FE1_OD",
         "feed_rate" = "FE1_FEEDRATEA",
         "do" = "FE1_DO",
         "red_light" = "FE1_REDLIGHT",
         "blue_light" = "FE1_BLUELIGHT",
         "co2_flow" = "FE1_CO2FLOW") %>%
  mutate(clock = as.character(clock))


# cult1: Select time span. 84 h before 1st sampling day (5 aug, 19:40) to 11 aug, 21:40
cult1 <- cult1 %>% filter(time >= 3098.194 & time <= 3244.193)

# cult2: Select time span. 84 h before 1st sampling day (15 dec, 21:00) to 21 dec, 21:40
cult2 <- cult2 %>% filter(time >= 83.36)


# Cult1: Set zero time point at start of first day of sampling (9 aug, 07:40)
t0 <- cult1 %>% filter(date == "08/09/17" & clock == "07:40:00") %>% select(time) %>% as.numeric 
cult1 <- cult1 %>% mutate(time = time - t0)

# Cult2: Set zero time point at start of first day of sampling (19 dec, 09:00)
t0 <- cult2 %>% filter(date == "12/19/18" & clock == "09:00:00") %>% select(time) %>% as.numeric 
cult2 <- cult2 %>% mutate(time = time - t0)


# Convert light intensities from lux to umol photons m-2 s-1.
# Replace neg. values in blue light -> 0 and values = 21.306 in red -> 0
cult1 <- cult1 %>% mutate(red_light = red_light*0.018+21.306,
                          blue_light =  blue_light*0.007-13.379) %>% 
  mutate(red_light = case_when(red_light == 21.306 ~ 0, TRUE ~ red_light),
         blue_light = case_when(blue_light < 0 ~ 0, TRUE ~ blue_light)
         )

cult2 <- cult2 %>% mutate(red_light = red_light*0.018+21.306,
                          blue_light =  blue_light*0.007-13.379) %>% 
  mutate(red_light = case_when(red_light == 21.306 ~ 0, TRUE ~ red_light),
         blue_light = case_when(blue_light < 0 ~ 0, TRUE ~ blue_light)
         )


# Convert OD from 880 nm -> 730nm
ods1 <- cult1 %>% filter(
  ((clock == "06:40:00" | clock == "08:40:00" | clock == "13:40:00") & date == "08/09/17") |
    ((clock == "13:40:00" | clock == "18:40:00" | clock == "20:40:00") & date == "08/10/17") |
    ((clock == "06:40:00" | clock == "08:40:00") & date == "08/11/17")
  ) %>%
  select(date, clock, od) %>% 
  mutate(od730 = c(0.77, 0.762, 0.859, 0.830, 0.809, 0.759, 0.772, 0.740))
conv_factor1 <- ods1 %>% mutate(conv = od730/od) %>% summarise(conv_mean = mean(conv)) %>%
  as.numeric
# conv_factor1
# [1] 0.1247875
cult1 <- cult1 %>% mutate(od = od*conv_factor1)

ods2 <- cult2 %>% filter(
  ((clock == "08:00:00" | clock == "10:00:00" | clock == "20:00:00" | clock == "22:00:00") &
     date == "12/19/18") |
    ((clock == "20:00:00" | clock == "22:00:00") & date == "12/20/18") |
    ((clock == "08:00:00" | clock == "10:00:00") & date == "12/21/18")
) %>%
  select(date, clock, od) %>% 
  mutate(od730 = c(0.621, 0.624, 0.718, 0.702, 0.704, 0.677, 0.643, 0.652))
conv_factor2 <- ods2 %>% mutate(conv = od730/od) %>% summarise(conv_mean = mean(conv)) %>%
  as.numeric
# conv_factor2
# [1] 0.08987442
cult2 <- cult2 %>% mutate(od = od*conv_factor2)


# Convert feed rate to dilution rate (h-1). Replace negative values with 0.
cult1 <- cult1 %>%
  mutate(dil_rate = (feed_rate*4.105-1.866)/1600, .after = od) %>%
  mutate(dil_rate = case_when(dil_rate < 0 ~ 0, TRUE ~ dil_rate))

cult2 <- cult2 %>%
  mutate(dil_rate = (feed_rate*4.105-1.866)/1600, .after = od) %>%
  mutate(dil_rate = case_when(dil_rate < 0 ~ 0, TRUE ~ dil_rate))


# Remove noise from dilution rate
for(cycles in 1:2) {
  cult1 <- cult1 %>% mutate(dil_rate = zoo::rollapply(cult1$dil_rate, 480,
                                                    mean, trim = 0.4, na.pad = T))
}

for(cycles in 1:2) {
  cult2 <- cult2 %>% mutate(dil_rate = zoo::rollapply(cult2$dil_rate, 480,
                                                      mean, trim = 0.4, na.pad = T))
}


# Remove noise from OD
cult1 <- cult1 %>% mutate(od = zoo::rollapply(cult1$od, 2*60*8,
                                            mean, trim = 0.4, na.pad = T))
cult1 <- cult1 %>% mutate(od = zoo::rollapply(cult1$od, 360,
                                            mean, trim = 0, na.pad = T))

cult2 <- cult2 %>% mutate(od = zoo::rollapply(cult2$od, 2*60*8,
                                              mean, trim = 0.4, na.pad = T))
cult2 <- cult2 %>% mutate(od = zoo::rollapply(cult2$od, 360,
                                              mean, trim = 0, na.pad = T))


# Calculate OD rate of change for µ calculation (OD/h)
cult1 <- cult1  %>% mutate(od_chg = zoo::rollapply(cult1$od, 121, na.pad = T,
                                                   function(x) {(x[121] - x[1]) / (120/120)}
))
cult2 <- cult2  %>% mutate(od_chg = zoo::rollapply(cult2$od, 61, na.pad = T,
                                                   function(x) {(x[61] - x[1]) / (60/60)}
))


# Calculate of growth rate = (od_chg + D*OD) / OD (h-1)
cult1 <- cult1  %>% mutate(mu = (od_chg + dil_rate*od) / od ) 

cult2 <- cult2  %>% mutate(mu = (od_chg + dil_rate*od) / od )


# Remove noise from growth rate
cult1 <- cult1 %>% mutate(mu_filt = zoo::rollapply(cult1$mu, 2*60*6,
                                               mean, trim = 0.4, na.pad = T)) %>%
  mutate(mu_filt = case_when(mu_filt < 0 ~ 0, TRUE ~ mu_filt))

cult2 <- cult2 %>% mutate(mu_filt = zoo::rollapply(cult2$mu, 2*60*6,
                                                   mean, trim = 0.4, na.pad = T)) %>%
  mutate(mu_filt = case_when(mu_filt < 0 ~ 0, TRUE ~ mu_filt))


# Calculate average growth rate
av_mu1 <- cult1 %>% filter(time > 0 & time < 48) %>% summarise(avg = mean(mu_filt)) %>% as.numeric
av_mu2 <- cult2 %>% filter(time > 0 & time < 24) %>% summarise(avg = mean(mu_filt)) %>% as.numeric


# Adjust DO values i each cultivation to make up for difference in probe calibration and
# change 100% saturation  value from 40/50 -> 100
cult1 <- cult1 %>% mutate(do = do + 40)
cult2 <- cult2 %>% mutate(do = do + 50)


# Reduce to two data points/h before plotting
cult1 <- cult1[seq(1,nrow(cult1),60),]
cult2 <- cult2[seq(1,nrow(cult2),30),]


# Convert data tibbles to long format and prepare for plotting
cult1_l <- cult1 %>% select(-co2_flow, -feed_rate, -mu) %>% pivot_longer(cols = od:mu_filt,
                                                                         names_to = "variable",
                                                                         values_to = "value") %>%
  filter(time >= -72 & time <= 50) %>%
  mutate(variable = factor(variable, levels = c("red_light", "blue_light", "od", "od_chg",
                                                "dil_rate", "mu_filt", "do")))

cult2_l <- cult2 %>% select(-co2_flow, -feed_rate, -mu) %>% pivot_longer(cols = od:mu_filt,
                                                                         names_to = "variable",
                                                                         values_to = "value") %>%
  filter(time >= -72 & time <= 50) %>%
  mutate(variable = factor(variable, levels = c("red_light", "blue_light", "od", "od_chg",
                                                "dil_rate", "mu_filt", "do")))

# Plots (Fig S2)
c1 <- ggplot(cult1_l, aes(x = time, y = value, color = variable)) +
  geom_point(size = 0.1) +
  facet_grid(variable~., scales = "free") + 
  scale_x_continuous(breaks = seq(-72, 50, 12), minor_breaks = F,
                     expand = expansion(mult = c(0.02, 0.02))) + 
  scale_y_continuous(minor_breaks = F,
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_color_manual(values = c("firebrick1", "royalblue2", "darkolivegreen",
                                "darkolivegreen4", "cyan4", "cyan3", "mediumpurple3")) +
  theme_bw() +
  theme(legend.position = "none",
        line = element_line(size = 0.4),
        panel.border = element_rect(size = 0.2, fill = NA),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())

# ggsave(file = "path",
#            plot = c1, height=9.5, width=15.8, units = "cm")

c2 <- ggplot(cult2_l, aes(x = time, y = value, color = variable)) +
  geom_point(size = 0.1) +
  facet_grid(variable~., scales = "free") + 
  scale_x_continuous(breaks = seq(-72, 50, 12), minor_breaks = F,
                     expand = expansion(mult = c(0.02, 0.02))) + 
  scale_y_continuous(minor_breaks = F,
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_color_manual(values = c("firebrick1", "royalblue2", "darkolivegreen",
                                "darkolivegreen4", "cyan4", "cyan3", "mediumpurple3")) +
  theme_bw() +
  theme(legend.position = "none",
        line = element_line(size = 0.4),
        panel.border = element_rect(size = 0.2, fill = NA),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank()) 

# ggsave(file ="path",
#        plot = c2, height=9.5, width=15.8, units = "cm")