library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(dir,"/vG_W_2016_data",sep = ""))
rand_demographic = data.frame()

for (i in  list.files(path = '.',pattern = "random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+.*"))
{
  conditions = read.table(i, sep = ",")
  conditions$seed = sub( "^.*s(\\d+).*",'\\1', i)
  conditions$change = sub( "^.*change_(\\d+).*",'\\1', i)
  conditions$amplitude = sub( "^.*amplitude_(\\d+)",'\\1', i)
  conditions$amplitude = sub( "(\\d+).csv",'\\1', conditions$amplitude)
  rand_demographic = rbind(conditions,rand_demographic)
}

n_columns = ncol(rand_demographic)
rand_demographic = rand_demographic[,-c(5: (n_columns - 3) )]
colnames(rand_demographic)= c("condition",
                              "active",
                              "spore",
                              "sporu", 
                              "seed",
                              "change_freq",
                              "amplitude")


rand_demographic <- pivot_longer(
  rand_demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)

rand_demographic$amplitude = as.factor(rand_demographic$amplitude)
rand_demographic$change_freq = as.factor(rand_demographic$change_freq)
rand_demographic$seed = as.factor(rand_demographic$seed) 
rand_demographic$condition = as.factor(rand_demographic$condition) 
rand_demographic$variable = as.factor(rand_demographic$variable)

new_rand_demo = rand_demographic %>%
  group_by(condition,change_freq) %>%
  mutate(
    "ratio_value" = value / max(value)
  )

new_rand_demo %>% 
  subset(variable == "spore") %>% 
  subset(as.numeric(seed) < 51) %>% 
  subset(as.numeric(amplitude) < 2) %>% 
  ggplot(aes(condition, seed, fill = ratio_value)) + 
  geom_tile(color = "black", size = 0.5) +
  facet_grid(change_freq  ~ amplitude)

