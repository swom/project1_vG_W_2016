library(broom)
library(dplyr)
library(gplots)
library(ggplot2)
library(pals)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
first = "C:/Users/p288427/Desktop/hd_rand_evo/test_first"
last = "C:/Users/p288427/Desktop/hd_rand_evo/test_last"

setwd(last)
rand_demographic = data.frame()

if(getwd() == first)
{
  pattern = "first_gen_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
}
if(getwd() == last)
{
  pattern = "last_gen_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
  
}

n_env = 1
for (i in  list.files(path = '.', pattern = pattern)
     )
  {
    if(file.size(i) <= 0) next()
    replicate = read.csv(i)
    replicate$seed = sub( "^.*sim_demographic_s(\\d+).*",'\\1', i);
    replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
    replicate$n_env = sub( "^.*cond_per_seq(\\d+).*",'\\1', i)
    replicate$seq = sub( "^.*seq_(\\d+).*",'\\1', i, perl = T)
    colnames(replicate) = colnames(rand_demographic)
    rand_demographic = rbind(replicate,rand_demographic)
  }



n_columns = ncol(rand_demographic)
colnames(rand_demographic)= c("condition",
                         "active",
                         "spore",
                         "sporu" ,
                         "n_timesteps",
                         sprintf("env_p_%s",seq(1:(n_columns - 5 - 4))),
                         "seed",
                         "change_freq",
                         "seq")

rand_demographic$seed = as.factor(rand_demographic$seed)
rand_demographic$change_freq = as.factor(rand_demographic$change_freq)
rand_demographic$condition = as.factor(rand_demographic$condition)
rand_demographic$n_timesteps = as.numeric(rand_demographic$n_timesteps)

last_rand_demographic = rand_demographic
save(last_rand_demographic, file = "last_test_rand_evo_demo.R")
####Create color plaettes####

n_colors = 100
rbg<- colorRampPalette(c("red", "blue", "green"))(n_colors)
cub_hel <- rev(cubehelix(n_colors))

###Create breaks
color_scale = seq(min(first_rand_demographic$spore),max(last_rand_demographic$spore),
                  length.out = n_colors/10)
###
first_rand_demographic$gen = as.factor("first")
last_rand_demographic$gen = as.factor("last")
test_demog = rbind(first_rand_demographic, last_rand_demographic)

ggplot(data =  test_demog) +
  geom_tile(aes(x = condition, y = seed, fill = spore))+
  scale_fill_gradientn("spore", colors = rbg, breaks = color_scale)+
  facet_grid(.  ~ gen)


new_rand_demo = rand_demographic %>%
  group_by(condition,change_freq) %>%
  mutate(
    "ratio_value" = value / max(value)
  )

new_rand_demo %>% 
  subset(variable == "spore") %>% 
  subset(as.numeric(seed) < 51) %>% 
  #subset(as.numeric(amplitude) < 2) %>% 
  ggplot(aes(condition, seed, fill = ratio_value)) + 
  geom_tile(color = "black", size = 0.5) +
  facet_grid(change_freq  ~ amplitude)

