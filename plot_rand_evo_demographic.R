library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pals)
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(dir,"/vG_W_2016_data/rand_evo",sep = ""))
demographic = data.frame()


for (i in  list.files(path = '.',pattern = "rand_evo_a3.000000cond_\\d+sim_demographic_s\\d+change_\\d+"))
{
  replicate = read.csv(i)
  replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
  replicate$change = sub( "^.*_(\\d+).*",'\\1', i)
  replicate$condition = sub( "^.*cond_(\\d+).*",'\\1', i, perl = T)
  colnames(replicate) = colnames(demographic)
  demographic = rbind(replicate,demographic)
}

n_columns = ncol(demographic)
demographic = demographic[,-c(6: (n_columns - 3) )]
colnames(demographic)= c("cycle",
                         "active",
                         "spore",
                         "sporu" ,
                         "n_timesteps",
                         "seed",
                         "change_freq",
                         "condition"
)

demographic <- pivot_longer(
  demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)
demographic$seed = as.factor(demographic$seed) 
demographic$change_freq = as.factor(demographic$change_freq) 
demographic$variable = as.factor(demographic$variable)
demographic$condition = as.factor(demographic$condition)
demographic$n_timesteps = as.factor(demographic$n_timesteps)

# create new columns for ratio of spore produced and starting production
demographic = demographic %>%
  subset(variable == "spore") %>% 
  group_by(condition) %>% 
  mutate(
    "ratio_value" = value / max(value)
  ) %>% 
  ungroup() %>% 
  group_by(condition, seed) %>% 
  mutate("ratio_start_production" = ratio_value[min(cycle)]) %>% 
  mutate("start_production" = value[min(cycle)])


ggplot(data = o) +
  geom_rect(data=o, aes(ymin=min(value), ymax= max(value + 1), xmin=min(as.numeric(cycle)),
                        xmax= max(as.numeric(cycle)), fill = ratio_start_production), alpha =0.5) + 
  geom_smooth(aes(cycle,value), method='lm', formula= y~x) +
  scale_fill_gradientn(colors = cubehelix(10)) +
  theme_classic()
facet_grid(seed ~ condition)

demographic %>% 
  subset(cycle == 1)  %>%
  ggplot(aes(condition, seed, fill = ratio_value)) + 
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradientn(colors = cubehelix(10))



