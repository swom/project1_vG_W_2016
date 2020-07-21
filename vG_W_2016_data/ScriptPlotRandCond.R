library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
rand_demographic = data.frame()
for (i in  list.files(path = '.',pattern = "random_cond_sim_demographic_s\\d+change_\\d+"))
{
  conditions = read.table(i, sep = ",")
  conditions$seed = sub( "^.*s(\\d+).*",'\\1', i);
  conditions$change = sub( "^.*_(\\d+).*",'\\1', i);
  colnames(conditions) = colnames(rand_demographic)
  rand_demographic = rbind(conditions,rand_demographic)
}

n_columns = ncol(rand_demographic)
rand_demographic = rand_demographic[,-c(5: (n_columns - 2) )]
colnames(rand_demographic)= c("condition",
                              "active",
                              "spore",
                              "sporu", 
                              "seed",
                              "change_freq")


rand_demographic <- pivot_longer(
  rand_demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)
rand_demographic$change_freq = as.factor(rand_demographic$change_freq)
rand_demographic$seed = as.factor(rand_demographic$seed) 
rand_demographic$condition = as.factor(rand_demographic$condition) 
rand_demographic$variable = as.factor(rand_demographic$variable)

  rand_demographic %>% 
  subset(variable == "spore") %>% 
  subset(change_freq == "1") %>% 
  ggplot(aes(condition, seed, fill = value)) + 
  geom_tile(color = "black", size = 0.5) 
