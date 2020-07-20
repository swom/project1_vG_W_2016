library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
rand_demographic = data.frame()
for (i in  list.files(path = '.',pattern = "random_cond_sim_demographic_s*"))
{
  conditions = read.table(i, sep = ",")
  conditions$seed = str_extract(i, pattern = "(\\d+)");
  colnames(conditions) = colnames(rand_demographic)
  rand_demographic = rbind(conditions,rand_demographic)
}

n_columns = ncol(rand_demographic)
rand_demographic = rand_demographic[,-c(5: (n_columns - 1) )]
colnames(rand_demographic)= c("condition", "active", "spore", "sporu" , "seed")


rand_demographic <- pivot_longer(
  rand_demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)
rand_demographic$seed = as.factor(rand_demographic$seed) 
rand_demographic$condition = as.factor(rand_demographic$condition) 
rand_demographic$variable = as.factor(rand_demographic$variable)

spore_rand_cond = rand_demographic %>% subset(variable == "spore")
ggplot(spore_rand_cond, aes(condition, seed, fill = value)) + 
  geom_tile(color = "black", size = 0.5) 
