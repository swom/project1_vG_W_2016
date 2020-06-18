library(tidyr)
library(stringr)
library(ggplot2)
demographic = data.frame(cycle = c(), active = c(), spore = c(), sporu = c(), seed = c())

for (i in  list.files(path = '.',pattern = "sim_demographic_s*"))
{
  replicate = read.csv(i)
  replicate$seed = str_extract(i, pattern = "(\\d+)");
  colnames(replicate) = colnames(demographic)
  demographic = rbind(replicate,demographic)
}

colnames(demographic)= c("cycle", "active", "spore", "sporu" , "seed")
demographic$seed = as.factor(demographic$seed) 

demographic <- pivot_longer(
  demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)
demographic$variable = as.factor(demographic$variable)

ggplot(demographic,aes(cycle,value, color = variable)) + 
  geom_point(alpha = 0.1) 

