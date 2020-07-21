library(tidyr)
library(stringr)
library(ggplot2)
demographic = data.frame()

for (i in  list.files(path = '.',pattern = "sim_demographic_s\\d+change_\\d+"))
{
  replicate = read.csv(i)
  replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
  replicate$change = sub( "^.*_(\\d+).*",'\\1', i)
  colnames(replicate) = colnames(demographic)
  demographic = rbind(replicate,demographic)
}

n_columns = ncol(demographic)
demographic = demographic[,-c(5: (n_columns - 2) )]
colnames(demographic)= c("cycle",
                         "active",
                         "spore",
                         "sporu" ,
                         "seed",
                         "change_freq")

demographic <- pivot_longer(
  demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)
demographic$seed = as.factor(demographic$seed) 
demographic$change_freq = as.factor(demographic$change_freq) 
demographic$variable = as.factor(demographic$variable)

ggplot(demographic,aes(cycle,value, color = variable)) + 
  geom_point(alpha = 0.1) 

