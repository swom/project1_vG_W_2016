library(ggplot2)
demographic = data.frame(cycle = c(), active = c(), spore = c(), sporu = c(), seed = c())

for (i in 1 : length(list.files(path = '.',pattern = "sim_demographic_s*")))
{
  replicate = read.csv(paste("sim_demographic_s",i,".csv", sep = ""))
  replicate$seed = i;
  demographic = rbind(replicate,demographic)
}

colnames(demographic)= c("cycle", "active", "spore", "sporu" , "seed")
demographic$seed = as.factor(demographic$seed) 
ggplot(demographic,aes(cycle,active)) + 
  geom_point()
