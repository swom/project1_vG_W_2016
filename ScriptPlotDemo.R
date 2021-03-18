library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")
setwd(evo_dir)

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
colnames(demographic)= c("cycle",
                         "active",
                         "spore",
                         "sporu" ,
                         sprintf("env_p_%s",seq(1:(n_columns - 4 - 2))),
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

#Check that for change of freq 0 environment par are always the same
check = 
  demographic %>% 
  subset(change_freq == 0)

#should return 1 row dataframe
check_u = 
  check %>% 
  subset(cycle == max(cycle)) %>% 
  select(starts_with("env_p")) %>% 
  unique()

#plot
demographic %>% 
  ggplot(aes(x = cycle, y = value, color = variable)) + 
  geom_smooth()  +
  facet_grid(change_freq ~ seed)


ggsave("C:/Users/p288427/Desktop/research presentation/first_evo_demo.pdf",
       width = 100,
       height = 20, 
       units = "cm",
       limitsize = F)



