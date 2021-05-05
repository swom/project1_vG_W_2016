library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(gganimate)

dir = "X:/project1_vG_W_2016"
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")

getwd()
setwd(evo_dir)
funders_success = data.frame()
for (i in list.files(path = '.',pattern = "funders_success_s42change_0.csv"))
{
  sim_run = read.csv(i, header = FALSE, sep = ',')
  funders_success = rbind(sim_run,funders_success)
}

#put some useful names to column
names(funders_success)[1]<-"cycle" 
names(funders_success)[2]<-"ID" 
names(funders_success)[length(names(funders_success))]<-"success" 

##take away symbols for saving
funders_success = funders_success %>% select(where(is.numeric),ID,cycle,success)

#calculate distance matrix for one cycle
subs = funders_success %>% subset(cycle > 0)

njs = list();
for(nth_cycle in unique(subs$cycle))
{
  f = list()
  f$tree = nj(subs %>%
    subset(cycle == nth_cycle) %>% 
    select(-c(cycle,ID,success))%>% 
    dist(method ="euclidean",
         diag = T,
         upper = T))
  f$cycle = nth_cycle
}


###Plot ###
  ggtree(nj(subs %>%
            select(-c(cycle,ID,success)) %>% 
            dist(method ="euclidean",
                 diag = T,
                 upper = T)), 
       layout = "circular") 

  
  
  trees <- lapply(c(10, 20, 40), rtree)
  class(trees) <- "multiPhylo"
  p  = ggtree( trees, aes(frame = .id))


