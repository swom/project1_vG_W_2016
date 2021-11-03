######Setting up############
library(tidyr)
library(tidyverse)
library(tidytree)
library(ggplot2)
library(ggtree)
library(ape)
library(gganimate)
library(gridExtra)
library(R.utils)
library(treeio)

dir = "X:/project1_vG_W_2016"
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")

#########Reading######
getwd()
setwd(evo_dir)
funders_success = data.frame()
pattern = "^before_last_pop_s\\d+change_0"
if(pattern == "^before_last_pop_s\\d+change_0") for (i in list.files(path = '.',pattern = pattern)) {
    
  last_pop = read.csv(i, header = FALSE, sep = ',')
  last_pop$seed = as.numeric(sub( "^.*_s(\\d+).*",'\\1', i));
  funders_success = rbind(last_pop, funders_success)
  
  } else for (i in list.files(path = '.',pattern = "funders_success_s42change_0.csv")) {
  
    nl <- countLines(i)
    sim_start = read.csv(i, header = FALSE, sep = ',', nrows = 100, skip = 1)
    sim_mid =  read.csv(i, header = FALSE, sep = ',', nrows = 100, skip = nl/2 + (nl/2) %% 100 )
    sim_end = read.csv(i, header = FALSE, sep = ',', nrows = 100, skip = nl - 100)
    
    sim_run = rbind(sim_start, sim_mid, sim_end)
    sim_run$seed = as.numeric(sub( "^.*_s(\\d+).*",'\\1', i));
    funders_success = rbind(sim_run,funders_success)
}

#put some useful names to column
names(funders_success)[1]<-"cycle" 
names(funders_success)[2]<-"ID" 
names(funders_success)[length(names(funders_success)) - 1]<-"success" 
funders_success$seed = as.numeric(funders_success$seed)

##take away symbols for saving
funders_success = funders_success %>% select(where(is.numeric),ID,cycle,success,seed)
funder = c(499,rep(0,ncol(funders_success) - 5),0,0,0,0)
#####Make the tree#######

### Add funder (a network with all weights == 0), that will be used as root
funders_success = rbind(funder, funders_success)
njs = list();
for(nth_cycle in unique(funders_success$cycle))
{
  f = list()
  tree = nj(funders_success %>%
                subset(seed < 10) %>% 
                remove_rownames() %>% 
                column_to_rownames("ID") %>% 
                select(-c(cycle,
                          seed,
                          success)) %>% 
                dist(method ="euclidean",
                     diag = T,
                     upper = T))
  f = as_tibble(tree) %>% left_join(funders_success %>% select(seed, ID, success, cycle) %>% rename(label = ID))
  tips_seed_2 = f %>% filter(seed == 0) %>% select(node)
  ff = as_tibble(root(as.treedata(f), as.vector(tips_seed_2$node))) %>% 
    left_join(f, by = "label") %>% 
    select(-c(node.y, parent.y, branch.length.y)) %>% 
    rename(parent = parent.x, node = node.x, branch.length = branch.length.x)
  njs[[length(njs) + 1]] = ff
}

#######Create the table of nodes that distinguish the different clades (seeds)#######
clades = array()
for (i in unique(ff$seed)){
  if(!is.na(i))
    {
      b = as.tibble(MRCA(as.treedata(ff), unlist(ff %>% filter(seed == i) %>% select(node)))) %>% 
        mutate(seed = i) %>% 
        rename(node = value)
      
     clades = rbind(clades, b)
  }
  else{
    clades = clades[-1,]
  }
}    

#####[need data from tot_dem in plot_rand_evo_demographic.R]#####
death_death = "C:/Users/p288427/Desktop/hd_rand_evo/death_death"
load( paste(death_death,"tot_dem.R",sep = "/"))
tot_dem$type = as.factor(tot_dem$type)
# extract avg_success for one evo run of a seed in one condition (called seq in this file)
tot_avgs = 
  tot_dem %>%
  rename(seq = condition) %>% 
  select(cycle, env_type, seed, seq, type, success) %>% 
  group_by(seq, seed, type, env_type) %>% 
  mutate(improvement = success[length(success)] - success[1]) %>% 
  mutate(start = success[1]) %>% 
  group_by(seq, seed, type) %>% 
  summarise(avg_start = mean(start), 
            avg_improvement = mean(improvement),
            avg_success = mean(success)) 
tot_avgs$type = tot_avgs$type %>% recode_factor(sel_spores = "sel")

tot_first_last_evo_suc = total %>% 
  filter( type != "eden") %>% 
  mutate(
    seq = as.factor(seq),
    seed = as.factor(seed),
    type = as.factor(type)) %>% 
  inner_join(tot_avgs) %>% 
  rename(avg_pre_adap = avg_start)

####continue plotting tree####

clades = clades %>% 
  mutate(seed = as.factor(seed))%>%
  left_join(
    tot_first_last_evo_suc %>% 
      filter(type == "seqex") %>% 
      group_by(seed) %>% 
      summarise(avg_suc = mean(avg_success))
    )

p1 = ggtree(as.treedata(ff),
            aes(color = as.factor(seed)),
            layout = "circular") +
  geom_treescale() +
  geom_rootpoint() +
  # geom_tippoint(aes(label = seed)) +
  geom_highlight(data=clades, aes(node=node, fill= avg_suc)) +
  scale_fill_gradientn("success", colors = rbg)





