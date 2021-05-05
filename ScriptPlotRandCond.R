library(broom)
library(dplyr)
library(gplots)
library(ggplot2)
library("ggpubr")
library(pals)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
first_s = "C:/Users/p288427/Desktop/hd_rand_evo/test_first"
last_s = "C:/Users/p288427/Desktop/hd_rand_evo/test_last"
first_death = "C:/Users/p288427/Desktop/hd_rand_evo/test_first_death"
last_death = "C:/Users/p288427/Desktop/hd_rand_evo/test_last_death"
first_eden = "C:/Users/p288427/Desktop/hd_rand_evo/test_first_eden"
last_eden = "C:/Users/p288427/Desktop/hd_rand_evo/test_last_eden"
folders = c(first_s, last_s, first_death, last_death, first_eden, last_eden)
types = c("normal","eden", "death")

####reading####
for(j in folders){
  setwd(j)
  rand_demographic = data.frame()
  
  if(getwd() == first_s){
    
    pattern = "first_gen_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
  } else if(getwd() == last_s){
    
    pattern = "last_gen_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
  } else if(getwd() == first_death){
    
    pattern = "first_gen_death_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
  } else if(getwd() == last_death){
    
    pattern = "last_gen_death_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
  } else if(getwd() == first_eden){
    
    pattern = "first_gen_eden_death_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
  } else if(getwd() == last_eden){
    
    pattern = "last_gen_eden_death_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
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
                                "seq",
                                "cycle")
  
  rand_demographic$seed = as.factor(rand_demographic$seed)
  rand_demographic$change_freq = as.factor(rand_demographic$change_freq)
  rand_demographic$condition = as.factor(rand_demographic$condition)
  rand_demographic$n_timesteps = as.numeric(rand_demographic$n_timesteps)
  
  
  rand_demographic =  rand_demographic %>% 
    group_by(seed, condition) %>% 
    mutate("success" = spore * 10 + (sporu * 0 + active))
  
  if(getwd() == first_s){
    
    first_rand_demographic = rand_demographic
    save(first_rand_demographic, file = "first_test_rand_evo_demo.R")
    
  } else if (getwd() == last_s){
    
    last_rand_demographic = rand_demographic
    save(last_rand_demographic, file = "last_test_rand_evo_demo.R")
    
  }else if (getwd() == first_death){
    
    first_death_rand_demographic = rand_demographic
    save(first_death_rand_demographic, file = "first_death_test_rand_evo_demo.R")
    
  } else if (getwd() == last_death){
    
    last_death_rand_demographic = rand_demographic
    save(last_death_rand_demographic, file = "last_death_test_rand_evo_demo.R")
    
  }else if (getwd() == first_eden){
    
    first_eden_rand_demographic = rand_demographic
    save(first_eden_rand_demographic, file = "first_eden_test_rand_evo_demo.R")
    
  } else if (getwd() == last_eden){
    
    last_eden_rand_demographic = rand_demographic
    save(last_eden_rand_demographic, file = "last_eden_test_rand_evo_demo.R")
    
  }
  print(j)
  
}
####loading####

for(j in folders){
  setwd(j)
  if(getwd() == first_s){
    load("first_test_rand_evo_demo.R")
    
  } else if (getwd() == last_s){
    
    load("last_test_rand_evo_demo.R")
    
  }else if (getwd() == first_death){
    
    load("first_death_test_rand_evo_demo.R")
    
  } else if (getwd() == last_death){
    
    load("last_death_test_rand_evo_demo.R")
    
  }else if (getwd() == first_eden){
    
    load("first_eden_test_rand_evo_demo.R")
    
  } else if (getwd() == last_eden){
    
    load("last_eden_test_rand_evo_demo.R")
    
  }
}
####Create color plaettes####

n_colors = 100
rbg<- colorRampPalette(c("red", "blue", "green"))(n_colors)
cub_hel <- rev(cubehelix(n_colors))

######Plot first and last####
# types = c("death")

for( i in types)
{
  if(i == "normal")
  {
    first = first_rand_demographic
    last = last_rand_demographic
  } else if ( i == "eden")
  {
    first = first_eden_rand_demographic
    last = last_eden_rand_demographic
  } else if ( i == "death")
  {
    first = first_death_rand_demographic
    last = last_death_rand_demographic
  }
  
  ###Create breaks
  color_scale_spore = seq(min(first$spore),max(last$spore),
                          length.out = n_colors/10)
  
  color_scale_success = seq(min(first$success),max(last$success),
                            length.out = n_colors/10)
  ###uniting in a single data_frame
  first$gen = as.factor("first")
  last$gen = as.factor("last")
  test_demog = rbind(first, last)
  ###saving united dataframes
  if(i == "normal"){
    
    setwd(first_s)
    seqex_test_demog = test_demog
    save(seqex_test_demog, file = "seqex_test_demog.R")
    
  } else if (i == "death"){
    
    setwd(first_death)
    death_test_demog = test_demog
    save(death_test_demog, file = "death_test_demog.R") 
  } else if (i == "eden"){
    
    setwd(first_eden)
    eden_test_demog = test_demog
    save(eden_test_demog, file = "eden_test_demog.R") 
  }
  ####Check envs are the same####
  
  ggplot(data = test_demog) +
    geom_tile(aes(condition, seed, fill = as.factor(env_p_3) )) +
    facet_grid(gen  ~ .) +
    theme(legend.position = "none") 
  
  ####Plot dem####
  sp = ggplot(data =  test_demog) +
    geom_tile(aes(x = condition, y = seed, fill = spore))+
    scale_fill_gradientn("spore", colors = rbg, breaks = color_scale_spore)+
    ggtitle(paste(i,"_spore",sep = ""))+
    facet_grid(.  ~ gen)
  
  print(sp)
  
  su = ggplot(data =  test_demog) +
    geom_tile(aes(x = condition, y = seed, fill = success))+
    scale_fill_gradientn("success", colors = rbg, breaks = color_scale_success)+
    ggtitle(paste(i,"_success",sep = ""))+
    facet_grid(.  ~ gen)
  
  print(su)
  
  scaled = ggplot(data =  test_demog %>%  
                    pivot_longer(c(spore, success)) %>%
                    subset(seed != "9") %>% 
                    group_by(name) %>% 
                    mutate(value = scale(value))) +
    geom_tile(aes(x = condition, y = seed, fill = value))+
    scale_fill_gradientn("value", colors = rbg)+
    ggtitle(paste("scaled_all_",i, sep = ""))+
    facet_grid(name  ~ gen)
  
  print(scaled)
}

####Old plotting####
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



####Correlation tests####
####1-Loading####
for( i in types)
{
  if(i == "normal")
  {
    
    setwd(first_s)
    load("seqex_test_demog.R")
    test_demog = seqex_test_demog
    
  } else if ( i == "eden")
  {
    
    setwd(first_eden)
    load("eden_test_demog.R")
    test_demog = eden_test_demog
    
  } else if ( i == "death")
  {
  
    setwd(first_death)
    load("death_test_demog.R")
    test_demog = death_test_demog
  
  }
 
  cor_df = test_demog %>% 
    select(gen,seed,success, condition) %>% 
    group_by(seed, gen) %>% 
    summarise(avg_suc = mean(success)) %>% 
    pivot_wider(names_from = gen, values_from = avg_suc)
  
  print(paste("#############",i,"########################################"))

  #scatter_plot
  sp = ggscatter(cor_df, x = "first", y = "last", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            main = i)
  print(sp)
  # Shapiro-Wilk normality test for mpg
  sf = shapiro.test(cor_df$first) 
  print(sf)
  # Shapiro-Wilk normality test for wt
  sl = shapiro.test(cor_df$last) 
  print(sl)
  # first
  pf = ggqqplot(cor_df$first, ylab = "first", main = i)
  print(pf)
  # wt
  pl = ggqqplot(cor_df$last, ylab = "last", main = i)
  print(pl)
  #Pearson
  p = cor.test(cor_df$first, cor_df$last, 
           method = "pearson")
  print(p)
  #Spearman rho (rank test)
  s = cor.test(cor_df$first, cor_df$last, 
           method = "spearman")
  print(s)
  print(paste("#####################################################"))
  print(paste("#####################################################"))
  print("")
}
####Free style####
death_test_demog %>%
  select(seed, condition, env_p_3) %>% 
  ggplot()+
  geom_tile(aes(x = condition, y = seed, fill = as.factor(env_p_3)))+
  theme(legend.position = "none") 

save_local = local_run_first_gen
local_run_first_gen = s9
local_run_first_gen$seed = 123456789;
local_run_first_gen$change = 0
local_run_first_gen$n_env = 50
local_run_first_gen$seq = 1
local_run_first_gen$gen = "first"

n_columns = ncol(local_run_first_gen)
colnames(local_run_first_gen)= c("condition",
                                 "active",
                                 "spore",
                                 "sporu" ,
                                 "n_timesteps",
                                 sprintf("env_p_%s",seq(1:(n_columns - 5 - 4))),
                                 "seed",
                                 "change_freq",
                                 "seq",
                                 "cycle")

local_run_first_gen =  local_run_first_gen %>% 
  group_by(seed, condition) %>% 
  mutate("success" = spore * 10 + (sporu * 0 + active)) 

s = rbind(death_test_demog,
 local_run_first_gen %>% 
  subset(condition > 1) %>% 
  mutate(condition = as.factor(condition)) %>% 
   mutate(seed = as.factor(seed)) %>% 
   mutate(change_freq =as.factor(change_freq)) %>% 
   mutate(seq = as.character(seq)) %>% 
   mutate(cycle = as.character(cycle))
)

ggplot(s %>% subset(gen = "first") %>%  subset(seed != "123456789")) +
  geom_tile(aes(x = condition, y = seed, fill = success)) +
  scale_fill_gradientn("success", colors = rbg)

pattern = "first_gen_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s9_change_\\d+_amplitude_\\d+"
getwd()
list.files()
s9 = read.csv("first_gen_death_rand_evo_extreme_a3.000000seq_1cond_per_seq50random_cond_sim_demographic_s9_change_0_amplitude_3.000000.csv")
