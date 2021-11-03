####Setting up#####
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
library(tidyselect)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")
hd_rand_evo = "C:/Users/p288427/Desktop/hd_rand_evo"
hd_rand_evo_no_upt = "C:/Users/p288427/Desktop/hd_rand_evo/nouptake"
sequence = "C:/Users/p288427/Desktop/hd_rand_evo/sequence"
sequence_extr = "C:/Users/p288427/Desktop/hd_rand_evo/sequence_extr"
death_eden = "C:/Users/p288427/Desktop/hd_rand_evo/death_eden"
death_death = "C:/Users/p288427/Desktop/hd_rand_evo/death_death"
sel_spores = "C:/Users/p288427/Desktop/hd_rand_evo/sel_spores"

#####read data####
for( j in c(sel_spores))
{
  print(j)
  setwd(j)
  demographic = data.frame()
  
  n_env = 1
  if(getwd() == sel_spores){
    pattern = "sel_sp_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+sim_demographic_s\\d+change_0"
    for (i in  list.files(path = '.',
                          pattern = pattern ))
    {
      if(file.size(i) <= 0) next()
      replicate = read.csv(i)
      replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
      replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
      replicate$n_env = sub( "^.*cond_per_seq(\\d+).*",'\\1', i)
      replicate$condition = sub( "^.*seq_(\\d+).*",'\\1', i, perl = T)
      colnames(replicate) = colnames(demographic)
      demographic = rbind(replicate,demographic)
    }
  } else if(getwd() == death_eden){
    pattern = "eden_death_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+sim_demographic_s\\d+change_0"
    for (i in  list.files(path = '.',
                          pattern = pattern ))
    {
      if(file.size(i) <= 0) next()
      replicate = read.csv(i)
      replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
      replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
      replicate$n_env = sub( "^.*cond_per_seq(\\d+).*",'\\1', i)
      replicate$condition = sub( "^.*seq_(\\d+).*",'\\1', i, perl = T)
      colnames(replicate) = colnames(demographic)
      demographic = rbind(replicate,demographic)
    }
  }else if(getwd() == death_death){
    pattern = "death_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+sim_demographic_s\\d+change_0"
    for (i in  list.files(path = '.',
                          pattern = pattern ))
    {
      if(file.size(i) <= 0) next()
      replicate = read.csv(i)
      replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
      replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
      replicate$n_env = sub( "^.*cond_per_seq(\\d+).*",'\\1', i)
      replicate$condition = sub( "^.*seq_(\\d+).*",'\\1', i, perl = T)
      colnames(replicate) = colnames(demographic)
      demographic = rbind(replicate,demographic)
    }
  }else if(getwd() == sequence){
    for (i in  list.files(path = '.',
                          pattern = "rand_evo_a3.000000seq_\\d+cond_per_seq\\d+sim_demographic_s\\d+change_\\d+"))
    {
      if(file.size(i) <= 0) next()
      replicate = read.csv(i)
      replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
      replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
      replicate$n_env = sub( "^.*cond_per_seq(\\d+).*",'\\1', i)
      replicate$condition = sub( "^.*seq_(\\d+).*",'\\1', i, perl = T)
      colnames(replicate) = colnames(demographic)
      demographic = rbind(replicate,demographic)
    }
  }else if(getwd() == sequence_extr){
    for (i in  list.files(path = '.',
                          pattern = "rand_evo_extreme_a3.000000seq_\\d+cond_per_seq\\d+sim_demographic_s\\d+change_\\d+"))
    {
      if(file.size(i) <= 0) next()
      replicate = read.csv(i)
      replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
      replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
      replicate$n_env = sub( "^.*cond_per_seq(\\d+).*",'\\1', i)
      replicate$condition = sub( "^.*seq_(\\d+).*",'\\1', i, perl = T)
      colnames(replicate) = colnames(demographic)
      demographic = rbind(replicate,demographic)
    }
  }else if(getwd() == hd_rand_evo || getwd() == hd_rand_evo_no_upt
           || getwd() == rand_evo_dir || getwd() == evo_dir){
    
    for (i in  list.files(path = '.',
                          pattern = "rand_evo_a3.000000cond_\\d+sim_demographic_s\\d+change_\\d+"))
    {
      if(file.size(i) <= 0) next()
      replicate = read.csv(i)
      replicate$seed = sub( "^.*s(\\d+).*",'\\1', i);
      replicate$change = sub( "^.*_(\\d+).*",'\\1', i)
      replicate$condition = sub( "^.*cond_(\\d+).*",'\\1', i, perl = T)
      colnames(replicate) = colnames(demographic)
      demographic = rbind(replicate,demographic)
    }
  }
  
  
  n_columns = ncol(demographic)
  colnames(demographic)= c("cycle",
                           "active",
                           "spore",
                           "sporu" ,
                           "n_timesteps",
                           sprintf("env_p_%s",seq(1:(n_columns - 5 - 4))),
                           "seed",
                           "change_freq",
                           "n_env",
                           "condition")
  
  
  demographic$seed = as.factor(demographic$seed)
  demographic$change_freq = as.factor(demographic$change_freq)
  demographic$condition = as.factor(demographic$condition)
  demographic$n_timesteps = as.numeric(demographic$n_timesteps)
  
  # create new columns for ratio of spore produced and starting production
  demographic = demographic %>%
    group_by(seed,condition,cycle, change_freq) %>%
    mutate(total_n = sum(c_across(c(spore, sporu, active)))) %>% 
    ungroup() %>%
    group_by(condition,cycle) %>%
    mutate("ratio_value" = spore / max(spore)) %>%
    ungroup() %>%
    group_by(condition, seed) %>%
    mutate("change_value" = spore / spore[min(cycle)] - 1) %>%
    mutate("ratio_start_production" = ratio_value[min(cycle)]) %>%
    mutate("ratio_end_production" = ratio_value[max(cycle)]) %>%
    mutate("delta_rv_start_end" = ratio_end_production - ratio_start_production) %>%
    mutate("start_production" = spore[min(cycle)])  %>%
    mutate(success = spore * 10 + (sporu * 0 + active)) %>% 
    ungroup() %>%
    group_by(condition) %>%
    mutate("standardized_delta_rv_start_end" = delta_rv_start_end / max(delta_rv_start_end)) %>%
    mutate("overall_r_value" = spore / max(spore)) %>% 
    mutate(env_type = as.factor(as.numeric(as.factor(env_p_3)))) %>%  
    ungroup()
  
  coeffs = demographic %>% 
    group_by(condition, seed, env_type) %>%
    do(coeffs = coefficients(lm(spore ~ cycle, data = .))) %>% 
    mutate(slope = as.numeric(coeffs[2][[1]])) %>% 
    mutate(intercept = as.numeric(coeffs[1][[1]]))%>%
    ungroup()
  
  coeffs_success = demographic %>% 
    group_by(condition, seed, env_type) %>%
    do(coeffs_success = coefficients(lm(success ~ cycle, data = .))) %>% 
    mutate("slope_success" = as.numeric(coeffs_success[2][[1]])) %>% 
    mutate("intercept_success" = as.numeric(coeffs_success[1][[1]]))
  
  demographic = demographic %>% left_join(coeffs)
  demographic = demographic %>% left_join(coeffs_success)
  
  if(getwd() == death_eden) {
    
    death_eden_demographic  = demographic
    save(death_eden_demographic, file = "death_eden_rand_evo_demo.R")
    
  } else if(getwd() == sel_spores) {
    
    sel_spore_demographic  = demographic
    save(sel_spore_demographic, file = "sel_spores_rand_evo_demo.R")
    
  }  else if(getwd() == death_death) {
    
    death_death_demographic  = demographic
    save(death_death_demographic, file = "death_death_rand_evo_demo.R")
    
  } else if(getwd()  == sequence) {
    
    seq_demographic  = demographic
    save(seq_demographic, file = "seq_rand_evo_demo.R")
    
  } else if(getwd()  == sequence_extr) {
    
    seqex_demographic  = demographic
    save(seqex_demographic, file = "seqex_rand_evo_demo.R")
    
  }
}
####load hd_rand ####
setwd(hd_rand_evo)
load(file = "hd_rand_evo_demo.R")

setwd(hd_rand_evo_no_upt)
load(file = "hd_rand_evo_demo_no_upt_r.R")

setwd(sequence)
load("seq_rand_evo_demo.R")

setwd(sequence_extr)
load("seqex_rand_evo_demo.R")

setwd(death_eden)
load("death_eden_rand_evo_demo.R")

setwd(death_death)
load("death_death_rand_evo_demo.R")

setwd(sel_spores)
load("sel_spores_rand_evo_demo.R")

dem = sel_spore_demographic
####saving all dfs binded toghther####
tot_dem = rbind(seqex_demographic %>% mutate(type = "seqex"),
                death_eden_demographic %>% mutate(type = "eden"),
                sel_spore_demographic %>% mutate(type = "sel_spores"),
                death_death_demographic %>% mutate(type = "death"))
setwd(death_death)
save(tot_dem, file ="tot_dem.R")
####loading totatl dataframe####
setwd(death_death)
load("tot_dem.R")
tot_dem$type = as.factor(tot_dem$type)
dem = tot_dem %>% subset(type == "death")
####Check quantiles of change_value####
#i.e. how much the population improved its score production from the beginning

#looking at last generation

#overall quantiles
dem %>%
  group_by(seed, condition) %>%
  subset(cycle == max(cycle)) %>% 
  ungroup() %>% 
  summarise(qntl = quantile(change_value))

#quantiles per seed
qs = dem %>%
  group_by(seed, condition) %>%
  subset(cycle == max(cycle)) %>% 
  ungroup() %>% 
  group_by(seed) %>% 
  summarise(qntl = quantile(change_value))

#quantiles per condition
qc = dem %>%
  group_by(seed, condition) %>%
  subset(cycle == max(cycle)) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  summarise(qntl = quantile(change_value))

#looking at average
avg_change = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_change = mean(change_value))

#qunatiles per seed
av_qs = avg_change %>% 
  group_by(seed) %>% 
  summarise(qntl = quantile(avg_change))

#qunatiles per seed
av_qc = avg_change %>% 
  group_by(condition) %>% 
  summarise(qntl = quantile(avg_change)[4])

#create col that signals in which ntile the value is per condition
avg_q_change = 
  avg_change %>% 
  group_by(condition) %>% 
  mutate(qntl = ntile(avg_change, 4)) %>% 
  ungroup()

ggplot(data = avg_q_change)+
  geom_tile(aes(condition, seed, fill = qntl)) +
  scale_fill_gradientn("quantile",colours = rbg)

####make sure each cycle is either 125 timesteps or around 10000 inds####
c_check = dem %>% 
  select(cycle,n_timesteps,total_n,seed,condition)

#should give back empty tibble
c_check %>% 
  group_by(seed, condition,cycle) %>% 
  subset(n_timesteps != 125) %>% 
  subset(total_n < 9000)

####plot mean n_timesteps for condition seed combo####
#summarise conditions and timestep
cond_tmstps = dem %>% 
  group_by(condition,seed) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(starts_with("env_")|starts_with("n_time"),
         total_n,
         condition,
         seed) 

ggplot(data = cond_tmstps)+
  geom_tile(aes(x = condition, y = seed, fill = n_timesteps))


####check whihc environmental parameters induce####

#select only relevant conditions
plot_tmstps = cond_tmstps %>% 
  select(env_p_4,env_p_3,env_p_10,env_p_15,env_p_19, env_p_22,n_timesteps, total_n) %>%
  pivot_longer(env_p_4 : env_p_19, names_to = "variable")

ggplot(data = plot_tmstps)+
  geom_point(aes(x = value , y = total_n))+
  geom_smooth(aes(x = value, y = total_n))+
  facet_grid(.  ~ variable )

#for no_uptake check that var of uptake is 0
dem %>% 
  subset(as.numeric(seed) > 31) %>% 
  summarise(var_upt = sd(env_p_12))

dem$env_p_12 = as.factor(dem$env_p_12) 
ggplot(data = dem %>%  subset()) +
  geom_tile(aes(condition, seed, fill = env_p_12))

####subset only those simulation where pop cap is never reached####
####(always 125 timesteps)

dem_all_tmstps = dem %>% 
  subset(n_timesteps == 125 ) %>% 
  group_by(condition,seed) %>% 
  filter(n() == 499) %>% 
  ungroup()  

dem = dem_all_tmstps

#should give back empty tibble
dem_all_tmstps %>% 
  group_by(condition, seed) %>% 
  summarise(var = sd(n_timesteps), count = n()) %>% 
  subset(var != 0 & count != 499)

#plot
ggplot(data = pivot_longer(dem_all_tmstps %>% subset(as.numeric(seed) < 15 & as.numeric(condition) < 20),
                           c(spore, active, sporu, total_n),
                           names_to = "variable"))+
  geom_line(aes(x = cycle, y = value, color = variable)) +
  ylim(0,8000) +
  facet_grid(seed  ~ condition)

####subset only those simulation where pop cap is reached####
####(always 125 timesteps)

dem_cap = dem %>% 
  group_by(condition,seed) %>% 
  filter(sum(n_timesteps) != 125 * 499 ) %>% 
  ungroup()  

###n_row dem_cap +n_row dem_all_tmstps should be equal to n_row dem
nrow(dem_cap) + nrow(dem_all_tmstps) == nrow(dem)

ggplot(data = pivot_longer(dem_cap %>% subset(as.numeric(seed) < 30 & as.numeric(condition) < 25),
                           c(spore, active, sporu, total_n),
                           names_to = "variable"))+
  geom_line(aes(x = cycle, y = value, color = variable)) +
  ylim(0,11000) +
  facet_grid(seed  ~ condition)

####make sure env conditions are the same####

#make sure conditions stay the same throughout a simulation
#should return 1 row
dem %>% 
  select(starts_with("env_p") | condition | seed) %>% 
  group_by(condition, seed) %>%
  summarise_all(funs(n_distinct(.))) %>%
  ungroup() %>% 
  select(starts_with("env_p")) %>% 
  unique()

#get all different conditions existing in demo with seed for reference
demographic_conditions_cycle =
  dem %>% 
  select(starts_with("env_p") | condition | cycle | seed) %>% 
  group_by(condition) %>% 
  subset(!duplicated(env_p_3))


#Plot the different environments#
#(made to verify that one condition should have same param)
dem_c = demographic_conditions_cycle %>% 
  group_by(condition) %>% 
  mutate(env_type = seq(1:length(condition))) %>% 
  ungroup() %>% 
  select(-c(cycle,condition,seed))

dem_c$env_type = as.factor(dem_c$env_type)

dem_d = dem %>% 
  left_join(dem_c)

dem_d$env_type = as.factor(dem_d$env_type)

#for normal
ggplot(data = dem_d %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = env_type))

#for sequence
ggplot(data = dem_d) +
  geom_rect(aes(xmin=cycle,
                xmax=cycle+1,
                ymin=min(ratio_value),
                ymax=max(ratio_value), 
                fill = env_type)) +
  facet_grid(seed ~ condition) +
  theme(legend.position = "none")


ggsave("../../research presentation/env_seqex_cub.pdf",
       width = 500,
       height = 300, 
       units = "cm",
       limitsize = F)

####plot populations that do not reach cap####
ggplot(data = dem %>%  subset(cycle == max(cycle))) +
  geom_tile(aes(condition,
                seed,
                fill = n_timesteps)
  )
####select top 20 seeds for spore production####
top_20 = data.frame()

setwd(evo_dir)

for (i in  list.files(path = evo_dir,
                      pattern = "sim_demographic_s\\d+change_\\d+.*"))
{
  replicate = read.csv(i)
  replicate = replicate[,1:4]
  replicate$seed = sub( "^.*s(\\d+).*",'\\1', i)
  replicate$change = sub( "^.*change_(\\d+).*",'\\1', i)
  colnames(replicate) = colnames(top_20)
  top_20 = rbind(replicate, top_20)
}

colnames(top_20) = c("cycle",
                     "active",
                     "spore",
                     "sporulating",
                     "seed",
                     "change_freq")

top_20 = top_20 %>% 
  subset(cycle == 499) %>% 
  subset(change_freq == 0) %>% #demographic is only for change_freq == 0
  slice_max(spore, n = 20) %>% 
  select(seed) %>% 
  ungroup()

top_20$seed = as.factor(top_20$seed)
### Reduce demographic only to the best
top_dem = top_20 %>% left_join(demographic ) %>% drop_na()

####Create color plaettes####

n_colors = 100
rbg<- colorRampPalette(c("red", "blue", "green"))(n_colors)
cub_hel <- rev(cubehelix(n_colors))

####make sections for color legend in clust####

#the vector needs to be longer the the number of colors picked in the color palettes
####for overall_ratio_value
ov_qntl= quantile(dem$overall_r_value)
sections_ov_value = unique(
  c( seq(ov_qntl[1], ov_qntl[2], length = as.integer(n_colors/3 + 4)),
     seq(ov_qntl[2], ov_qntl[4], length = as.integer(n_colors/3 + 1) ),
     seq(ov_qntl[4], ov_qntl[5], length = as.integer(n_colors/3) )))

####for ratio_value cluster
rv_qntl = quantile(dem$ratio_value)
sections_rv_value = unique(
  c( seq(rv_qntl[1], rv_qntl[2], length = as.integer(n_colors/3 + 4)),
     seq(rv_qntl[2], rv_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(rv_qntl[4], rv_qntl[5], length = as.integer(n_colors/3 ))))

####for ratio_value trend plot
rv_qntl = quantile(dem$ratio_value)
sections_rv_value_plot = unique(
  c( seq(rv_qntl[1], rv_qntl[2], length = as.integer(n_colors/3 - 2)),
     seq(rv_qntl[2], rv_qntl[4], length = as.integer(n_colors/3 - 1)),
     seq(rv_qntl[4], rv_qntl[5], length = as.integer(n_colors/3 + 8))))

####for change value
c_qntl = quantile(dem$change_value)
sections_c_value = unique(
  c( seq(c_qntl[1], c_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(c_qntl[2], c_qntl[4], length = as.integer(n_colors/3 + 4)),
     seq(c_qntl[4], c_qntl[5], length = as.integer(n_colors/3))))


####for ratio_value cluster
avg_ratio = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_ratio = mean(ratio_value))
rv_avg_qntl = quantile(avg_ratio$avg_ratio)
sections_avg_rv_value = unique(
  c( seq(rv_avg_qntl[1], rv_avg_qntl[2], length = as.integer(n_colors/3 + 2)),
     seq(rv_avg_qntl[2], rv_avg_qntl[4], length = as.integer(n_colors/3 + 3)),
     seq(rv_avg_qntl[4], rv_avg_qntl[5], length = as.integer(n_colors/3 ))))

####for value
v_qntl = quantile(dem$spore)
sections_v_value = unique(
  c( seq(v_qntl[1], v_qntl[2], length = as.integer(n_colors/3 + 3)),
     seq(v_qntl[2], v_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(v_qntl[4], v_qntl[5], length = as.integer(n_colors/3))))


####Plotting spore/active/sporu lines + success on backgr####

env_change_time = dem$cycle[which(dem$env_type != dplyr::lag(dem$env_type))]
ggplot(data = dem %>% filter(as.numeric(seed) < 20, as.numeric(condition) < 20) %>% 
         mutate( spores = spore) %>%  
         pivot_longer(c(spore, sporu, active, total_n))) +
  geom_rect(aes(xmin= cycle - 1,
                xmax= cycle ,
                ymin= 0,
                ymax= max(value)  + 1,
                fill =  success), alpha =0.5) +
  scale_fill_gradientn("success",colors = rbg) +
  geom_vline(xintercept = env_change_time, linetype = "dashed", size = 0.1) +
  geom_line(aes(cycle,value, color = name)) +
  facet_grid(seed ~ condition)

ggsave("C:/Users/p288427/Desktop/research presentation/death_death_plots/death_ind+suc_20seeds_20conditions.pdf",
       width = 500,
       height = 300, 
       units = "cm",
       limitsize = F)
####Plotting spore value every last cycle before condition change####
#(useful only when change of cond)

env_duration = (max(dem$cycle) + 1) / length(levels(dem$env_type))

ggplot(data = dem  %>% subset(cycle %% env_duration == env_duration - 1)) +
  geom_rect(aes(xmin=cycle - env_duration,
                xmax=cycle,
                ymin=min(min(spore)),
                ymax=max(max(spore)), 
                fill = env_type,
  ), alpha = 0.5) +
  geom_line(aes(x = cycle, y = spore), col = "black") +
  facet_grid(seed ~ condition)

ggsave("../../research presentation/before_change_values.pdf",
       width = 500,
       height = 300, 
       units = "cm", 
       limitsize = F)
####Plotting improvement value and ratio value####

ggplot(data = dem %>% pivot_longer(c(change_value,ratio_value))) +
  geom_rect(data = dem, 
            aes(xmin=cycle,
                xmax=cycle+1,
                ymin=min(min(change_value),min(ratio_value)),
                ymax=max(max(change_value),max(ratio_value)), 
                fill = ratio_value)) +
  geom_line(aes(cycle, value, color = name)) +
  scale_fill_gradientn("Ratio",colors = rbg, breaks = sections_rv_value_plot) +
  facet_grid(seed ~ condition)

ggsave("../../research presentation/improvement+ratio_seqex_rbg.pdf",
       width = 500,
       height = 300, 
       units = "cm",
       limitsize = F)

####Plotting same as above or 3 best/worst in 3-6 most interesting envs####

#find best 3 and worst seeds
avg_ratio_pop = 
  dem %>%
  group_by(seed, condition) %>% 
  summarise(avg_ratio = mean(ratio_value)) %>% 
  ungroup() %>% 
  group_by(seed) %>% 
  summarise(avg_avg_ratio = mean(avg_ratio)) 

best = avg_avg_ratio_pop %>% 
  slice_max(avg_avg_ratio , n = 3)
worst = avg_avg_ratio_pop %>% 
  slice_min(avg_avg_ratio, n = 3)

best_worst_avg = rbind(best, worst) %>%
  mutate(seed = factor(seed)) %>%
  left_join(dem)

ggplot(data = best_worst_avg) +
  geom_rect(aes(xmin=cycle,
                xmax=cycle + 1,
                ymin=min(min(spore)),
                ymax=max(max(spore)), 
                fill = ratio_value,
  ), alpha = 0.5) +
  scale_fill_gradientn("Ratio",colors = rbg, breaks = sections_rv_value_plot) +
  geom_line(aes(x = cycle, y = spore), col = "pink", size = 0.2) +
  geom_rect(data = best_worst_avg %>%
              subset(cycle %% env_duration == env_duration - 1),
            aes(xmin=cycle - env_duration,
                xmax=cycle,
                ymin=min(min(spore)),
                ymax=max(max(spore))),
            size = 0.1,
            fill = NA,
            colour = "black") +
  facet_grid(seed ~ condition)

ggsave("../../research presentation/best_worst_seqex_rbg.pdf",
       width = 750,
       height = 75, 
       units = "cm",
       limitsize = F)

#find most mutable environments
sd_avg_ratio_cond = 
  dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_ratio = mean(ratio_value)) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  summarise(sd_avg_ratio = sd(avg_ratio))

best_sd = sd_avg_ratio_cond %>% 
  slice_max(sd_avg_ratio, n = 3)
worst_sd = sd_avg_ratio_cond %>% 
  slice_min(sd_avg_ratio, n = 3)

####Plot final change value####
ggplot(data = dem ) +
  geom_tile(aes(seed, condition, fill = change_value)) +
  scale_fill_gradientn("Final change value", colors = cubehelix(20))


####Plot quantiles of avg change####
avg_q_change = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_change = mean(change_value))%>% 
  group_by(condition) %>% 
  mutate(qntl = ntile(avg_change, 5)) %>% 
  ungroup()

ggplot(data = avg_q_change)+
  geom_tile(aes(condition, seed, fill = qntl)) +
  scale_fill_gradientn("quantile",colours = rbg)

####Plotting diff between rvalue end vs start normalized####


ggplot(data = dem%>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition,
                seed,
                fill = ratio_end_production - ratio_start_production),
            color = "black", size = 0.5) + 
  scale_fill_gradientn("Delta", colors = rev(cubehelix(10))) 
ggsave("../research presentation/delta.png", width = 50, height = 30, units = "cm")

###Plotting start or end value####
###ratio_value

ggplot(data = dem %>% subset(cycle == min(cycle))) +
  geom_tile(aes(condition, seed, fill = ratio_value),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))
ggsave("../research presentation/start_rv", device = "png",
       width = 50, height = 30, units = "cm")

ggplot(data = dem %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = ratio_value),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))
ggsave("../research presentation/end_rv",device = "png",
       width = 50, height = 30, units = "cm")



###value
ggplot(data = demographic %>% subset(cycle == min(cycle))) +
  geom_tile(aes(condition, seed, fill = spore),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))

ggsave("../research presentation/start_v",device = "png",
       width = 50, height = 30, units = "cm")

ggplot(data = demographic %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = spore),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))

ggsave("../research presentation/end_v",device = "png",
       width = 50, height = 30, units = "cm")

####heatmap cluster for ratio_production of spores####

#start
clust_dem_start = dem %>% 
  subset(cycle == min(cycle)) %>% 
  select(seed, ratio_value, condition) %>% 
  spread(condition, ratio_value) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

#col and row clustering at the start
start_clust_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_start)))))
start_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_start)))))

#end
clust_dem_end = dem %>% 
  # subset(as.numeric(seed) > 30) %>% 
  subset(cycle == max(cycle)) %>% 
  select(seed, ratio_value, condition) %>% 
  spread(condition, ratio_value) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

#col and row clustering at the end
end_clust_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_end)))))
end_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_end)))))

#clust_Start
#plot start with clust_start
heatmap.2(as.matrix(clust_dem_start) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = rbg, 
          breaks = sections_rv_value,
          main = "start_c_start")

#plot end with clust_start
heatmap.2(as.matrix(clust_dem_end) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = rbg, 
          breaks = sections_rv_value,
          main = "end_c_start")

#clust_end
#plot start with clust_end
heatmap.2(as.matrix(clust_dem_start) ,
          Rowv = end_clust_row,
          Colv = end_clust_col,
          scale = "none",
          trace = "none",
          col = rbg, 
          breaks = sections_rv_value,
          main = "start_c_end")

#plot end with clust_end
heatmap.2(as.matrix(clust_dem_end) ,
          Rowv = end_clust_row,
          Colv = end_clust_col,
          scale = "none",
          trace = "none",
          col = rbg, 
          breaks = sections_rv_value,
          main = "end_c_end")

#plot in loop for all cycle and save images
for(cycle_n in min(dem$cycle) : max(dem$cycle))
{
  clust_dem = dem %>% 
    subset(cycle == cycle_n) %>% 
    select(seed, ratio_value, condition) %>% 
    spread(condition, ratio_value) %>% 
    column_to_rownames(var = "seed") %>% 
    drop_na()
  
  
  jpeg(paste("Clust_rv_cycle_all_tmstps",cycle_n,".png", sep = "_")
       ,width = 700,
       height = 700)
  
  #by keeping clustering structure fixed
  heatmap.2(as.matrix(clust_dem) ,
            Rowv = start_clust_row,
            Colv = start_clust_col,
            scale = "none",
            trace = "none",
            col = rbg, 
            breaks = sections_rv_value,
            main = cycle_n)
  
  dev.off() 
}


####heatmap cluster for avg ratio value####
#avg rv
clust_dem_rv_avg = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_ratio = mean(ratio_value))%>% 
  spread(condition, avg_ratio) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

#col clustering for avg value
avg_rv_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_rv_avg)))))
#row clustering for avg ratio_value
avg_rv_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_rv_avg)))))

####for ratio_value cluster
avg_ratio = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_ratio = mean(ratio_value))
rv_avg_qntl = quantile(avg_ratio$avg_ratio)
sections_avg_rv_value = unique(
  c( seq(rv_avg_qntl[1], rv_avg_qntl[2], length = as.integer(n_colors/3 + 2)),
     seq(rv_avg_qntl[2], rv_avg_qntl[4], length = as.integer(n_colors/3 + 2)),
     seq(rv_avg_qntl[4], rv_avg_qntl[5], length = as.integer(n_colors/3 ))))

#clust avg
heatmap.2(as.matrix(clust_dem_rv_avg),        
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          scale = "none",
          trace = "none",
          col = rbg, 
          breaks = sections_avg_rv_value,
          main = "avg_ratio_value")

####heatmap cluster for avg spore value####
avg_spore_v = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_spore = mean(spore))

clust_avg_spore = avg_spore_v %>% 
  spread(condition, avg_spore) %>% 
  column_to_rownames("seed")

#col and row clustering at the start
avg_spore_col = as.dendrogram(hclust(dist(t(as.matrix(clust_avg_spore)))))
avg_spore_row = as.dendrogram(hclust(dist((as.matrix(clust_avg_spore)))))

heatmap.2(as.matrix(clust_avg_spore),
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          scale = "none",
          trace = "none",
          main = "avg_spore_value")

####heatmap cluster for avg success value####
avg_success_v = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_success = mean(success))

clust_avg_success = avg_success_v %>% 
  spread(condition, avg_success) %>% 
  column_to_rownames("seed")

#col and row clustering at the start
avg_success_col = as.dendrogram(hclust(dist(t(as.matrix(clust_avg_success)))))
avg_success_row = as.dendrogram(hclust(dist((as.matrix(clust_avg_success)))))

heatmap.2(as.matrix(clust_avg_success),
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          scale = "none",
          trace = "none",
          main = "avg_success_value",
          breaks = sections_avg_start_sc_value_plot)

####heatmap cluster for overall_r_value(ratio_value rescaled to [0-1] over all cycle) of spores####
#start
clust_dem_ov_start = dem %>% 
  subset(cycle == min(cycle)) %>% 
  select(seed, overall_r_value, condition) %>% 
  spread(condition, overall_r_value) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

#col and row clustering at the start
start_clust_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_ov_start)))))
start_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_ov_start)))))

#end
clust_dem_ov_end = dem %>% 
  subset(cycle == max(cycle)) %>% 
  select(seed, overall_r_value, condition) %>% 
  spread(condition, overall_r_value) %>% 
  column_to_rownames(var = "seed")%>% 
  drop_na()

#col and row clustering at the end
end_clust_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_ov_end)))))
end_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_ov_end)))))


#clust_Start
#plot start with clust_start
heatmap.2(as.matrix(clust_dem_ov_start) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = rbg,
          breaks = sections_ov_value)



#plot end with clust_start
heatmap.2(as.matrix(clust_dem_ov_end) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = rbg,
          main = "end_c_start",
          breaks = sections_ov_value)

#clust_end
#plot start with clust_end
heatmap.2(as.matrix(clust_dem_ov_start) ,
          Rowv = end_clust_row,
          Colv = end_clust_col,
          scale = "none",
          trace = "none",
          col = cub_hel,
          main = "start_c_end",
          breaks = sections_ov_value)

#plot end with clust_end
heatmap.2(as.matrix(clust_dem_ov_end) ,
          Rowv = end_clust_row,
          Colv = end_clust_col,
          scale = "none",
          trace = "none",
          col = cub_hel,
          main = "end_c_end",
          breaks = sections_ov_value)

#plot in loop for all cycle and save images
for(cycle_n in min(dem$cycle) : max(dem$cycle))
{
  clust_dem = dem %>% 
    subset(cycle == cycle_n) %>% 
    select(seed, overall_r_value, condition) %>% 
    spread(condition, overall_r_value) %>% 
    column_to_rownames(var = "seed") %>% 
    drop_na()
  
  
  jpeg(paste("Clust_ov_start_all_tmstps",cycle_n,".png", sep = "_")
       ,width = 700,
       height = 700)
  
  #by keeping clustering structure fixed
  heatmap.2(as.matrix(clust_dem) ,
            Rowv = start_clust_row,
            Colv = start_clust_col,
            scale = "none",
            trace = "none",
            col = rbg, 
            main = cycle_n,
            breaks = sections_ov_value)
  
  dev.off() 
}


####heatmap cluster for value of spores####
#start
clust_dem_v_start = dem %>% 
  subset(cycle == min(cycle)) %>% 
  select(seed, spore, condition) %>% 
  spread(condition, spore) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

#col and row clustering at the start
start_clust_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_v_start)))))
start_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_v_start)))))

#end
clust_dem_v_end = dem %>% 
  subset(cycle == max(cycle)) %>% 
  select(seed, spore, condition) %>% 
  spread(condition, spore) %>% 
  column_to_rownames(var = "seed")%>% 
  drop_na()

#col and row clustering at the end
end_clust_col = as.dendrogram(hclust(dist(t(as.matrix(clust_dem_v_end)))))
end_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_dem_v_end)))))


#clust_Start
#plot start with clust_start
heatmap.2(as.matrix(clust_dem_v_start) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = rbg,
          main = "start_c_start",
          breaks = sections_v_value)



#plot end with clust_start
heatmap.2(as.matrix(clust_dem_v_end) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = rbg,
          main = "end_c_start",
          breaks = sections_v_value)

#clust_end
#plot start with clust_end
heatmap.2(as.matrix(clust_dem_v_start) ,
          Rowv = end_clust_row,
          Colv = end_clust_col,
          scale = "none",
          trace = "none",
          col = cub_hel,
          main = "start_c_end",
          breaks = sections_v_value)

#plot end with clust_end
heatmap.2(as.matrix(clust_dem_v_end) ,
          Rowv = end_clust_row,
          Colv = end_clust_col,
          scale = "none",
          trace = "none",
          col = cub_hel,
          main = "end_c_end",
          breaks = sections_v_value)

#plot in loop for all cycle and save images
for(cycle_n in min(dem$cycle) : max(dem$cycle))
{
  clust_dem = dem %>% 
    subset(cycle == cycle_n) %>% 
    select(seed, spore, condition) %>% 
    spread(condition, spore) %>% 
    column_to_rownames(var = "seed") %>% 
    drop_na()
  
  
  jpeg(paste("Clust_v_all_tmstps",cycle_n,".png", sep = "_")
       ,width = 700,
       height = 700)
  
  #by keeping clustering structure fixed
  heatmap.2(as.matrix(clust_dem) ,
            Rowv = start_clust_row,
            Colv = start_clust_col,
            scale = "none",
            trace = "none",
            col = rbg, 
            main = cycle_n,
            breaks = sections_v_value)
  
  dev.off() 
}

####heatmap cluster for delta in ranking####
clust_delta = dem %>% 
  subset(cycle == max(cycle)) %>% 
  select(seed, delta_rv_start_end, condition) %>% 
  spread(condition, delta_rv_start_end) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()


heatmap.2(as.matrix(clust_delta) , scale = "none",
          col = rbg)

#plot it next to start ratio_value
par(mfrow=c(1,2))
layout(matrix(c(1,2), 1, 2, byrow = FALSE))

#by keeping clustering structure fixed
heatmap.2(as.matrix(clust_dem_start) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          col = cub_hel, 
          main = "start_ratio_value")
#by keeping clustering structure fixed
heatmap.2(as.matrix(clust_delta) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          col = cub_hel, 
          main = "delta_ratio_value")

####heatmap cluster for avg change value#####
clust_change = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_change = mean(change_value))%>% 
  ungroup()  %>% 
  select(condition, seed, avg_change) %>% 
  spread(condition, avg_change) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

heatmap.2(as.matrix(clust_change) ,
          scale = "none",
          trace = "none",
          col = rbg,
          breaks = sections_c_value)


####heatmap cluster for quantiles of avg change value#####

clust_change_qntl = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_change = mean(change_value))%>% 
  ungroup() %>% 
  group_by(condition) %>% 
  mutate(qntl = ntile(avg_change, 5)) %>% 
  ungroup()  %>% 
  select(condition, seed, qntl) %>% 
  spread(condition, qntl) %>% 
  column_to_rownames(var = "seed") %>% 
  drop_na()

heatmap.2(as.matrix(clust_change_qntl) ,
          scale = "none",
          trace = "none",
          col = rbg)

####Heatmap cluster for slope of spore production taking last cycle before change for each env_type#####
env_duration = (max(dem$cycle) + 1) / length(levels(dem$env_type))

clust_slope = dem %>% 
  subset(cycle %% env_duration == env_duration - 1) %>% 
  group_by(condition, seed) %>%
  do(coeff = coefficients(lm(spore ~ cycle, data = .))[2][[1]]) %>% 
  mutate(coeff = as.numeric(coeff)) %>% 
  ungroup() %>% 
  spread(condition, coeff) %>% 
  column_to_rownames(var = "seed") 


#create breaks for bef_ch_slope
bef_ch_sl_qntl = quantile(unlist(clust_slope))
sections_bf_ch_sl_value_plot = unique(
  c( seq(bef_ch_sl_qntl[1], bef_ch_sl_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(bef_ch_sl_qntl[2], bef_ch_sl_qntl[4], length = as.integer(n_colors/3 + 2)),
     seq(bef_ch_sl_qntl[4], bef_ch_sl_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_slope),
          scale = "none",
          trace = "none",
          col = rbg,
          breaks = sections_bf_ch_sl_value_plot,
          main = "v_before_change_slope")


####Heatmap cluster for avg of SLOPE of SPORE production in each env_type#####


#summarise average slope over all env_types in a condition
avg_slope = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_slope = mean(slope))

#create dataframe for slope
clust_avg_slope = avg_slope%>% 
  select(condition, seed, avg_slope) %>% 
  spread(condition, avg_slope) %>% 
  column_to_rownames(var = "seed") 

#pre-cluster rows so to use it in clustering of intercept
avg_sl_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_avg_slope)))))

#create breaks for avg_slope
avg_sl_qntl = quantile(avg_slope$avg_slope)
sections_avg_sl_value_plot = unique(
  c( seq(avg_sl_qntl[1], avg_sl_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(avg_sl_qntl[2], avg_sl_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(avg_sl_qntl[4], avg_sl_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_avg_slope),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_avg_sl_value_plot,
          main = "avg_slope")
####Heatmap cluster for avg of INTERCEPT of SPORE production in each env_type#####

#summarise average intercept over all env_types in a condition
avg_intercept = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_intercept = mean(intercept))

#create dataframe for intercept
clust_avg_intercept = avg_intercept%>% 
  select(condition, seed, avg_intercept) %>% 
  spread(condition, avg_intercept) %>% 
  column_to_rownames(var = "seed") 

#create breaks for avg_intercept
avg_int_qntl = quantile(avg_intercept$avg_intercept)
sections_avg_int_value_plot = unique(
  c( seq(avg_int_qntl[1], avg_int_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(avg_int_qntl[2], avg_int_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(avg_int_qntl[4], avg_int_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_avg_intercept),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_avg_int_value_plot,
          main = "avg_intercept")


####Heatmap cluster for avg of STARTING SPORE production in each env_type#####

#summarise average initial spore quantity over all env_types in a condition
avg_start_sp = dem %>% 
  select(spore, cycle, condition, seed, env_type ) %>% 
  group_by(condition, seed, env_type) %>%
  slice(which.min(cycle)) %>% 
  ungroup() %>% 
  group_by(condition, seed) %>% 
  summarise(avg_start_sp = mean(spore))

#create dataframe for intercept
clust_avg_start_sp = avg_start_sp%>% 
  select(condition, seed, avg_start_sp) %>% 
  spread(condition, avg_start_sp) %>% 
  column_to_rownames(var = "seed") 

#create breaks for avg_intercept
avg_start_sp_qntl = quantile(avg_start_sp$avg_start_sp)
sections_avg_start_sp_value_plot = unique(
  c( seq(avg_start_sp_qntl[1], avg_start_sp_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(avg_start_sp_qntl[2], avg_start_sp_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(avg_start_sp_qntl[4], avg_start_sp_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_avg_start_sp),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_avg_start_sp_value_plot,
          main = "avg_start_sp")


####Heatmap cluster for the SLOPE OF SLOPE for SPORE PROD for each env type in each seq####

env_duration = (max(dem$cycle) + 1) / length(levels(dem$env_type))

###slope of slopes
slope_of_slopes = dem %>%
  subset(cycle %% env_duration == env_duration - 1) %>% 
  group_by(condition, seed) %>%
  do(coeffs = coefficients(lm(slope ~ cycle, data = .))) %>% 
  mutate(slope_of_slope = as.numeric(coeffs[2][[1]])) %>% 
  ungroup()

clust_slope_slope = slope_of_slopes %>% 
  select(condition, seed, slope_of_slope) %>% 
  spread(condition, slope_of_slope) %>% 
  column_to_rownames(var = "seed")

#pre-cluster rows so to use it in clustering of intercept
sl_sl_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_slope_slope)))))

#create breaks for slope_of_slopes
sl_sl_qntl = quantile(slope_of_slopes$slope_of_slope)
sections_sl_sl_value_plot = unique(
  c( seq(sl_sl_qntl[1], sl_sl_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(sl_sl_qntl[2], sl_sl_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(sl_sl_qntl[4], sl_sl_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_slope_slope),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          #breaks = sections_sl_sl_value_plot,
          main = "slope_of_slopes")


####heatmap cluster for SLOPE OF INTERCEPT  for SPORE PROD ####

env_duration = (max(dem$cycle) + 1) / length(levels(dem$env_type))

slope_of_intercept = dem %>%
  subset(cycle %% env_duration == env_duration - 1) %>% 
  group_by(condition, seed) %>%
  do(coeffs = coefficients(lm(intercept ~ cycle, data = .))) %>% 
  mutate(slope_of_intercept = as.numeric(coeffs[2][[1]])) %>% 
  ungroup()

clust_slope_intercept = slope_of_intercept %>% 
  select(condition, seed, slope_of_intercept) %>% 
  spread(condition, slope_of_intercept) %>% 
  column_to_rownames(var = "seed")

#create breaks for slope_of_intercept
sl_int_qntl = quantile(slope_of_intercept$slope_of_intercept)
sections_sl_int_value_plot = unique(
  c( seq(sl_int_qntl[1], sl_int_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(sl_int_qntl[2], sl_int_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(sl_int_qntl[4], sl_int_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_slope_intercept),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_sl_int_value_plot,
          main = "slope_of_intercepts")
####Heatmap cluster for avg of SLOPE of SUCCESS in each env_type#####

#summarise average slope over all env_types in a condition
avg_slope_success = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_slope_success = mean(slope_success))

#create dataframe for slope
clust_avg_slope_success = avg_slope_success%>% 
  select(condition, seed, avg_slope_success) %>% 
  spread(condition, avg_slope_success) %>% 
  column_to_rownames(var = "seed") 

#pre-cluster rows so to use it in clustering of intercept
avg_sl_clust_row_scs = as.dendrogram(hclust(dist((as.matrix(clust_avg_slope_success)))))

#create breaks for avg_slope
avg_sl_scs_qntl = quantile(avg_slope_success$avg_slope_success)
sections_avg_sl_scs_value_plot = unique(
  c( seq(avg_sl_scs_qntl[1], avg_sl_scs_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(avg_sl_scs_qntl[2], avg_sl_scs_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(avg_sl_scs_qntl[4], avg_sl_scs_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_avg_slope_success),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_avg_sl_scs_value_plot,
          main = "avg_slope_scs")

####Heatmap cluster for avg of INTERCEPT of SUCCESS in each env_type#####

#summarise average intercept over all env_types in a condition
avg_intercept_success = dem %>% 
  group_by(seed, condition) %>% 
  summarise(avg_intercept_success = mean(intercept_success))

#create dataframe for intercept
clust_avg_intercept_success = avg_intercept_success%>% 
  select(condition, seed, avg_intercept_success) %>% 
  spread(condition, avg_intercept_success) %>% 
  column_to_rownames(var = "seed") 

#create breaks for avg_intercept
avg_int_scs_qntl = quantile(avg_intercept_success$avg_intercept_success)
sections_avg_int_scs_value_plot = unique(
  c( seq(avg_int_scs_qntl[1], avg_int_scs_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(avg_int_scs_qntl[2], avg_int_scs_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(avg_int_scs_qntl[4], avg_int_scs_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_avg_intercept_success),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_avg_int_scs_value_plot,
          main = "avg_intercept_success")

####Heatmap cluster for avg of STARTING SUCCESS production in each env_type#####

#summarise average initial spore quantity over all env_types in a condition
avg_start_sc = dem %>% 
  select(success, cycle, condition, seed, env_type ) %>% 
  group_by(condition, seed, env_type) %>%
  slice(which.min(cycle)) %>% 
  ungroup() %>% 
  group_by(condition, seed) %>% 
  summarise(avg_start_sc = mean(success))

#create dataframe for intercept
clust_avg_start_sc = avg_start_sc%>% 
  select(condition, seed, avg_start_sc) %>% 
  spread(condition, avg_start_sc) %>% 
  column_to_rownames(var = "seed") 

#create breaks for avg_intercept
avg_start_sc_qntl = quantile(avg_start_sc$avg_start_sc)
sections_avg_start_sc_value_plot = unique(
  c( seq(avg_start_sc_qntl[1], avg_start_sc_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(avg_start_sc_qntl[2], avg_start_sc_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(avg_start_sc_qntl[4], avg_start_sc_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_avg_start_sc),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_avg_start_sc_value_plot,
          main = "avg_start_sc")


####Heatmap cluster for the SLOPE OF SLOPE for SUCCESS for each env type in each seq####

env_duration = (max(dem$cycle) + 1) / length(levels(dem$env_type))

###slope of slopes
slope_of_slopes_success = dem %>%
  subset(cycle %% env_duration == env_duration - 1) %>% 
  group_by(condition, seed) %>%
  do(coeffs = coefficients(lm(slope_success ~ cycle, data = .))) %>% 
  mutate(slope_of_slope_success = as.numeric(coeffs[2][[1]])) %>% 
  ungroup()

clust_slope_slope_success = slope_of_slopes_success %>% 
  select(condition, seed, slope_of_slope_success) %>% 
  spread(condition, slope_of_slope_success) %>% 
  column_to_rownames(var = "seed")

#pre-cluster rows so to use it in clustering of intercept
sl_sl_clust_row = as.dendrogram(hclust(dist((as.matrix(clust_slope_slope_success)))))

#create breaks for slope_of_slopes
sl_sl_qntl_scs = quantile(slope_of_slopes_success$slope_of_slope_success)
sections_sl_sl_scs_value_plot = unique(
  c( seq(sl_sl_qntl_scs[1], sl_sl_qntl_scs[2], length = as.integer(n_colors/3 + 1)),
     seq(sl_sl_qntl_scs[2], sl_sl_qntl_scs[4], length = as.integer(n_colors/3 + 1)),
     seq(sl_sl_qntl_scs[4], sl_sl_qntl_scs[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_slope_slope_success),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_sl_sl_scs_value_plot,
          main = "slope_of_slopes_success")

#plot slopes of spore production slope in each env_type over the entire seq 
ggplot(data = dem %>%
         left_join(seq_slopes) %>%
         subset(cycle %% env_duration == env_duration - 1)) +
  geom_point(aes(x = cycle, y = slope)) +
  facet_grid(seed ~ condition)

ggsave("../../research presentation/slope_values.pdf",
       width = 500,
       height = 300, 
       units = "cm",
       limitsize = F)
####Heatmap cluster for SLOPE OF INTERCEPT  for SUCCESS ####

env_duration = (max(dem$cycle) + 1) / length(levels(dem$env_type))

slope_of_intercept_scs = dem %>%
  subset(cycle %% env_duration == env_duration - 1) %>% 
  group_by(condition, seed) %>%
  do(coeffs = coefficients(lm(intercept_success ~ cycle, data = .))) %>% 
  mutate(slope_of_intercept_scs = as.numeric(coeffs[2][[1]])) %>% 
  ungroup()

clust_slope_intercept_scs = slope_of_intercept_scs %>% 
  select(condition, seed, slope_of_intercept_scs) %>% 
  spread(condition, slope_of_intercept_scs) %>% 
  column_to_rownames(var = "seed")

#create breaks for slope_of_intercept
sl_int_qntl_scs = quantile(slope_of_intercept_scs$slope_of_intercept_scs)
sections_sl_int_scs_value_plot = unique(
  c( seq(sl_int_qntl_scs[1], sl_int_qntl_scs[2], length = as.integer(n_colors/3 + 1)),
     seq(sl_int_qntl_scs[2], sl_int_qntl_scs[4], length = as.integer(n_colors/3 + 1)),
     seq(sl_int_qntl_scs[4], sl_int_qntl_scs[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_slope_intercept_scs),
          scale = "none",
          trace = "none",
          Colv = start_clust_col,
          Rowv = avg_sl_clust_row,
          col = rbg,
          breaks = sections_sl_int_scs_value_plot,
          main = "slope_of_intercepts_success")

####Heatmap cluster for total sum of spore deltas from first to lst cycle in one env_type####  
total_improvement = dem %>%
  group_by(condition, seed, env_type) %>% 
  summarise(improvement = spore[length(spore)] - spore[1]) %>% 
  summarise(total_impr = sum(improvement))

ggplot(data =  total_improvement) +
  geom_tile(aes(x = seed, y = condition, fill = total_impr)) +
  scale_fill_gradientn("tot_impr",colors = rbg)

clust_tot_impr = total_improvement %>%
  group_by(condition, seed, env_type) %>% 
  summarise(improvement = spore[length(spore)] - spore[1]) %>% 
  summarise(total_impr = sum(improvement)) %>% 
  spread(condition, total_impr) %>% 
  column_to_rownames(var = "seed")

#create breaks for tot_impr
tot_impr_qntl = quantile(total_improvement$total_impr)
sections_tot_impr_value_plot = unique(
  c( seq(tot_impr_qntl[1], tot_impr_qntl[2], length = as.integer(n_colors/3 + 1)),
     seq(tot_impr_qntl[2], tot_impr_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(tot_impr_qntl[4], tot_impr_qntl[5], length = as.integer(n_colors/3 + 2))))

heatmap.2(as.matrix(clust_tot_impr),
          scale = "none",
          trace = "none",
          col = rbg,
          breaks = sections_tot_impr_value_plot,
          main = "total_impr")

####Plot the avg slope of spore production over diff env_types in a seq####
seq_slopes %>% 
  group_by(seed, condition) %>% 
  summarise(avg_slope = mean(slope)) %>% 
  ggplot() +
  geom_tile(aes(seed, condition, fill = as.numeric(avg_slope))) +
  scale_fill_gradientn("Ratio",colors = rbg) 

####Plotting variance of production over time for conditions####
var_demo = dem %>% 
  ungroup() %>% 
  group_by(cycle, condition) %>% 
  summarise(val = sd(value), rat = sd(ratio_value))

ggplot(data = var_demo) +
  geom_smooth(aes(cycle,val)) +
  facet_grid( . ~ condition)

ggplot(data = var_demo) +
  geom_smooth(aes(cycle,rat)) +
  facet_grid( . ~ condition)

####Plotting CORRELATION of avg_slope over success#####

avg_vals = tot_dem %>% 
  group_by(seed, condition, type) %>% 
  summarise(avg_slope = mean(slope_success),
            avg_intercept = mean(intercept_success),
            avg_success = mean(success))

avg_start_suc = tot_dem %>% 
  group_by(type, condition, seed, env_type) %>%
  slice(which.min(cycle)) %>% 
  ungroup() %>% 
  group_by(type, condition, seed) %>% 
  summarise(avg_start_suc = mean(success))

cor_dem =  avg_vals %>%  
  left_join(avg_start_suc) %>% 
  pivot_longer(c(avg_slope, avg_intercept, avg_start_suc))

ggplot(cor_dem %>% 
         subset(type != "eden") %>% 
         pivot_wider(names_from = name,values_from = value) %>%
         group_by(type) %>% 
         mutate(across(c(avg_slope, avg_intercept, avg_success, avg_start_suc), scale)) %>% 
         pivot_longer(c(avg_success, avg_slope, avg_intercept, avg_start_suc)))+
  geom_tile(aes(x = condition, y = seed,  fill = value ))+
  scale_fill_gradientn(colors = rbg) +
  facet_grid(type ~ name)

cdd = cor_dem %>% 
  pivot_wider(names_from = name,values_from = value) %>%
  mutate(across(c(avg_slope, avg_intercept, avg_success, avg_start_suc), scale))

ggscatter(cdd %>%
            # subset(seed != "9") %>% #taking out the outlier
            pivot_longer(c(avg_slope, avg_start_suc, avg_intercept))
          # %>% 
          # subset(name != "avg_intercept") %>% 
          # subset(name != "avg_start_suc")
          ,
          x = "value", 
          y = "avg_success", 
          add = "reg.line",
          conf.int = TRUE, 
          cor.coef = TRUE,
          cor.method = "pearson",
          main = "pearson") +
  facet_grid(name ~ type)
####Plotting CORRELATION of avg_improvment over success#####

tot_dem_deltas_and_starts = tot_dem %>%
  select(cycle, env_type, seed, condition, type, success) %>% 
  group_by(condition, seed, type, env_type) %>% 
  mutate(improvement = success[length(success)] - success[1]) %>% 
  mutate(start = success[1]) %>% 
  group_by(condition, seed, type) %>% 
  summarise(avg_start = mean(start), 
            avg_improvement = mean(improvement),
            avg_success = mean(success)) 

my.formula <- y ~ x
data = tot_dem_deltas_and_starts %>% filter(type == "death") %>% 
  pivot_longer(c(avg_start, avg_improvement)) %>%
  filter(name == "avg_improvement")
ggscatter(data = data,
          x = "value", 
          y = "avg_success",
          color = "type",
          add = "reg.line",
          add.params = list(size = 0.5),
          fullrange = T,
          size = 0.01,
          conf.int = TRUE, 
          cor.coef = TRUE,
          cor.method = "pearson",
          # main = "pearson"
          ) +
  stat_poly_eq( data = data,
                formula = my.formula, 
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE,
                label.y = "bottom",
                label.x = "right") +    
  theme(legend.position = "none")+
  facet_grid(. ~ .)

####Testing if improvement grows over time###
tot_improv_start_regr = 
  tot_dem %>%
  select(cycle, env_type, seed, condition, type, success) %>% 
  group_by(condition, seed, type, env_type) %>% 
  mutate(improvement = success[length(success)] - success[1]) %>% 
  mutate(start = success[1]) %>% 
  slice_min(cycle)

coeffs_impr_start_regr = 
  tot_improv_start_regr%>% 
  group_by(seed,condition, type) %>% 
  do(coeffs_success = coefficients(lm(improvement ~ cycle, data = .))) %>% 
  mutate("slope_improvement" = as.numeric(coeffs_success[2][[1]]))

ggplot(coeffs_impr_start_regr %>%  subset(type == "death")) + 
  geom_tile(aes(x = condition, y = seed,
                fill = slope_improvement))+
  scale_fill_gradientn("slope of improvement",colors = rbg) +
  facet_grid(. ~ type )

ggscatter(data = tot_improv_start_regr %>%  subset(type == "death"),
          x = "cycle", 
          y = "improvement",
          color = "type",
          add = "reg.line",
          add.params = list(size = 0.5),
          fullrange = T,
          size = 0.01,
          conf.int = TRUE, 
          cor.coef = TRUE,
          cor.method = "pearson",
          # main = "pearson"
) +
  stat_poly_eq( data = tot_improv_start_regr %>%  subset(type == "death"),
                formula = my.formula, 
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE,
                label.y = "bottom",
                label.x = "right") +    
  theme(legend.position = "none")+
  facet_grid(condition ~ seed)
###Testing if slope of improv correlates with final success
ggscatter(tot_dem_deltas_and_starts %>%
            left_join(coeffs_impr_start_regr %>%
                        select(-coeffs_success)),
          x = "", 
          y = "avg_success", 
          add = "reg.line",
          conf.int = TRUE, 
          cor.coef = TRUE,
          cor.method = "pearson",
          main = "pearson") +
  facet_grid(name ~ type)

####Testing if SUCCESS_SLOPE changes over time(evolves)####
slope_of_slopes = dem %>%
  group_by(condition, seed, env_type) %>% 
  slice_min(cycle) %>%
  group_by(condition, seed) %>%
  do(model = lm(slope_success ~ cycle, data = .))
unserialize(slope_of_slopes$model)
unlist(slope_of_slopes$model)