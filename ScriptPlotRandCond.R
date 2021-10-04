####Setting up####
library(broom)
library(dplyr)
library(gplots)
library(ggplot2)
library("ggpubr")
library(ggpmisc)
library(ggh4x)
library(pals)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)
library(ggstatsplot)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
first_s = "C:/Users/p288427/Desktop/hd_rand_evo/test_first"
last_s = "C:/Users/p288427/Desktop/hd_rand_evo/test_last"
first_death = "C:/Users/p288427/Desktop/hd_rand_evo/test_first_death"
last_death = "C:/Users/p288427/Desktop/hd_rand_evo/test_last_death"
first_eden = "C:/Users/p288427/Desktop/hd_rand_evo/test_first_eden"
last_eden = "C:/Users/p288427/Desktop/hd_rand_evo/test_last_eden"
first_sel = "C:/Users/p288427/Desktop/hd_rand_evo/test_first_sel"
last_sel = "C:/Users/p288427/Desktop/hd_rand_evo/test_last_sel"
folders = c(first_s, last_s, 
            first_death, last_death, 
            first_eden, last_eden,
            first_sel, last_sel)
types = c("normal","eden", "death", "sel")
# types = c("death","normal","sel")

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
  else if(getwd() == first_sel){
    
    pattern = "first_gen_sel_sp_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
  } else if(getwd() == last_sel){
    
    pattern = "last_gen_sel_sp_rand_evo_extreme_a\\d+.000000seq_\\d+cond_per_seq\\d+random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+"
    
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
                                "cond_per_seq",
                                "seq")
  
  rand_demographic$seed = as.numeric(rand_demographic$seed)
  rand_demographic$change_freq = as.factor(rand_demographic$change_freq)
  rand_demographic$condition = as.factor(rand_demographic$condition)
  rand_demographic$n_timesteps = as.numeric(rand_demographic$n_timesteps)
  
  
  rand_demographic =  rand_demographic %>% 
    group_by(seed, condition, seq) %>% 
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
    
  }else if (getwd() == first_sel){
    
    first_sel_rand_demographic = rand_demographic
    save(first_sel_rand_demographic, file = "first_sel_test_rand_evo_demo.R")
    
  } else if (getwd() == last_sel){
    
    last_sel_rand_demographic = rand_demographic
    save(last_sel_rand_demographic, file = "last_sel_test_rand_evo_demo.R")
    
  }
  print(j)
  
}
####loading separate dataframe####

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
    
  }else if (getwd() == first_sel){
    
    load("first_sel_test_rand_evo_demo.R")
    
  } else if (getwd() == last_sel){
    
    load("last_sel_test_rand_evo_demo.R")
    
  }
}

#### Unite dataframes ####
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
  } else if ( i == "sel")
  {
    first = first_sel_rand_demographic
    last = last_sel_rand_demographic
  }
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
  } else if (i == "sel"){
    
    setwd(first_sel)
    sel_test_demog = test_demog
    save(sel_test_demog, file = "sel_test_demog.R") 
  }
}

eden_test_demog$type = "eden"
a = eden_test_demog %>% 
  select(gen, seed, success, condition, type, seq) %>% 
  pivot_wider(names_from = gen, values_from = success)

death_test_demog$type = "death"
b = death_test_demog %>% 
  select(gen, seed, success, condition, type, seq) %>% 
  pivot_wider(names_from = gen, values_from = success)

seqex_test_demog$type = "seqex"
c = seqex_test_demog %>% 
  select(gen, seed, success, condition, type, seq) %>% 
  pivot_wider(names_from = gen, values_from = success)

sel_test_demog$type = "sel"
d = sel_test_demog %>% 
  select(gen, seed, success, condition, type, seq) %>% 
  pivot_wider(names_from = gen, values_from = success)

# total = rbind(a, b, c) %>% 
#   subset(seed != "9") %>%
#   group_by(seed, type, seq) %>% 
#   summarise(avg_first = mean(first), avg_last = mean(last)) 
# levels(total$type) = c("seqex","eden","death")
#   
total = rbind(a, b, c, d) %>% 
  group_by(seed, type, seq) %>% 
  summarise(avg_first = mean(first), avg_last = mean(last)) 
save(total, file = "C:/Users/p288427/Desktop/hd_rand_evo/test_first/total_test.R") 

####loading united dataframe####

for(j in types){
  if(j == "normal"){
    
    setwd(first_s)
    load("seqex_test_demog.R")
    
  } else if (j == "death"){
    
    setwd(first_death)
    load("death_test_demog.R")
    
  }else if (j == "eden"){
    
    setwd(first_eden)
    load("eden_test_demog.R")
    
  } else if (j == "sel"){
    
    setwd(first_sel)
    load("sel_test_demog.R")
    
  }
}

load("C:/Users/p288427/Desktop/hd_rand_evo/test_first/total_test.R")
####Create color plaettes####
n_colors = 100
rbg<- colorRampPalette(c("red", "blue", "green"))(n_colors)
cub_hel <- rev(cubehelix(n_colors))

####Plot first and last####
for( i in types)
{
  if(i == "normal")
  {
    test_demog = seqex_test_demog
  } else if ( i == "eden")
  {
    test_demog = eden_test_demog
  } else if ( i == "death")
  {
    test_demog = death_test_demog
  }
  
  ####Check envs are the same####
  
  ggplot(data = test_demog) +
    geom_tile(aes(condition, seed, fill = as.factor(env_p_3) )) +
    facet_grid(gen  ~ .) +
    theme(legend.position = "none") 
  
  ####Plot dem####
  sp = ggplot(data =  test_demog) +
    geom_tile(aes(x = condition, y = seed, fill = spore))+
    scale_fill_gradientn("spore", colors = rbg)+
    ggtitle(paste(i,"_spore",sep = ""))+
    facet_grid(.  ~ gen)
  
  # print(sp)
  
  su = ggplot(data =  test_demog) +
    geom_tile(aes(x = condition, y = seed, fill = success))+
    scale_fill_gradientn("success", colors = rbg)+
    ggtitle(paste(i,"_success",sep = ""))+
    facet_grid(.  ~ gen)
  
  # print(su)
  
  scaled = ggplot(data =  test_demog %>%  
                    pivot_longer(c(spore, success)) %>%
                    subset(seed != "9") %>% 
                    group_by(name) %>% 
                    mutate(value = scale(value))) +
    geom_tile(aes(x = condition, y = seed, fill = value))+
    scale_fill_gradientn("value", colors = rbg)+
    ggtitle(paste("scaled_all_",i, sep = ""))+
    facet_nested(seq + name  ~ gen)
  
  print(scaled)
}


####Plotting Heatmap clusters####
for(i in types)
{ 
  if(i == "normal"){
    test_demog = seqex_test_demog
  } else if ( i == "eden"){
    test_demog = eden_test_demog
  } else if ( i == "death"){
    test_demog = death_test_demog
  }else if ( i == "sel"){
    test_demog = sel_test_demog
  }
  #create breaks
  qntl = quantile(test_demog$success)
  sections = unique(
    c( seq(qntl[1], qntl[2], length = as.integer(n_colors/3 + 1)),
       seq(qntl[2], qntl[4], length = as.integer(n_colors/3 + 1)),
       seq(qntl[4], qntl[5], length = as.integer(n_colors/3 + 2))))
  
  clust_test_demog_first = 
    test_demog %>% 
    subset(gen == "first") %>% 
    subset(seed != 9) %>% 
    subset(seq == 1) %>%
    select(seed, condition, success) %>% 
    spread(condition, success)  %>% 
    column_to_rownames(var = "seed") %>% 
    select(-seq)
  
  start_col_test = as.dendrogram(hclust(dist(t(as.matrix(clust_test_demog_first)))))
  start_row_test = as.dendrogram(hclust(dist((as.matrix(clust_test_demog_first)))))
  
  heatmap.2(as.matrix(clust_test_demog_first ),
            trace = "none",
            Colv = start_col_test,
            Rowv = start_row_test,
            col = rbg,
            # breaks = sections,
            main = paste("first_",i))
  
  clust_test_demog_last = 
    test_demog %>% 
    subset(gen == "last") %>% 
    subset(seed != 9) %>% 
    subset(seq == 1) %>%
    select(seed, condition, success) %>% 
    spread(condition, success) %>% 
    column_to_rownames(var = "seed") %>% 
    select(-seq)
  
  heatmap.2(as.matrix(clust_test_demog_last),
            trace = "none",
            Colv = start_col_test,
            Rowv = start_row_test,
            col = rbg,
            # breaks = sections,
            main = paste("last_",i))
}
####Boxplots#####
  total %>% 
  filter(type == "sel") %>% 
  pivot_longer(c(avg_first,avg_last)) %>% 
  # group_by(seed) %>% 
  # filter(any(gen == "first" & success > quantile(test_demog$success, 0.9))) %>%
  ggplot(aes(x=name, y=value, fill=name)) +
  geom_violin() +
  geom_boxplot(width = 0.5) +
  facet_grid( . ~ .)

####Correlation tests####
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
    
  }else if ( i == "sel")
  {
    
    setwd(first_sel)
    load("sel_test_demog.R")
    test_demog = sel_test_demog
    
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
  # Shapiro-Wilk normality test for "first"
  sf = shapiro.test(cor_df$first) 
  print(sf)
  # Shapiro-Wilk normality test for "last"
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
##### Overall Plot######


#scatter_plot for seed average
my.formula <- y ~ x
data = total %>% 
  # filter(type == "sel") %>%
  filter(type != "eden") %>% 
  rename(first = avg_first, last = avg_last)

data = sel_test_demog %>% rbind(death_test_demog) %>% rbind(seqex_test_demog) %>% 
  pivot_wider(id_cols = c(condition, seed, seq, type), names_from = gen, values_from = success)

ggscatter(data = data,
          x = "first",
          y = "last",
          color = "type",
          size = 0.01,
          fullrange=T,
          #label = "seed",
          add = "reg.line", 
          conf.int = F, 
          cor.coef = T,
          cor.method = "spearman") +
  stat_poly_eq( data = data,
                formula = my.formula, 
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE,
                label.y = "bottom",
                label.x = "right") +    
  #add line for 1 to 1 cor
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position = "none") +
  facet_grid(type ~ .)

###Plotting percentages of tests that are above the x=y line####
percentage = total %>% 
  group_by(type, seq) %>% 
  summarise(above_line = sum(avg_first < avg_last)/length(avg_first)) %>% 
  mutate(below_line = 1 - above_line) %>% 
  pivot_longer(c(above_line, below_line)) %>%
  group_by(seq, type) %>% 
  mutate(ypos = cumsum(value) - 0.5 * value )

###Plotiing corr between start level and how much above/below the equality line you are####

#scatter_plot for seed average
my.formula <- x ~ y
ggscatter(data = total %>%
            mutate(improvement_after_2nd_round = avg_last - avg_first),
          x = "avg_first",
          y = "improvement_after_2nd_round",
          color = "type",
          fullrange=TRUE,
          #label = "seed",
          add = "reg.line", 
          conf.int = F, 
          cor.coef = TRUE,
          cor.method = "spearman",
          main = "first last delta equality correlations")  +
  #add line for 1 to 1 cor
  facet_grid(seq ~ type)

####Plotting correlation btw avg_last and avg_success over evo####
#[need data from tot_dem in plot_rand_evo_demographic.R]
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

my.formula <- y ~ x
ggscatter(data = tot_first_last_evo_suc %>% filter(type != "seqex"),
          x = "avg_success",
          y = "avg_first",
          size = 0.01,
          fullrange=T,
          #label = "seed",
          add = "reg.line", 
          conf.int = T, 
          cor.coef = T,
          cor.method = "spearman") +
  stat_poly_eq( data = tot_first_last_evo_suc,
                formula = my.formula, 
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE,
                label.y = "bottom",
                label.x = "right") +    
  #add line for 1 to 1 cor
  #geom_abline(intercept = 0, slope = 1) +
  theme(legend.position = "none") +
  facet_grid(type ~ .)

ggplot(data = tot_first_last_evo_suc %>% 
         pivot_longer(c(avg_first, avg_last, avg_pre_adap, avg_improvement, avg_success))) +
  geom_violin(aes(x = name, y = value, fill = name))+
  facet_grid(type ~ seed)
