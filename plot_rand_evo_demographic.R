library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(pals)
library(gplots)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")
hd_rand_evo = "C:/Users/p288427/Desktop/hd_rand_evo"


#####read data####
setwd(hd_rand_evo)
demographic = data.frame()


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

n_columns = ncol(demographic)
colnames(demographic)= c("cycle",
                         "active",
                         "spore",
                         "sporu" ,
                         "n_timesteps",
                         sprintf("env_p_%s",seq(1:(n_columns - 5 - 3))),
                         "seed",
                         "change_freq",
                         "condition")



demographic$seed = as.factor(demographic$seed)
demographic$change_freq = as.factor(demographic$change_freq)
demographic$condition = as.factor(demographic$condition)
demographic$n_timesteps = as.numeric(demographic$n_timesteps)

# create new columns for ratio of spore produced and starting production
demographic = demographic %>%
  group_by(seed,condition,cycle) %>%
  mutate(total_n = sum(c_across(c(spore, sporu, active)))) %>% 
  ungroup() %>%
  group_by(condition,cycle) %>%
  mutate(
    "ratio_value" = spore / max(spore)
  ) %>%
  ungroup() %>%
  group_by(condition, seed) %>%
  mutate("ratio_start_production" = ratio_value[min(cycle)]) %>%
  mutate("ratio_end_production" = ratio_value[max(cycle)]) %>%
  mutate("delta_rv_start_end" = ratio_end_production - ratio_start_production) %>%
  mutate("start_production" = spore[min(cycle)])  %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate("standardized_delta_rv_start_end" = delta_rv_start_end / max(delta_rv_start_end)) %>%
  mutate("overall_r_value" = spore / max(spore)) %>% 
ungroup()

hd_demographic  = demographic
save(hd_demographic, file = "hd_rand_evo_demo.R")

####load hd_rand condition from 0 : 49####
#object name : hd_demographic
setwd(hd_rand_evo)
load(file = "hd_rand_evo_demo.R")

dem = hd_demographic %>% subset(change_freq == "0")

####make sure each cycle is either 125 timesteps or around 10000 inds####
c_check = dem %>% 
  select(cycle,n_timesteps,total_n,seed,condition)

#should give back empty tibble
c_check %>% 
  group_by(seed, condition,cycle) %>% 
  subset(n_timesteps != 125) %>% 
  subset(total_n < 9000)

####check whihc environmental parameters induce####

#summarise conditions and timestep
cond_tmstps = dem %>% 
  group_by(condition) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(starts_with("env_")|starts_with("n_time"),total_n) 

ggplot(data = cond_tmstps)+
  geom_tile(aes(x = condition, y = seed, fill = n_timesteps))

#select only relevant conditions
plot_tmstps = cond_tmstps %>% 
  select(env_p_4,env_p_3,env_p_10,env_p_15,env_p_19, env_p_22,n_timesteps, total_n) %>%
  pivot_longer(env_p_4 : env_p_19, names_to = "variable")

ggplot(data = plot_tmstps)+
  geom_point(aes(x = value , y = total_n))+
  geom_smooth(aes(x = value, y = total_n))+
  facet_grid(.  ~ variable )




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

dem_d = dem %>% 
  left_join(dem_c)

dem_d$env_type = as.factor(dem_d$env_type)

ggplot(data = dem_d %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = env_type))

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

n_colors = 20
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

####for ratio_value
rv_qntl = quantile(dem$ratio_value)
sections_rv_value = unique(
  c( seq(rv_qntl[1], rv_qntl[2], length = as.integer(n_colors/3 + 4)),
     seq(rv_qntl[2], rv_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(rv_qntl[4], rv_qntl[5], length = as.integer(n_colors/3 ))))

####for value
v_qntl = quantile(dem$spore)
sections_v_value = unique(
  c( seq(v_qntl[1], v_qntl[2], length = as.integer(n_colors/3 + 4)),
     seq(v_qntl[2], v_qntl[4], length = as.integer(n_colors/3 + 1)),
     seq(v_qntl[4], v_qntl[5], length = as.integer(n_colors/3))))


###Plotting ratio-value beginning and end plus points of value####

ggplot(data = dem_all_tmstps) +
  # geom_rect(aes(ymin=min(spore),
  #               ymax= max(spore)  + 1,
  #               xmin=min(cycle),
  #               xmax= max(cycle) / 2,
  #               fill = ratio_start_production), alpha =0.5) + 
  # geom_rect(aes(ymin=min(spore),
  #               ymax= max(spore)  + 1,
  #               xmin= max(cycle) / 2,
  #               xmax= max(cycle),
  #               fill = ratio_end_production), alpha =0.5) +
  # geom_point(aes(cycle,spore)) +
  # scale_fill_gradientn("Ratio",colors = cubehelix(10)) +
  geom_line(aes(cycle,spore)) +
  facet_grid(seed ~ condition)

ggsave("../research presentation/rv_s_e_p_only_complete_cycle.png",
       width = 50,
       height = 30, 
       units = "cm")

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
          breaks = sections_v_value)



#plot end with clust_start
heatmap.2(as.matrix(clust_dem_v_end) ,
          Rowv = start_clust_row,
          Colv = start_clust_col,
          scale = "none",
          trace = "none",
          col = cub_hel,
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

####### Fitting logistic regression model to data and adding it to dataframe####

fit_logistic <-
  function(value, cycle) {
    out <- tryCatch(
      {
        nls(value ~ SSlogis(cycle, Asym, xmid, scal))
      },
      error=function(cond){
        second_fit <- tryCatch(
          {
            #Estimation of starting parameter, and formula for model taken from here: 
            # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
            # This works for datapoints that describe the ascending part of the logisitc function already past the flexion point
            c.0 = min(value) * 0.5
            model.0 = lm(log(value - c.0) ~ cycle)
            start = list(a=exp(coef(model.0)[1]), b=coef(model.0)[2], c=c.0)
            nlc = nls.control(maxiter = 1000, minFactor = 1/10^20)
            nls(value ~ a * exp(b * cycle) + c, start = start, control = nlc)
          },    
          error=function(cond) 
          {
            return(NA)
          }
          
        )
        return(second_fit)
      }
    )    
    return(out)
  }

me = 
  demographic %>% 
  group_by(condition, seed) %>% 
  mutate( "fit" = list(fit_logistic(value, cycle))) %>% 
  ungroup()%>% 
  drop_na()

mef = me %>%   subset(cycle == 1) 

me_nested =
  me %>% 
  group_by(condition, seed) %>% 
  nest()

mef = me %>% 
  select(c(fit)) %>% 
  mutate( "fit_prediction" = list(predict(fit[[1]],seq(1,499))))


###plotting three params of regressions####
###NOT WORKING!!!!!!
ggplot(data = me) +
  geom_line(aes(x = cycle, y = predict(fit, seq(1,499), data = me %>% group_by(seed, condition)))) + 
  facet_grid(seed ~ condition)


###working on a single plot######

o = dem_all_tmstps %>% 
  subset(condition ==  "9" ) %>% 
  subset(seed == "90")

coef(o$fit[[1]])[2]
#Estimation of starting parameter, and formula for model taken from here: 
# https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
# This works for datapoints that describe the ascending part of the logisitc function already past the flexion point
c.0 = min(of$value) * 0.5
model.0 = lm(log(value - c.0) ~ cycle, data = of)
start = list(a=exp(coef(model.0)[1]), b=coef(model.0)[2], c=c.0)
nlc = nls.control(maxiter = 1000, minFactor = 1/10^20)
model = nls(value ~ a * exp(b * cycle) + c, data = of, start = start, control = nlc)
coeffit = as.data.frame(coef(model))

# this work for other datapoints(full sigmoid or first half)
model1 <- nls(value ~ SSlogis(cycle, Asym, xmid, scal), data = of)


ggplot(data = pivot_longer(o,c(spore, active, sporu, total_n), names_to = "variable"))+
  geom_line(aes(x = cycle, y = value, color = variable))


lines(
  seq(1,499),
  predict(o$fit[[1]], seq(1,499))
)

#plot asymptote
lines(of$cycle, 
      rep(exp(coef(model)[1]),length(of$cycle))
)

#plot midpoint
abline(v = exp(coef(model)[2]))
