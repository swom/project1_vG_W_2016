library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(pals)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")
hd_rand_evo = "C:/Users/p288427/Desktop/hd_rand_evo"
hd_rand_evo_1 = "C:/Users/p288427/Desktop/hd_rand_evo_1"

setwd(rand_evo_dir)
demographic = data.frame()


for (i in  list.files(path = '.',pattern = "rand_evo_a3.000000cond_\\d+sim_demographic_s\\d+change_\\d+"))
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
                          sprintf("env_p_%s",seq(1:n_columns - 5 - 3)),
                         "seed",
                         "change_freq",
                         "condition")


demographic <- pivot_longer(
  demographic,
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)

demographic$seed = as.factor(demographic$seed)
demographic$change_freq = as.factor(demographic$change_freq)
demographic$variable = as.factor(demographic$variable)
demographic$condition = as.factor(demographic$condition)
demographic$n_timesteps = as.factor(demographic$n_timesteps)

# create new columns for ratio of spore produced and starting production
demographic = demographic %>%
  subset(variable == "spore") %>%
  group_by(condition,cycle) %>%
  mutate(
    "ratio_value" = value / max(value)
  ) %>%
  ungroup() %>%
  group_by(condition, seed) %>%
  mutate("ratio_start_production" = ratio_value[min(cycle)]) %>%
  mutate("ratio_end_production" = ratio_value[max(cycle)]) %>%
  mutate("delta_rv_start_end" = ratio_end_production - ratio_start_production) %>%
  mutate("start_production" = value[min(cycle)])  %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate("standardized_delta_rv_start_end" = delta_rv_start_end / max(delta_rv_start_end)) %>%
  ungroup() %>%
  subset(condition != 50)

# save(demographic, file = "rand_evo_demo.R")

setwd(rand_evo_dir)
load(file = "rand_evo_demo.R")

####load hd_rand condition from 0 : 49###
setwd(hd_rand_evo)
load(file = "hd_rand_evo_demo.R")

hd_demographic = hd_demographic %>% subset(condition != 0)
hd_demographic$seed = as.factor(as.numeric(hd_demographic$seed) + 900) 

####load hd_rand condition from 1 : 50###
setwd(hd_rand_evo_1)
load(file = "hd_rand_evo_demo_1.R")

hd_demographic_1$seed = as.factor(as.numeric(hd_demographic_1$seed) + 900) 
hd_demographic_1 %>% mutate("delta_rv_start_end" = ratio_end_production - ratio_start_production)
                            


demographic %>%
  mutate("delta_rv_start_end" = ratio_end_production - ratio_start_production)

demo = rbind(demographic, hd_demographic_1)

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

###Plotting ratio-value beginning and end plus points of value####

ggplot(data = demographic  %>% subset(change_freq == 0)) +
  geom_rect(aes(ymin=min(value),
                ymax= max(value)  + 1,
                xmin=min(cycle),
                xmax= max(cycle) / 2,
                fill = ratio_start_production), alpha =0.5) + 
  geom_rect(aes(ymin=min(value),
                ymax= max(value)  + 1,
                xmin= max(cycle) / 2,
                xmax= max(cycle),
                fill = ratio_end_production), alpha =0.5) +
  # geom_point(aes(cycle,value)) +
  scale_fill_gradientn("Ratio",colors = cubehelix(10)) +
  geom_line(aes(cycle,value)) +
  facet_grid(seed ~ condition)


####Plotting diff between rvalue end vs start normalized####
ggplot(data = demographic %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = delta_rv_start_end),
            color = "black", size = 0.5) + 
  scale_fill_gradientn("Delta", colors = rev(cubehelix(10))) 

###Plotting start or end value####
###ratio_value

ggplot(data = demo %>% subset(cycle == min(cycle))) +
  geom_tile(aes(condition, seed, fill = ratio_value),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))

ggplot(data = demographic %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = ratio_value),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))



###value
ggplot(data = demographic %>% subset(cycle == min(cycle))) +
  geom_tile(aes(condition, seed, fill = value),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))

ggplot(data = demographic %>% subset(cycle == max(cycle))) +
  geom_tile(aes(condition, seed, fill = value),
            color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))

####heatmap cluster for production of spores####
col<- colorRampPalette(c("red", "blue", "green"))(20)

#start
clust_dem_start = demographic %>% 
  subset(cycle == min(cycle)) %>% 
  select(seed, ratio_value, condition) %>% 
  spread(condition, ratio_value) %>% 
  column_to_rownames(var = "seed")


heatmap(as.matrix(clust_dem_start) , scale = "none",
        col = col)


#end
clust_dem_end = demographic %>% 
  subset(cycle == max(cycle)) %>% 
  select(seed, ratio_value, condition) %>% 
  spread(condition, ratio_value) %>% 
  column_to_rownames(var = "seed")


heatmap(as.matrix(clust_dem_end) , scale = "none",
        col = col)


####heatmap cluster for delta in ranking####
col<- colorRampPalette(c("red", "blue", "green"))(20)

#start
clust_delta = demographic %>% 
  subset(cycle == max(cycle)) %>% 
  select(seed, delta_rv_start_end, condition) %>% 
  spread(condition, delta_rv_start_end) %>% 
  column_to_rownames(var = "seed")


heatmap(as.matrix(clust_delta) , scale = "none",
        col = rev(cubehelix(10)))
####Plotting variance of production over time for conditions####
var_demo = demographic %>% 
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

o = me %>% 
  subset(condition ==  "42" ) %>% 
  subset(seed == "1")

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

#plotting data points
plot(rep(10000, 499) ~ seq(1,499), ylim = c(0,6000) )
#lines
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
