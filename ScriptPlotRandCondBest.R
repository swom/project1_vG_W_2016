library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pals)
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(paste(dir,"vG_W_2016_data/evo", sep = "/"))
setwd("C:/Users/p288427/Desktop/500cycles old sim_par")
rand_demographic = data.frame()

for (i in  list.files(path = '.',pattern = "best_ind_random_cond_sim_demographic_s\\d+_change_\\d+_amplitude_\\d+.*"))
{
  conditions = read.table(i, sep = ",")
  conditions$seed = sub( "^.*s(\\d+).*",'\\1', i)
  conditions$change = sub( "^.*change_(\\d+).*",'\\1', i)
  conditions$amplitude = sub( "^.*amplitude_(\\d+)",'\\1', i)
  conditions$amplitude = sub( "(\\d+).csv",'\\1', conditions$amplitude)
  rand_demographic = rbind(conditions,rand_demographic)
}


n_columns = ncol(rand_demographic)
#store env param in a different data frame
env_params = rand_demographic[,c(5: (n_columns - 3) )]
#Remove the params, and name relevant variables
rand_demographic_slim = rand_demographic[,-c(5: (n_columns - 3) )]
colnames(rand_demographic_slim)= c("condition",
                                   "active",
                                   "spore",
                                   "sporu", 
                                   "seed",
                                   "change_freq",
                                   "amplitude")
#reattach the parameters
rand_demographic = cbind2(rand_demographic_slim,env_params)
rand_demographic <- pivot_longer(
  rand_demographic, 
  cols = c("active", "spore", "sporu"),
  names_to = "variable"
)

rand_demographic$amplitude = as.factor(rand_demographic$amplitude)
rand_demographic$change_freq = as.factor(rand_demographic$change_freq)
rand_demographic$seed = as.factor(rand_demographic$seed) 
rand_demographic$condition = as.factor(rand_demographic$condition) 
rand_demographic$variable = as.factor(rand_demographic$variable)


#Plot spore production on heatmap

f = 
  rand_demographic %>%
  subset(variable == "spore") %>% 
  group_by(amplitude, condition) %>% 
  mutate(
    "ratio_value" = value / max(value)
  )

f %>% 
  subset(variable == "spore") %>% 
  ggplot(aes(condition, seed, fill = ratio_value)) + 
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradientn(colors = rev(cubehelix(10)))+
  facet_grid(change_freq  ~ amplitude)

#Now only look at chnage_freq 0 and amplitude 3
# and calculate the standard deviation
# for relative spor production
# by environment and by population


f0a3 = f %>% 
  subset(change_freq == "0") %>% 
  subset(amplitude == "1.500000") %>% 
  group_by(condition) %>% 
  mutate( sd_cond = sd(ratio_value)) %>% 
  mutate(avg_cond = mean(ratio_value)) %>% 
  ungroup() %>% 
  group_by(seed) %>% 
  mutate(sd_pop = sd(ratio_value)) %>%
  mutate(avg_pop = mean(ratio_value)) %>% 
  ungroup()

#Find 6 environments
#3 with relatively constant production: high, medium, low
#3 with variable production: high, medium, low
# AND
#Find 12 populations
#3 with relatively constant production: high, medium, low
#3 with variable production: high, medium, low

#find in which of the 3 quantile(0.33,0.66,1) 
# is each condition's stand dev in relative spore production
f0a3 = f0a3 %>% 
  group_by(condition) %>% 
  mutate(sd_cond_quantile = 
           if(sd_cond < quantile(f0a3$sd_cond,0.33))
           {1}
         else if(sd_cond < quantile(f0a3$sd_cond,0.66))
         {2}
         else
         {3}) %>% 
  ungroup() %>% 
#find in which of the 3 quantile(0.33,0.66,1) 
# is each population's stand dev in relative spore production  
  group_by(seed) %>% 
  mutate(sd_pop_quantile = 
           if(sd_pop < quantile(f0a3$sd_pop,0.33))
           {1}
         else if(sd_pop < quantile(f0a3$sd_pop,0.66))
         {2}
         else
         {3})

#find in which of the 3 quantile(0.33,0.66,1) 
# is each condition's average relative spore production
f0a3 = f0a3 %>% 
  group_by(condition) %>% 
  mutate(rel_prod_quantile_cond = 
           if(avg_cond < quantile(f0a3$avg_cond,0.33))
           {1}
         else if(avg_cond < quantile(f0a3$avg_cond, 0.66))
         {2}
         else
         {3}) %>% 
  ungroup() %>% 
  #find in which of the 3 quantile(0.33,0.66,1) 
  # is each population's average relative spore production  
  group_by(seed) %>% 
  mutate(rel_prod_quantile_pop = 
           if(avg_pop < quantile(f0a3$avg_pop,0.33))
           {1}
         else if(avg_pop < quantile(f0a3$avg_pop,0.66))
         {2}
         else
         {3}) %>% 
  ungroup()

#####create Jordi's color palette#####
drsimonj_colors <- c(
  `red`        = "#d11141",
  `green`      = "#00b159",
  `blue`       = "#00aedb"
)
#' Function to extract drsimonj colors as hex codes
#'
#' @param ... Character names of drsimonj_colors 
#'
drsimonj_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (drsimonj_colors)
  
  drsimonj_colors[cols]
}
#combine colors into palette
drsimonj_palettes <- list(
  `main`  = drsimonj_cols("red", "blue", "green"),
  `cool`  = drsimonj_cols("red", "green")
  
)
#' Return function to interpolate a drsimonj color palette
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette()
#'
drsimonj_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- drsimonj_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}
#' Color scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_drsimonj <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- drsimonj_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("drsimonj_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_drsimonj <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- drsimonj_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("drsimonj_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}

#####plot##### 
a =f0a3 %>%  
  # subset(as.numeric(seed) < 51) %>% 
  ggplot(aes(condition, seed, fill = ratio_value)) + 
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradientn(colours = rev(cubehelix(10))) +
facet_grid(change_freq  ~ amplitude)

#plot the different sd quantiles for pops
b =f0a3 %>%   
  subset(as.numeric(seed) < 51) %>% 
  ggplot(aes(condition, seed, fill = sd_pop_quantile)) + 
  geom_tile(color = "black", size = 0.5) 

#plot the different sd quantiles for conditions
c =f0a3 %>%   
  subset(as.numeric(seed) < 51) %>% 
  ggplot(aes(condition, seed, fill = sd_cond_quantile)) + 
  geom_tile(color = "black", size = 0.5) 

#plot the different avg_ratio_Value quantiles for pops
d = f0a3 %>%  
  subset(as.numeric(seed) < 51) %>% 
  ggplot(aes(condition, seed, fill = rel_prod_quantile_pop)) + 
  geom_tile(color = "black", size = 0.5) 

#plot the different avg_ratio_Value quantiles for conditions
e = f0a3 %>% 
  subset(as.numeric(seed) < 51) %>% 
  ggplot(aes(condition, seed, fill = rel_prod_quantile_cond)) + 
  geom_tile(color = "black", size = 0.5) 

ggarrange(b,c,d,e,
          labels = c("sd_q_pop", "sd_q_cond",
                     "avg_q_pop", "avg_q_cond"),
          ncol = 2, nrow = 2,
          align = "v")
