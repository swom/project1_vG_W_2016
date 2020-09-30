library(tidyr)
library(dplyr)
library(stringr)
library(rgl)
library(misc3d)
all_reac_norms = data.frame()

#for (i in  list.files(path = '.',pattern = "reaction_norm_best_ind_s\\d+_f\\d+.csv"))
#{
#  reac_norm = read.table(i, sep = ",")
#  #Make file smaller otherwise size is too big
#  reac_norm = reac_norm[ as.numeric(row.names(reac_norm)) %% 2 == 0,]
#  reac_norm$seed = sub( "^.*s(\\d+).*",'\\1', i)
#  reac_norm$change = sub( "^.*f(\\d+).*",'\\1', i)
#  all_reac_norms = rbind(reac_norm,all_reac_norms)
#  
###}
#all_reac_norms$seed = as.factor(all_reac_norms$seed##)
#all_reac_norms$change = as.factor(all_reac_norms$change)##
#green <- rgb( 0, 255, 0, max = 255, alpha = 0.5, names = "t")
#blue <- rgb(0, 0, 255, max = 255, alpha = 0, names = "o")

#all_reac_norms$col = as.factor(ifelse(all_reac_norms$V4 == 0, blue,transparent))
#all_reac_norms = all_reac_norms[ as.numeric(row.names(all_reac_norms)) %% 2 == 0,]
#save(all_reac_norms,file = "all_reac_norms.R")##

load("all_reac_norms.R")

###Plotting the reaction Norms

#Making colors for graph



library(plotly)

axx <- list(
  nticks = 4,
  range = c(0,20)
)

axy <- list(
  nticks = 4,
  range = c(0,20)
)

axz <- list(
  nticks = 4,
  range = c(0,20)
)

for (i in levels(all_reac_norms$seed)) {
  for (j in levels(all_reac_norms$change)) {
    fig <-plot_ly(all_reac_norms %>% subset(seed == i & change == j ), 
                   x = ~V1, y = ~V2, z = ~V3, color = ~V4, colors = c(blue,green), alpha = 0.1
                   )
    fig <- fig %>% layout(title = paste("RN_",i,"_",j, sep = ""))
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'Energy'),
                                       yaxis = list(title = 'Food'),
                                       zaxis = list(title = 'Metabolite')))
    fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
    fig
    htmlwidgets::saveWidget(fig, paste("RN_",i,"_",j,".html", sep = ""))
  }
}


###Calculating distances and making phenograms

