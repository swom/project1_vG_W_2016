library(tidyr)
library(dplyr)
library(stringr)
library(rgl)
library(misc3d)
all_reac_norms = data.frame()

for (i in  list.files(path = '.',pattern = "reaction_norm_best_ind_s1_f0.csv"))
{
  reac_norm = read.table(i, sep = ",")
  reac_norm$seed = sub( "^.*s(\\d+).*",'\\1', i)
  reac_norm$change = sub( "^.*f(\\d+).*",'\\1', i)
  all_reac_norms = rbind(reac_norm,all_reac_norms)
  
}

#save(all_reac_norms,file = "all_reac_norms.R")
load("all_reac_norms.R")


#Making colors for graph
transparent <- rgb( 0, 255, 0, max = 255, alpha = 0.5, names = "blue50")
blue <- rgb(0, 0, 255, max = 255, alpha = 0, names = "o")

reac_norm$col = ifelse(reac_norm$V4 == 0, blue,transparent)
small_reac_norm = reac_norm[ as.numeric(row.names(reac_norm)) %% 2 == 0,]

points3d( 
  x= small_reac_norm$V1, y= small_reac_norm$V2, z= small_reac_norm$V3,
  col = small_reac_norm$col
)
box3d(col = "gray")

library(akima)
require(rgl)
xyz = small_reac_norm %>% subset(V4 == 1)
xyz = xyz[,1:3] 
mesh <- as.mesh3d(xyz, triangles = FALSE, col = "red")
mesh$vb
mesh$ib
open3d()
shade3d(mesh)



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

fig <- plot_ly(small_reac_norm, 
               x = ~V1, y = ~V2, z = ~V3, color = ~V4, colors = c(transparent,blue), alpha = 0.1)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Energy'),
                                   yaxis = list(title = 'Food'),
                                   zaxis = list(title = 'Metabolite')))
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

fig
