library(tidyr)
library(dplyr)
library(stringr)
library(rgl)
library(misc3d)
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(dir,"/vG_W_2016_data",sep = ""))

all_reac_norms = data.frame()

for (i in  list.files(path = '.',pattern = "reaction_norm_best_ind_s\\d+_f\\d+.csv"))
{
  reac_norm = read.table(i, sep = ",")
 #Make file smaller otherwise size is too big
 reac_norm = reac_norm[ as.numeric(row.names(reac_norm)) %% 2 == 0,]
 reac_norm$seed = sub( "^.*s(\\d+).*",'\\1', i)
 reac_norm$change = sub( "^.*f(\\d+).*",'\\1', i)
 all_reac_norms = rbind(reac_norm,all_reac_norms)

}
all_reac_norms$seed = as.factor(all_reac_norms$seed)
all_reac_norms$change = as.factor(all_reac_norms$change)
green <- rgb( 0, 255, 0, max = 255, alpha = 0.5, names = "t")
blue <- rgb(0, 0, 255, max = 255, alpha = 0, names = "o")

all_reac_norms$col = as.factor(ifelse(all_reac_norms$V4 == 0, blue,transparent))
all_reac_norms = all_reac_norms[ as.numeric(row.names(all_reac_norms)) %% 2 == 0,]
colnames(all_reac_norms)[1:4] = c("Energy", "Metabolite", "Food", "Sporulation")
save(all_reac_norms,file = "all_reac_norms.R")##

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

for( j in levels(all_reac_norms$change))
{
#subset for a certain frequency of change
freq = subset(all_reac_norms, change == j)[,c("Sporulation","seed")]
#resize factors in seed column
freq$seed = factor(freq$seed)

#create empty data frame
reac_norms_freq = data.frame()

#bind sporulation columns for each seed as rows in reac_norms_freq
for (i in levels(freq$seed)) {
  reac = freq %>% 
  subset(seed == i)
  reac_norms_freq = rbind(reac[,1] , reac_norms_freq)
}

#rename rows so that it corresponds to seed
row.names(reac_norms_freq) = levels(freq$seed)

#calculate distance matrix between rows(sporulation responses)
#using euclidean method
#and also normalizing it so values are betwen 0 and 1
d = dist(reac_norms_freq, method ="euclidean", diag = T, upper = T)
dd = d/sqrt(length(reac_norms_freq))

#create title for graph
title = paste("freq",j, sep = "_")

#plot dendrograms using hierarchical clustering
pdf(paste( "hd",title, sep = "_"),
    width = 8.27,
    height = 11.69)
plot(hclust(d), main = title )
dev.off()

pdf(paste( "hd_averaged",title, sep = "_"),
    width = 8.27,
    height = 11.69)
plot(hclust(dd), main = title)
dev.off()

#using neighbour joining alghoritm(same as jordi?)
library(ape)
pdf(paste( "ud",title, sep = "_"),
    width = 8.27,
    height = 11.69)
plot(nj(d),"u", main = title)
dev.off()

pdf(paste( "ud_averaged",title, sep = "_"),
    width = 8.27,
    height = 11.69)
plot(nj(dd),"u", main = title)
dev.off()

#For other options see: 
#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
#using ggplot?
#using extendeddendr
}
