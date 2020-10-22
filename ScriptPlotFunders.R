library(tidyr)
library(tidyverse)
library(ggplot2)
library(stringi)
library(rlist)
library(igraph)
library(ggraph)
library(networkD3)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(dir,"/vG_W_2016_data",sep = ""))

########Plot Philogenesis###############
# funders_success = data.frame()
# for (i in list.files(path = '.',pattern = "funders_success_s1.csv"))
# {
#   print(i)
#   sim_run = read.csv(i, header = FALSE, sep = ',')
#   funders_success = rbind(sim_run,funders_success)
# }
# 
# 
# #simple plot to see if there are any super fit individuals
# simple_frame = funders_success[,c(1,length(funders_success))]
# ggplot(simple_frame, aes(x = V1, y = V54))+
#   geom_point()+
#   geom_smooth()
# 
# #we take only last generation
# last_gen = as_tibble(funders_success[funders_success$V1 == max(funders_success$V1) - 1,2])
# last_gen$value = as.character(last_gen$value) 
# 
# for (i in 1:length(last_gen$value)) {
#   last_gen$value[i] = 
#     substring(last_gen$value[i],
#               first = 2,
#               last = nchar(last_gen$value[i]) - 2)
# }
# 
# get_parts <- function(x){
#   
#   parts <- c(unlist(strsplit(x, split = " ")))
#   
#   parts
# }
# 
# test = get_parts(last_gen$value[1])
# 
# for (i in 2 : (nrow(last_gen))) {
#   test = qpcR:::rbind.na(
#     test,
#     get_parts(last_gen$value[i])
#     )
# }
# test = as.data.frame(test)
# colnames(test) = paste("level",c(1:ncol(test)), sep = "")
# 
# edge_list <- test %>% select("level1","level2") %>% unique %>% rename(from="level1", to="level2")
# edge_list$from = paste("level", toString(1), "_", edge_list$from, sep="")
# edge_list$to = paste("level", toString(2), "_", edge_list$to, sep="")
# for(j in 2:(ncol(test) - 1))
# {
#   edges <- 
#     test %>% select(colnames(test)[j], colnames(test)[j + 1]) %>% unique  %>% rename(from=colnames(test)[j], to=colnames(test)[j + 1])
#   edges$from = paste("level", toString(j), "_",edges$from, sep="")
#   edges$to = paste("level", toString(j + 1), "_", edges$to, sep="")
#   edge_list = rbind(edge_list, edges)
# }
# 
# 
# saveRDS(edge_list,"edge_list_S1_last_gen")
loaded_edge_list = readRDS("edge_list_S1_last_gen")

mygraph <- graph_from_data_frame( loaded_edge_list )

#saveRDS(mygraph,file = "phylogenesis_s1")
phylo = readRDS("phylogenesis_s1")

#Testing different possible layouts
ggraph(mygraph, layout = 'tree') + 
  geom_edge_diagonal() +
  geom_node_point() +
  theme_graph()

ggraph(mygraph, layout = 'tree') + 
  geom_edge_elbow() + 
  coord_fixed()

ggraph(mygraph, layout = 'circlepack') + 
  geom_node_circle(aes(fill = depth), size = 0.25, n = 50) + 
  coord_fixed()

#Change layout so all nodes in same generation have x from 1:100
mylayout = create_layout(mygraph, 'tree')
mylayout$x = as.numeric(stringi::stri_match(layout$name, regex = "\\d+$" ))
 
ggraph(mylayout)+
  geom_edge_diagonal() +
  geom_node_point() +
  theme_graph()

#########Plot Networks###############
as_tibble(funders_success)

funder_network = funders_success[1,]
funder_network = funder_network[,-c(1,2)]
funder_network = funder_network[,1:(length(funder_network)- 6)]

as_tibble(funder_network)


layers = list()
node = list()
iterator = 1
last_node_iterator = 1
layer_threshold_iterator = 2 #important, because it skips input layer
node_threshold_iterator = 1
node_state_iterator = 1
for(value in funder_network)
{
  
  if(is.character(value) && length(layers) < 3)
  {
    if(value == "  ! " )
    {
      weights = as.vector(funder_network[,last_node_iterator : (iterator - 1)])
      names(weights) = c(paste0("w",1:length(weights)))
      node = list.append(node, weights)
      names(node)[length(node)] = paste("Node_", length(node), sep = "")
      print("node!")
      last_node_iterator = iterator + 1
    }
    else if(value == "  | ")
    { 
      print("layer!")
      layers = list.append(layers, node)
      names(layers)[length(layers)] = paste("Layer_", length(layers), sep = "")
      node = list()
      last_node_iterator = iterator + 1
    }
  }
  else if(is.character(value) && 
          length(layers) >= 3 && 
          layer_threshold_iterator <= 3)
  {
    if(value == "  | ")
    {
print("node_treshold")

      for( i in funder_network[,last_node_iterator : (iterator - 1)])
      {
        print(i)
        
        layers[[layer_threshold_iterator]][[node_threshold_iterator]]$threshold = i
        node_threshold_iterator = node_threshold_iterator + 1 
      }
      last_node_iterator = iterator + 1
      layer_threshold_iterator = layer_threshold_iterator + 1
      node_threshold_iterator = 1
    }
  }
  else if(layer_threshold_iterator > 3)
  {
    layers[[2]][[node_state_iterator]]$state = value
    node_state_iterator = node_state_iterator + 1
    
  }
  iterator = iterator + 1
}

names(layers)= c("Layer1", "Layer2", "Layer2")
names(layers)[length(layers)] = paste("Layer_", length(layers), sep = "")
layers_df = as.data.frame(layers)
