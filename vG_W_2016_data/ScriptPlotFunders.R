library(tidyr)
library(ggplot2)
library(stringr)
library(rlist)
funders_success = data.frame()
list.files(path = '.',pattern = "funders_success_s*")

for (i in list.files(path = '.',pattern = "funders_success_s*"))
{
  print(i)
  sim_run = read.csv(i)
  funders_success = rbind(sim_run,funders_success)
}

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
