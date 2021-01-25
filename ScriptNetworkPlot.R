library(tidyr)
library(tidyverse)
library(ggplot2)
library(stringi)
library(rlist)
library(igraph)
library(ggraph)
library(networkD3)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
rand_evo_dir = paste(dir,"/vG_W_2016_data/rand_evo",sep = "")
evo_dir = paste(dir,"/vG_W_2016_data/evo",sep = "")

setwd(evo_dir)

funders = data.frame()

# for (i in  list.files(path = '.',pattern = "funders_success_s\\d+change_\\d+"))
# {
#   replicate = read.csv(i)
#   replicate$seed = sub( "^.*s(\\d+).*",'\\1', i)
#   replicate$change = sub( "^.*_(\\d+).*",'\\1', i)
#   colnames(replicate) = colnames(funders)
#   funders = rbind(replicate, funders)
# }

funders = read.csv("funders_success_s9change_0.csv")
funder= funders[1,]
funder_mod = funder[,3:(length(funders) - 6)]
funder_mod$seed = 9
funder_mod$change = 0

#assuming we know network architecture:
# 3-3-2, fully connected, with hidden layer fully connected to itself
n_nodes_l1 = 3
n_nodes_l2 = 3
n_nodes_l3 = 2

#get clean dataframe(no characters)
no_ch_funder = 
  funder_mod %>% select_if(is.numeric)

connections = data.frame()

#I2H connections
for( i in 1:(n_nodes_l1 * n_nodes_l2))
{
  ID_sender = 1 + ((i - 1) %/% n_nodes_l2)
  ID_receiver = n_nodes_l1 + if(i %% n_nodes_l2 == 0)  n_nodes_l2 else i %% n_nodes_l2
  weight = as.numeric(no_ch_funder[i])
  node = data.frame("ID1" = ID_sender, "ID2" = ID_receiver, "weight" = weight)
  connections = rbind(connections, node)
}

#H2H connections
#the offset from which we will start iterating in the weights vector
offset =  n_nodes_l1 * n_nodes_l2

for( i in 1:(n_nodes_l2 * n_nodes_l2))
{
  #add the number of nodes in previous layer to find ID
  ID_sender = n_nodes_l1 + 1 + ((i - 1) %/% n_nodes_l2)
  #same as before
  ID_receiver = n_nodes_l1 + if(i %% n_nodes_l2 == 0)  n_nodes_l2 else i %% n_nodes_l2
  #apply the offset to the iterator
  weight = as.numeric(no_ch_funder[offset + i ])
  node = data.frame("ID1" = ID_sender, "ID2" = ID_receiver, "weight" = weight)
  connections = rbind(connections, node)
}

#H2O connections
#calculate the offset from which we will start iterating in the weights vector
offset2 =  n_nodes_l1 * n_nodes_l2 + n_nodes_l2 * n_nodes_l2

for( i in 1:(n_nodes_l2 * n_nodes_l3))
{
  #add the number of nodes in previous layer to find ID
  ID_sender = n_nodes_l1 + 1 + ((i - 1) %/% n_nodes_l2)
  #same as before
  ID_receiver = n_nodes_l1 + n_nodes_l2 + if(i %% n_nodes_l3 == 0)  n_nodes_l3 else i %% n_nodes_l3
  #apply the offset to the iterator
  weight = as.numeric(no_ch_funder[offset2 + i ])
  node = data.frame("ID1" = ID_sender, "ID2" = ID_receiver, "weight" = weight)
  connections = rbind(connections, node)
}
