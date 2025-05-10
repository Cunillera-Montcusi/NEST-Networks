

##########  Package uploading   ##########
# Network representation
library(network)
library(ggnetwork)

# Geographic tools (distances based on lat/long data) 
library(geosphere)

# Plot and data management
library(viridis)
library(tidyverse)

# Network analysis
library(igraph)
library(sna)

##########    Data uploading    ##########
## Environmental data
load("data/env.Rdata")
env 

## Sites coordinates 
Netowrks <- load("data/Networks.RData") 
Lake_Coordinates <- Netowrks[[3]]# We select one specific scale number 3

# We define which sites correspond to which sites on the network
coincidence <- (nrow(Lake_Coordinates)-54):nrow(Lake_Coordinates)

########  Percolation network definition   ##########
Dist_matrix <- geosphere::distm(Lake_Coordinates) # Geosphere for lat/long data

# How to find the percolation distance? 
# Many threshold distances from which we ask the number of components (connected neighbors), the first distance that contains all neighbours is the percolation distance.
max.comp_gradient<-function(Dist_mat, precision){
  Range<-seq(min(Dist_mat),max(Dist_mat),,precision)
  out<-NULL 
  for(t in Range){
    C<-sna::component.dist(ifelse(Dist_mat>t,0,1))
    C.max<-max(C$csize)
    out<-rbind(out,c(t,C.max))
  }
  
  out <- out[which(out[,2]==ncol(Dist_mat))[1],]
  return(out)
}

# We use it for our casa
max.comp_gradient(Dist_mat = Dist_matrix,precision = 100)
# Obtained distance 72036.2  
Dist_percol <- ifelse(as.matrix(Dist_matrix)>72036.2,0,1)
diag(Dist_percol) <- 0 # We eensure the diagonal is 0 (distance from one site to itself is 0)

n<- network::network(Dist_percol, directed=F, diag=F) # One easy way to represent networks is the package network as it works directly wiht the matrix
ggplot(n, layout=as.matrix(Lake_Coordinates), # Once creaated, the object "n" can be treated within a ggplot functions 
       aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", linewidth=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue",fill="darkblue")+
  labs(x="",y="")+
  theme_void()


##########  Threshold distance   ##########

# Threshold distance
Distances <- c(10000,45000,72036.2,75000,120000)

# Now we will test different network threshold distances and we will assess their impacts on centrality. 
# Instead of using directly the distance matrix (package sna), we will create a graph object (package igraph) 

Clo_Netw <- list()
for (Dist_Netw in 1:length(Distances)) { # This is a loop to create a graph for each distance 
  Dist_percol_temp <- ifelse(as.matrix(Dist_matrix)>Distances[[Dist_Netw]],0,1) # We define the connectivity matrix 
  diag(Dist_percol) <- 0 # Make sure that the diagonal is 0
  Clo_Netw[[Dist_Netw]] <- igraph::graph.adjacency(Dist_percol_temp, mode = "undirected",diag = F) # We create a graph object from which to calculate centrality
}

# We create a data frame with where we calculate the centrality for all nodes of the different network, each created with a different threshold distance
Dataset_all <- data.frame("Clo_1000"=igraph::harmonic_centrality(Clo_Netw[[1]], mode = "all",normalized = T),
                          "Clo_45000"=igraph::harmonic_centrality(Clo_Netw[[2]], mode = "all",normalized = T),
                          "Clo_72000"=igraph::harmonic_centrality(Clo_Netw[[3]], mode = "all",normalized = T),
                          "Clo_75000"=igraph::harmonic_centrality(Clo_Netw[[4]], mode = "all",normalized = T),
                          "Clo_120000"=igraph::harmonic_centrality(Clo_Netw[[5]], mode = "all",normalized = T)
  ) %>% gather(key, value)%>% mutate(key=factor(key,levels=c("Clo_1000","Clo_45000","Clo_72000","Clo_75000","Clo_120000"))) 

# We do the same but now for each of the sampled lakes using "vids" we can calculate centrality for specific nodes 
Dataset_sites <- data.frame("Clo_1000"=igraph::harmonic_centrality(Clo_Netw[[1]], mode = "all",vids = coincidence,normalized = T),
                            "Clo_45000"=igraph::harmonic_centrality(Clo_Netw[[2]], mode = "all",vids = coincidence,normalized = T),
                            "Clo_72000"=igraph::harmonic_centrality(Clo_Netw[[3]], mode = "all",vids = coincidence,normalized = T),
                            "Clo_75000"=igraph::harmonic_centrality(Clo_Netw[[4]], mode = "all",vids = coincidence,normalized = T),
                            "Clo_120000"=igraph::harmonic_centrality(Clo_Netw[[5]], mode = "all",vids = coincidence,normalized = T)
) %>% gather(key, value) %>% mutate(key=factor(key,levels=c("Clo_1000","Clo_45000","Clo_72000","Clo_75000","Clo_120000"))) 

# We plot the histogram of centrality values 
Plot_Threshold <- Dataset_all %>% 
  ggplot()+
  geom_histogram(aes(x=value), bins=100,colour="black",fill="grey10")+
  geom_vline(data = Dataset_sites, aes(xintercept =value) , colour="grey90",alpha=0.1)+
  theme_classic()+
  facet_grid(key~ .,scales = "free") 


# Scale

# We will repeat the same process with different networks that are basically making the distance buffer larger.
# This implies that for each network we are increasing the number of lakes that we consider

Clo_Netw <- list()
Netw_plot <- list()
Clo_Coinci <- list()
for (Dist_Netw in 1:5) {
  # We repeat the same process that we did at the beginning
  Lake_Coordinates <- Netowrks[[Dist_Netw]] # We select the coordinates of the different networks
  Dist_matrix <- geosphere::distm(Netowrks[[Dist_Netw]]) # Calculate the distance matrix for each network
  Dist_percol_temp <- ifelse(as.matrix(Dist_matrix)>100000,0,1) # Set the distance from which they are connected, the same for all now
  diag(Dist_percol_temp) <- 0 # Ensure diagonal is 0
  
  Clo_Netw[[Dist_Netw]] <- igraph::graph.adjacency(Dist_percol_temp, mode = "undirected",diag = F) # We create a graph object to calculate centrality
  
  # Some of the networks are big and take a bit of time depending on the computer. This part will create a network plot for each of them if activated.
  #n<- network::network(Dist_percol_temp, directed=F, diag=F)
  #Netw_plot[[Dist_Netw]]<-ggplot(n, layout=as.matrix(Lake_Coordinates),aes(x = x, y = y, xend = xend, yend = yend))+
  #                               geom_edges( color = "grey60", linewidth=0.1, alpha=0.4) +
  #                               geom_nodes(shape=21,alpha=.75,colour="blue",fill="darkblue")+
  #                               labs(x="",y="")+
  #                               theme_void()
  
  # We store the position of the sampled lakes
  Clo_Coinci[[Dist_Netw]]<- (nrow(Lake_Coordinates)-54):nrow(Lake_Coordinates)
  
}

# We create a data frame with where we calculate the centrality for all nodes of the different network, each created from a different scale
Dataset_all <- bind_rows(data.frame(key="Clo_Smallest","value"=igraph::harmonic_centrality(Clo_Netw[[5]], mode = "all",normalized = T)),
                         data.frame(key="Clo_Small","value"=igraph::harmonic_centrality(Clo_Netw[[4]], mode = "all",normalized = T)),
                         data.frame(key="Clo_Mid","value"=igraph::harmonic_centrality(Clo_Netw[[3]], mode = "all",normalized = T)),
                         data.frame(key="Clo_Large","value"=igraph::harmonic_centrality(Clo_Netw[[2]], mode = "all",normalized = T)),
                         data.frame(key="Clo_Largest","value"=igraph::harmonic_centrality(Clo_Netw[[1]], mode = "all",normalized = T))
  ) %>% mutate(key=factor(key,levels=c("Clo_Smallest","Clo_Small","Clo_Mid","Clo_Large","Clo_Largest")))

# We do the same but now for each of the sampled lakes using "vids" we can calculate centrality for specific nodes 
Dataset_sites <- data.frame("Clo_Smallest"=igraph::harmonic_centrality(Clo_Netw[[5]], mode = "all",vids = Clo_Coinci[[5]],normalized = T),
                            "Clo_Small"=igraph::harmonic_centrality(Clo_Netw[[4]], mode = "all",vids = Clo_Coinci[[4]],normalized = T),
                            "Clo_Mid"=igraph::harmonic_centrality(Clo_Netw[[3]], mode = "all",vids = Clo_Coinci[[3]],normalized = T),
                            "Clo_Large"=igraph::harmonic_centrality(Clo_Netw[[2]], mode = "all",vids = Clo_Coinci[[2]],normalized = T),
                            "Clo_Largest"=igraph::harmonic_centrality(Clo_Netw[[1]], mode = "all",vids = Clo_Coinci[[1]],normalized = T)
) %>% gather(key, value) %>% mutate(key=factor(key,levels=c("Clo_Smallest","Clo_Small","Clo_Mid","Clo_Large","Clo_Largest")))

# We plot the histogram of centrality values 
Plot_scale <- Dataset_all %>% 
  ggplot()+
  geom_histogram(aes(x=value), bins=100,colour="black",fill="grey10")+
  geom_vline(data = Dataset_sites, aes(xintercept =value) , colour="grey90",alpha=0.1)+
  theme_classic()+
  facet_grid(key~ .,scales = "free") 


Plot_Threshold

Plot_scale


