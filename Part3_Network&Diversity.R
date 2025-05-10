
# Hey!
# Welcome to the NEST workshop morning session 2! 
# I am still David Cunillera Montcusí and I am probably talking to you right now. So focus and look at me! 

# Package sna is one of the two packages that we can use to calculate network metrics and create them
library(sna)
library(bipartite)
library(igraph)

# Network representation
library(ggnetwork)
library(gridExtra)

# Plot and data management
library(viridis)
library(tidyverse)

# We charge the example dataset (we will use the "Ponds" examaple again)
xy.Ponds <- read.csv2("data/Ponds.csv") %>% mutate(UTM_x=as.numeric(UTM_x),UTM_y=as.numeric(UTM_y))
xy.Ponds <- xy.Ponds[1:200,] # We select a smalles subset of 200 sites just as a matter of considering few 

### Distance matrix
Dist_matrix <- as.matrix(dist(xy.Ponds[,3:4])) # As they are UTM no need to geosphere
# For lat and long geosphere::distm()

# Percolation distance is your friend for local & regional levels
Dist_percol <- ifelse(as.matrix(Dist_matrix)>2800,0,1) # Calculated using max.comp_gradient function
diag(Dist_percol) <- 0

n<- network(Dist_percol, directed=F, diag=F)
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", linewidth=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue",fill="darkblue")+
  labs(x="",y="")+
  theme_void()


# Coalescent simulations
### We begin by preparing all the data that we need to run the simulation

# Function to coalescent building
source("data/function_to_coalescent.R")

# This function has 3 major elements to consider:
# 1) Distance matrix & distance matrix
# 2) Communities size
# 3) Environmental filters

# We are going to only modify the first element along this exercise but first, we need to define these elements. 

## 1) Distance matrix
# The model considers distances as a negative kernel of probability.
# So, instead of changing distances above percolation to "0" as before we artificially increase the distance of those sites.
M.distance=ifelse(as.matrix(Dist_matrix)>2800,500000,Dist_matrix)

# Dispersal distance is considered here as the distance at which probability of dispersal is 0.5. 
# So "higher" or "lower" do not exclude other dispersal abilities. 
# They push "overall connectivity" towards higher or lower connections.
# We set the dispersal ability at the percolation distance.
dispersal <- 2800 # In our exmple, we set  dispersal  at the percolation distance

## 2) Community size
# In our example all sites will have the same size
J.habitat <- rep(100,nrow(Dist_matrix))

## 3) Enviornmental filter is set to the same for all species
# This table summarises the interaction site x species witha "performance" value 
# which defines how "good" (1) or "bad" (0) is the performance of that species in that site
Species_Filter <- matrix(nrow=200,ncol = nrow(Dist_matrix), data = 0.99)

# Other parameters of the model
id_NOmodule <- rep(1,nrow(Dist_matrix)) # Modules if we want some sites to belong to the same module. 
pool_200 <- rep(1,nrow(Species_Filter)) # Distribution of the species pool #rlnorm(n = 200,5,1) 


# We create an output objects for each iteration
Output <- data.frame(); a <- NULL 
for (it in 1:3) { # We repeat 3 times this function to obtain 3 communities that we will later average
  output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200, # Regional species pool 
                                    Js=J.habitat, # Habitat size
                                    filter.env=Species_Filter, # Environmental filters
                                    M.dist=M.distance,  # Distance matrix
                                    D50=dispersal, # Dispersal 
                                    id.obs=1:nrow(M.distance)) # Parameter to focus on determined sites
a <- rbind(a,output[[1]])
}
Final_out <- resume.out(a) # Function to collapse and obtain results from eaach iteration

# We create the graph object to calculate centrality metrics
g <- igraph::graph.adjacency(Dist_percol)

out <- data.frame("Disp"=round(as.numeric(dispersal),0),
                  "Clo"=igraph::closeness(g,normalized = T), # Graph corresponding to the network
                  "S_Sim"=Final_out$Median[6:(nrow(Dist_percol)+5)], # Median richness of the iterations
                  "B_Sim"=Final_out$Median[(nrow(Dist_percol)+6):((nrow(Dist_percol)*2)+5)])# Median beta of the iterations
out %>%  
  ggplot(aes(x=Clo,y=S_Sim))+
  geom_point()+
  geom_smooth(method = "loess")+
  theme_classic()


# Dispersal driven patterns ####
# What if we cange the dispersal threshold? 

# We will repeat the same process for different dispersal abilities in order to capture regional structure effect on diversity patterns.
Output <- data.frame()
dispersal <- c(10,100,1000,2800,5000,10000,100000)
for (disp in 1:length(dispersal)) { # Loop repeaating this process for each dispersal ability
  a <- NULL # We create an output object for each iteration
  for (it in 1:3) { # We repeat 3 times the same process
    output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200, 
                                      Js=J.habitat, 
                                      filter.env=Species_Filter,
                                      M.dist=M.distance, 
                                      D50=dispersal[disp],
                                      id.obs=1:nrow(M.distance))
    a <- rbind(a,output[[1]])
    }
  Final_out <- resume.out(a)
  
  # We create the graph object to calculate centrality metrics
  g <- igraph::graph.adjacency(Dist_percol)
  
  out <- data.frame("Disp"=round(as.numeric(dispersal[disp]),0),
                    "Clo"=igraph::closeness(g,normalized = T),
                    "S_Sim"=Final_out$Median[6:(nrow(Dist_percol)+5)],
                    "B_Sim"=Final_out$Median[(nrow(Dist_percol)+6):((nrow(Dist_percol)*2)+5)])
  Output <- bind_rows(Output,out) 
}

# We plot and represent the obtained results
n<- network(Dist_percol, directed=F, diag=F) # Create a network
n %v% "Rich" <- Output %>% filter(Disp==2800) %>% pull(S_Sim) # Add factor corresponding to richness

# Plot the three plots together
gridExtra::grid.arrange(
  Output %>%  
    ggplot(aes(x=Disp,y=S_Sim,colour=Clo))+
    geom_jitter()+
    scale_x_log10()+scale_color_viridis(option="H")+scale_y_continuous(limits=c(0,65))+
    theme_classic(),
  
  Output %>%  
    ggplot(aes(x=Disp,y=B_Sim,colour=Clo))+
    geom_jitter()+
    scale_x_log10()+scale_color_viridis(option="H")+
    theme_classic(),

  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),
       aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( colour="grey60", linewidth=0.1, alpha=0.4) +
  geom_nodes(aes(color = Rich),shape=16,size=3)+
  scale_colour_viridis()+
  labs(x="",y="")+
  theme_void(),

ncol=2,heights=c(1,2))


# Landscape Loss ####
# WHAT IF? 
# Our landscape losses some sites and we want to assess the structural impacts that this would have on diversity?
# We lose sites and edit the corresponding datasets
lost_sites <- sample(1:200,size = 60)
Dist_percol_loss <- Dist_percol[-lost_sites,-lost_sites]
M.distance_loss <- M.distance[-lost_sites,-lost_sites]

# We repeat the exact same exercise but with the new network structure thaat has lost some sites
Output_loss <- data.frame()
dispersal <- c(10,100,1000,2800,5000,10000,100000)
for (disp in 1:length(dispersal)) {
  a <- NULL # We create an output object for each iteration
  for (it in 1:3) { # We repeat 3 times the same process
    output <- Coalescent.exp.Kernel.J(Meta.pool=pool_200,m.pool=0.001, 
                                      Js=J.habitat[-lost_sites], 
                                      filter.env=Species_Filter[,-lost_sites],
                                      M.dist=M.distance_loss, 
                                      D50=dispersal[disp],m.max=1,
                                      id.obs=1:nrow(M.distance_loss))
    a <- rbind(a,output[[1]])}
  Final_out <- resume.out(a)
  
  # We create the graph object to calculate centrality metrics
  g <- igraph::graph.adjacency(Dist_percol_loss)
  
  out <- data.frame("Disp"=round(as.numeric(dispersal[disp]),0),
                    "Clo_Loss"=igraph::closeness(g,normalized = T),
                    "S_Sim_Loss"=Final_out$Median[6:(nrow(Dist_percol_loss)+5)],
                    "B_Sim_Loss"=Final_out$Median[(nrow(Dist_percol_loss)+6):((nrow(Dist_percol_loss)*2)+5)])
  Output_loss <- bind_rows(Output_loss,out) 
}

n<- network(Dist_percol_loss, directed=F, diag=F)
n %v% "Rich_loss" <- Output_loss %>% filter(Disp==2800) %>% pull(S_Sim_Loss)

# Plot the three plots together
gridExtra::grid.arrange(
  Output_loss %>%  
    ggplot(aes(x=Disp,y=S_Sim_Loss,colour=Clo_Loss))+
    geom_jitter()+
    scale_x_log10()+scale_color_viridis(option="H")+
    scale_y_continuous(limits=c(0,65))+
    theme_classic(),
  
  Output_loss %>%  
    ggplot(aes(x=Disp,y=B_Sim_Loss,colour=Clo_Loss))+
    geom_jitter()+
    scale_x_log10()+scale_color_viridis(option="H")+
    theme_classic(),
  
  ggplot(n, layout=as.matrix(xy.Ponds[-lost_sites,3:4]),
         aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( colour="grey60", linewidth=0.1, alpha=0.4) +
    geom_nodes(aes(color = Rich_loss),shape=16,size=3)+
    scale_colour_viridis()+
    labs(x="",y="")+
    theme_void(),
  
  ncol=2,heights=c(1,2))



# Scenario comparisons - LogRatios ####
Original <- Output %>% filter(Disp==2800) %>% dplyr::select(-Disp) # We isolate the percolation distance for the original site
Lost <- Output_loss %>% filter(Disp==2800)%>% dplyr::select(-Disp) # We isolate the percolation distance for the loss site

n<- network(Dist_percol_loss, directed=F, diag=F) # We create a network with the information 
n %v% "Log_Ratio" <- bind_cols(Original[-lost_sites,],Lost) %>% 
                     mutate(Log_Ratio=log(S_Sim_Loss/S_Sim),2) %>% # We define the logRatio between the loss scenario and the S_Sim scenario
                     pull(Log_Ratio)

ggplot(n, layout=as.matrix(xy.Ponds[-lost_sites,3:4]),
       aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(colour="grey60", linewidth=0.1) +
  geom_nodes(aes(color = Log_Ratio),shape=16,size=3)+
  scale_color_gradient2(low = "white",mid = "orange",high = "red",midpoint = 0.2)+
  labs(x="",y="")+
  theme_void()




