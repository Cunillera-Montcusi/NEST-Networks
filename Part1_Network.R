
# Hey!
# Welcome to the NEST workshop! 
# I am David Cunillera Montcusí and I am probably talking to you right now. So focus and look at me! 

# Network-focused mostly used packages 
library(sna)
library(igraph)
library(bipartite)

#* Bonus for some specific structures (Mimimum spanning tree)
library(ape)

# Network representation
library(ggnetwork)
library(gridExtra)

# Plot and data management
library(viridis)
library(tidyverse)

# We charge the example dataset
xy.Ponds <- read.csv2("data/Ponds.csv") %>% mutate(UTM_x=as.numeric(UTM_x),UTM_y=as.numeric(UTM_y))
xy.Ponds <- xy.Ponds[1:200,] # We select a smalles subset of 200 sites just as a matter of considering few 

### Distance matrix
Dist_matrix <- as.matrix(dist(xy.Ponds[,3:4]))
# For lat and long geosphere::distm()

# Percolation distance is your friend for local & regional levels. Calculated using max.comp_gradient function (see later)
Dist_percol <- ifelse(as.matrix(Dist_matrix)>2800,0,1) # We define the threshold distance at which to nodes will be considered as connected. 
diag(Dist_percol) <- 0 # We set the diagonal to 0 (distance to yourself is 0).

n<- network(Dist_percol, directed=F, diag=F) # One easy way to represent networks is the package network as it works directly wiht the matrix
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]), # Once creaated, the object "n" can be treated within a ggplot functions
       aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", linewidth=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue",fill="darkblue")+
  labs(x="",y="")+
  theme_void()


### MST network - Minimum spanning tree
# Not everything is summarised by the Percolation network. Different network structures can be defined. 
MST.Dist <- mst(Dist_matrix) # We create the matrix linking sites along the MST
n<- network(MST.Dist, directed=F, diag=F)
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(shape=21,alpha=.75,colour="blue")+
  labs(x="",y="")+
  theme_void()

# Packages and usage ####
# sna
# works directly with the matrix
Dist_percol
sna::degree(Dist_percol)

# igraph
g <- igraph::graph.adjacency(Dist_percol)
igraph::degree(g)

# Warning, if working with the two packages at the same time, make sure that functions correspond to one another. 

# Comparisons between networks ####
# We use diameter as a descriptor
# For MST 
# Now we use igraph, which requires an extra step to create a graph object
g <- igraph::graph.adjacency(MST.Dist) # Create a graph

igraph::diameter(g) # Calculate diameter
id.1 <- as.vector(igraph::get.diameter(g)) # Get the sites ID that define the diameter

# We create a variable to use for plotting defining if the site belongs or not to the diameter
Diam_ID <- rep("NO",nrow(xy.Ponds))
Diam_ID[id.1] <- "YES"

n<- network(MST.Dist, directed=F, diag=F)
n %v% "Diameter" <- Diam_ID
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Diameter),shape=21, alpha=.75)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",panel.background=element_blank())

# We repeat the same procedure but now using the percolation distance and compare the two outputs
g <- igraph::graph.adjacency(Dist_percol)
igraph::diameter(g)
id.2 <- as.vector(igraph::get.diameter(g))

Diam_ID <- rep("NO",nrow(xy.Ponds))
Diam_ID[id.2] <- "YES"

n<- network(Dist_percol, directed=F, diag=F)
n %v% "Diameter" <- Diam_ID
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Diameter),shape=21, alpha=.75)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",panel.background=element_blank())


# Network descriptors ####

# We will use the percolation network to calculate these values. We first define the graph object. 
g <- igraph::graph.adjacency(Dist_percol)

# Average path length 
igraph::average.path.length(g)

# Connectivity (relative fill of the network)
igraph::edge_density(g)

# Modularity ####
modul <- igraph::spinglass.community(g) # Different modularity algorithms exists to determine modules and network modularity 
# +Information: 

memb <- modul$membership

# Plot modules 
n<- network(Dist_percol, directed=F, diag=F)
n %v% "Memb_ID" <- LETTERS[memb] # We add the factor with the different groups 
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Memb_ID),shape=21, alpha=.75, size=4)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank())

# Let's test it with another dataset Predator-Prey network
load("data/Cor_troph.RData")
# Similarity between predators on consumed prays
cor.trof
d.t.aac<-as.matrix((1-vegan::vegdist(t(cor.trof), method = "bray", upper = TRUE)))
# We create the graph
g_trop <- igraph::graph_from_adjacency_matrix(d.t.aac, mode = "lower", weighted = TRUE, diag = FALSE)
# We assess the modularity (cluster_louvain allows to consider weights)
modul_trop<- igraph::cluster_louvain(g_trop, weights = igraph::E(g_trop)$weight)
memb_trop <- modul_trop$membership # We obtain the membership of each node to a determined module

# Let's plot it
net_trop<- network(d.t.aac, directed=F, diag=F)
net_trop %v% "Memb_ID" <- LETTERS[memb_trop]
ggplot(net_trop, aes(x, y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.4, color = "grey60", linewidth = (igraph::E(g_trop)$weight*5))  +
  geom_nodes(aes(fill=Memb_ID),shape = 21, size = 7, color = "black", alpha = 0.8) +
  scale_fill_viridis(discrete = T)+
  theme_blank()

# Backbone test
# g_trop_bkvn <- backbone::disparity(g_trop, alpha = 0.2, narrative = TRUE, class = "igraph")
# layout_g_trop_bkvn <- igraph::layout_with_fr(g_trop_bkvn)
# modul_trop_bkvn<- igraph::cluster_louvain(g_trop_bkvn, weights = igraph::E(g_trop_bkvn)$weight)
# igraph::V(g_trop_bkvn)$Module_Belonging <- modul_trop_bkvn$membership
# 
# ggnetwork::ggnetwork(igraph::simplify(g_trop_bkvn),layout=layout_g_trop_bkvn) %>% 
#   ggplot()+
#   ggnetwork::geom_edges(aes(x=x, xend=xend,y=y,yend=yend,alpha=weight),curvature = 0.1)+
#   ggnetwork::geom_nodes(aes(x=x,y=y,colour=Module_Belonging),size=3)+
#   scale_colour_viridis()+
#   theme_blank()


# Centrality metrics ####
# Calculated using sna package (based on matrix)

## Degree: number of links of each site
# UnWeigheted
sna::degree(Dist_percol,gmode="graph",cmode="freeman",ignore.eval=TRUE)
n %v% "Degree_Unw" <- sna::degree(Dist_percol,gmode="graph",cmode="freeman",ignore.eval=TRUE)

# Weighted degree
sna::degree((Dist_percol*Dist_matrix),gmode="graph", cmode="freeman",ignore.eval=FALSE)
n %v% "Degree_Weig" <- sna::degree((Dist_percol*Dist_matrix),gmode="graph", cmode="freeman",ignore.eval=FALSE)

# Plot the twor results one next to each other
gridExtra::grid.arrange(
ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Degree_Unw),shape=21, alpha=.75,size=5)+
  scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),

ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Degree_Weig),shape=21, alpha=.75,size=5)+
  scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
ncol=2)

# Betweenness
# Unweighted betweenness
sna::betweenness(Dist_percol,gmode="graph",cmode="undirected",ignore.eval=TRUE)
n %v% "Betwee_Unw" <- sna::betweenness(Dist_percol,gmode="graph",cmode="undirected",ignore.eval=TRUE)

# Unweighted betweenness
sna::betweenness((Dist_percol*Dist_matrix),gmode="graph",cmode="undirected",ignore.eval=FALSE)
n %v% "Betwee_Weig" <- sna::betweenness((Dist_percol*Dist_matrix),gmode="graph",cmode="undirected",ignore.eval=FALSE)

# Plot the twor results one next to each other
gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Betwee_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Betwee_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)


## Closenness
# Unweighted closenness
sna::closeness(Dist_percol,gmode="graph",ignore.eval=TRUE,cmode="undirected")
n %v% "Clos_Unw" <- sna::closeness(Dist_percol,gmode="graph",ignore.eval=TRUE,cmode="undirected")

# Weighted closenness
sna::closeness((Dist_percol*Dist_matrix),gmode="graph",ignore.eval=FALSE,cmode="undirected")
n %v% "Clos_Weig" <- sna::closeness((Dist_percol*Dist_matrix),gmode="graph",ignore.eval=FALSE,cmode="undirected")->clos_w

# Plot the twor results one next to each other
gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Clos_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=Clos_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)


##EIGENVECTOR centrality
# Unweighted closenness
sna::evcent(Dist_percol,gmode="graph",ignore.eval=TRUE)
n %v% "EgV_Unw" <- sna::evcent(Dist_percol,gmode="graph",ignore.eval=TRUE)

# Weighted closenness
sna::evcent(Dist_percol*Dist_matrix,gmode="graph",ignore.eval=TRUE)
n %v% "EgV_Weig" <- sna::evcent(Dist_percol*Dist_matrix,gmode="graph",ignore.eval=TRUE)

# Plot the twor results one next to each other
gridExtra::grid.arrange(
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=EgV_Unw),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  
  ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=EgV_Weig),shape=21, alpha=.75,size=5)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+theme_classic()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank()),
  ncol=2)


# Modularity and topological roles ####
# Bipartite is a package for network analysis that has some functions already available to identify topological roles-.
Mod.A <- bipartite::computeModules(Dist_percol)
ListMod <- bipartite::listModuleInformation(Mod.A)
bipartite::printoutModuleInformation(Mod.A)

roles_1 <- bipartite::czvalues(Mod.A, weighted=FALSE, level="lower")

# We plot the restults where we can identify different roles.
plot(roles_1$z~roles_1$c)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5)   # threshold of Olesen et al. 2007

# The same as before, we define the factors to plot and then plot the nework 
Memb_ID <- rep("NO", nrow(xy.Ponds))
Memb_ID[which(roles_1$c>0.62)] <- "Connector"
Memb_ID[max(roles_1$z)] <- "ModuleHub"

n<- network(Dist_percol, directed=F, diag=F)
n %v% "Memb_ID" <- Memb_ID

ggplot(n, layout=as.matrix(xy.Ponds[,3:4]),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges( color = "grey60", size=0.1, alpha=0.4) +
  geom_nodes(aes(fill=Memb_ID),shape=21, alpha=.75, size=4)+
  scale_fill_viridis(discrete = T,alpha = 1,begin = 1,end = 0)+
  labs(x="",y="")+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background=element_blank())



