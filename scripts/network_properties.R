#!/bin/R

# Network topology

# define arguments
args <- commandArgs(TRUE)
data_in <- args[1]

packs <- c("tidyverse", "sna", "igraph", "signnet")
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}

# data_in <- "~/Documents/LAO/LAO_DYNAMICS/networks/mg_enet_eco_pc.tsv"

net_tb <- read_tsv(data_in)

#for initial network of raw model results, save only basic numbers
if(grepl("eco", data_in) == F){
# percentage of positive and negative interactions
	net_tb %>% 
		group_by(Sign) %>% tally() %>% 
		mutate(percentage = n*100/sum(n)) -> sign_n
	write_tsv(sign_n, paste0(file_in, "_sign_n.tsv"))

	# net tb to igraph object
	edges <- net_tb %>% rename(from = A, to = B) %>% mutate(W = abs(W))
	nodes <- tibble(nodes = unique(c(edges$from, edges$to)))
	## create network object, igraph object
	g <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
	g <- set_edge_attr(g, "weight", value= edges$W)
	#is_weighted(g)
}	
if(grepl("eco", data_in) == T){
	net_tb %>% 
		mutate(Sign = case_when(interaction == "competition" ~ -1,
														interaction == "amensalism" ~ -1,
														interaction == "cooperation" ~ 1,
														interaction == "comensalism" ~ 1,
														interaction == "predation" ~ -1)) -> net_tb
	# percentage of types of interactions
	net_tb %>% 
		group_by(interaction) %>% tally() %>% 
		mutate(percentage = n*100/sum(n)) -> sign_n
	write_tsv(sign_n, paste0(data_in, "_sign_n.tsv"))
	
	
	# net tb to igraph object
	edges <- net_tb %>% rename(from = A, to = B) %>% mutate(W = abs(W))
	nodes <- tibble(nodes = unique(c(edges$from, edges$to)))
	## create network object, igraph object
	g <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
	g <- set_edge_attr(g, "weight", value= edges$W)
	g <- set_edge_attr(g, "sign", value= edges$Sign)
	#is_weighted(g)
}
	# network topology metrics
net_metrics <- tibble(nodes = vcount(g),
											edges = ecount(g),
											density = graph.density(g),
											diameter = diameter(g),
											transitivity = transitivity(g))
write_tsv(net_metrics, paste0(data_in, "_net_metrics.tsv"))

obj <- degree(g)
centralities <- tibble(node = names(obj),
											degree = degree(g),
											degree_w = graph.strength(g),
											betweenness = betweenness(g),
											eigen_centrality = eigen_centrality(g)$vector,
											subgraph_centrality = subgraph_centrality(g))
											

if(grepl("eco", data_in) == T){
# save cluster results
	cluster = cluster_louvain(g)
	# nodes membership, add to centralities table
	clst_tab_nodes = tibble(node = names(obj),
				 clstr_member = cluster$membership)
	centralities <- left_join(centralities, clst_tab_nodes, by = "node")
	# network cluster metrics, add to network metrics table
	clst_tab_net <- tibble(modularity = modularity(cluster))
	net_metrics <- cbind(net_metrics, clst_tab_net )
	
	# do blockmodelling calculations, several blocks
	# alpha = 0.5, equal penalization of negative inter group and positive intra group edges
	for (i in 5:31){
		model_blocks_i <- signed_blockmodel(g,k = i,alpha = 0.5, annealing = TRUE)
		
		# node membership info
		model_blocks_i_tb <- tibble(node = names(obj), block_member = model_blocks_i$membership)
		colnames(model_blocks_i_tb) <- c("node", paste0("block_member_N", i))
		centralities <- left_join(centralities, model_blocks_i_tb, by = "node")
		
		#visualization
		blocks_plot_i <- ggblock(g,model_blocks_i$membership,show_blocks = TRUE) + 
		theme(axis.text.y.left = element_text(size = 10, family = "Helvetica"),
					axis.text.x.bottom = element_text(angle = 90, size = 10, family = "Helvetica"))
		# save plot
		pdf(paste0(data_in, "_blocks_N", i, ".pdf"), height = 20, width = 20); print(blocks_plot_i); dev.off()
	}

}
# save tables
write_tsv(centralities, paste0(data_in, "_node_centralities.tsv"))
write_tsv(net_metrics, paste0(data_in, "_network_topology_metrics.tsv"))

