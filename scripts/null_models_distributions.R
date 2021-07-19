#!/bin/R

# Testing null models

# define arguments
args <- commandArgs(TRUE)
data_in <- args[1]

packs <- c("tidyverse", "sna", "igraph", "signnet")
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}

# data_in <- "~/Documents/LAO/LAO_DYNAMICS/networks/mg_enet_eco_pc.tsv"

net_tb <- read_tsv(data_in)

#for initial network of raw model results, save only basic numbers
if(grepl("eco", data_in) == F){
	# net tb to igraph object
	edges <- net_tb %>% rename(from = A, to = B) %>% mutate(W = abs(W))
	nodes <- tibble(nodes = unique(c(edges$from, edges$to)))
	## create network object, igraph object
	g <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
	g <- set_edge_attr(g, "weight", value= edges$W)
	#is_weighted(g)
}else{
	net_tb %>% 
		mutate(Sign = case_when(interaction == "competition" ~ -1,
														interaction == "amensalism" ~ -1,
														interaction == "cooperation" ~ 1,
														interaction == "comensalism" ~ 1,
														interaction == "predation" ~ -1)) -> net_tb
	# net tb to igraph object
	edges <- net_tb %>% rename(from = A, to = B) %>% mutate(W = abs(W))
	nodes <- tibble(nodes = unique(c(edges$from, edges$to)))
	## create network object, igraph object
	g <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
	g <- set_edge_attr(g, "weight", value= edges$W)
	g <- set_edge_attr(g, "sign", value= edges$Sign)
}

# Erdos-Renyi, Barabasi-Albert, and Stochastic-Block models
n <- vcount(g)
p <- 0.05
k_mean_erdos <- p*(n-1)
k_mean_bar <- mean(degree(g))
k_mean1_stock <- 0.3*(30-1)
k_mean2_stock <- 0.5*(41-1) 

prob_degrees_tb <- tibble(k = 1:n, 
													prob_observed = degree.distribution(g),
													prob_erdos = dpois(x=1:n,lambda=k_mean_erdos),
													prob_barabasi = dpois(x=1:n,lambda=k_mean_bar),
													prob_stock_1 = dpois(x=1:n,lambda=k_mean1_stock),
													prob_stock_2 = dpois(x=1:n,lambda=k_mean2_stock)) #%>% 


ggplot(prob_degrees_tb,aes(x=k)) +
	geom_point(aes(y=prob_observed), colour="black", shape=1) +
	geom_line(aes(y=prob_erdos), colour="blue") +
	geom_line(aes(y=prob_barabasi), colour="green") +
	geom_line(aes(y=prob_stock_1), colour="purple") +
	geom_line(aes(y=prob_stock_2), colour="purple") +
	theme_classic() +
	theme(legend.position = "right",
				axis.text.x = element_text(hjust = 0.95, vjust = 0.2, size = 12),
				axis.text.y = element_text(size = 12)) +
	xlab("Nodes degree distribution") +
	ylab("Probabilities") -> plot_tests; plot_tests
pdf(paste0(data_in, "_degre_distributions_null_models.pdf"), width = 16/2.56, height = 16/2.56, useDingbats = T); print(plot_tests); dev.off()
