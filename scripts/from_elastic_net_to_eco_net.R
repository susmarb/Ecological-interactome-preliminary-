#!/bin/R

# This script takes the coefficient tables from the elastic models and merge them to a single file in a SIF format (network)
# Then, if converst from directed network with + and - interactions to ecological network
# define arguments
args <- commandArgs(TRUE)
in_dir <- args[1]
omic <- args[2]
out_dir <- args[3]
win_y <- args[4]
pc_y <- args[5]
dir.create(out_dir, recursive = T, showWarnings = F)

packs <- c("tidyverse")
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}

# in_dir="~/Documents/LAO/LAO_DYNAMICS/elas_nets_mg"
# omic="mg"
# out_dir="~/Documents/LAO/LAO_DYNAMICS/networks"
# win_y="yes"
# pc_y="yes"


if(win_y == "no" & pc_y == "no"){
	enets_files <- list.files(in_dir, pattern = "coefficients.tsv", full.names = T)
}
if(win_y == "no" & pc_y == "yes"){
	enets_files <- list.files(in_dir, pattern = "coefficients_pc.tsv", full.names = T)
}

if(win_y == "yes" & pc_y == "no"){
	in_dir <- paste0(in_dir, "_timeWindows")
	enets_files <- list.files(in_dir, pattern = "coefficients.tsv", full.names = T)
}
if(win_y == "yes" & pc_y == "yes"){
	in_dir <- paste0(in_dir, "_timeWindows")
	enets_files <- list.files(in_dir, pattern = "coefficients_pc.tsv", full.names = T)
}


# from models to directed networks and save
if(win_y == "yes"){
	windows <- c("BF", "DR", "AF")
}else{windows <- ""}

for(win in windows){
	net_enets <- tibble(A = as.character(), B = as.character(), Sign = as.character(), W = as.double(), Rsq = as.double())
	enets_files_win <- enets_files[which(grepl(win, enets_files))]
for (file in enets_files_win){
	#print(file)
	enet_i <- read_tsv(file) %>% 
		filter(!grepl("ntercept", predictor_var)) %>% 
		mutate(Sign = ifelse(coefficient < 0, "Negative", "Positive")) %>% 
		rename(A = predictor_var, B = response_var, W = coefficient, Rsq = rsq_model) %>% 
		select(A, B, Sign, W, Rsq)
	net_enets <- rbind(net_enets, enet_i)
}

# histogram with Rsq distribution and its linked table
net_enets %>% select(B, Rsq) %>% unique() -> enet_rsq
enet_rsq %>% 
	ggplot(aes(x=Rsq)) +
	geom_histogram(color="black", fill="white", bins = 30) +
	theme_bw() -> hist_rsq

# filter interactions from models with Rsq above 0.5
net_enets %>% 
	filter(Rsq >= 0.5) -> net_enets

# save tables and plot
if(pc_y == "no"){
	write_tsv(enet_rsq, paste0(out_dir, "/", omic, win, "_enet_rsq.tsv"))
	pdf(paste0(out_dir, "/", omic, win, "_enet_histRsq.pdf"), height = 6, width = 6, useDingbats = T); print(hist_rsq); dev.off()
	write_tsv(net_enets, paste0(out_dir, "/", omic, win, "_enet.tsv"))
	}
if(pc_y == "yes"){
	write_tsv(enet_rsq, paste0(out_dir, "/", omic, win, "_enet_pc_rsq.tsv"))
	pdf(paste0(out_dir, "/", omic, win, "_enet_pc_histRsq.pdf"), height = 6, width = 6, useDingbats = T); print(hist_rsq); dev.off()
	write_tsv(net_enets, paste0(out_dir, "/", omic, win, "_enet_pc.tsv"))
}

# make ecological networks and save
net_enets %>% 
	select(-Rsq) %>% 
	rename(W_AB = W) %>% 
	mutate(interaction_AB = Sign) -> net_enets

net_enets_eco <- tibble(A=as.character(), B=as.character(),Sign=as.character(),
												W_AB=as.numeric(),Rsq=as.numeric(),interaction_AB=as.character(),
												interaction_BA=as.character(),
												W_BA=as.numeric(), interaction=as.character(),
												W=as.numeric())
for(i in 1:length(net_enets$A)){
	a <- net_enets$A[i]
	b <- net_enets$B[i]
	
	tmp_a <- net_enets %>% filter(A==b, B==a)
	# remove this line from the net_enets, to avoid duplicate relationships
	index2remove_i <- which(net_enets$A==b & net_enets$B==a)
	
net_enets_i <- net_enets %>% 
		slice(i) %>% 
		mutate(interaction_BA = tmp_a$Sign[1]) %>% 
		mutate(W_BA = tmp_a$W_AB[1]) %>% 
		mutate(interaction = case_when(
			interaction_AB == "Positive" & interaction_BA == "Positive" ~ "cooperation", 
			interaction_AB == "Positive" & interaction_BA == "Negative" ~ "predation",
			interaction_AB == "Negative" & interaction_BA == "Positive" ~ "predation",
			interaction_AB == "Negative" & interaction_BA == "Negative" ~ "competition",
			is.na(interaction_AB) & interaction_BA == "Negative" ~ "amensalism",
			interaction_AB == "Negative" & is.na(interaction_BA) ~ "amensalism",
			is.na(interaction_AB) & interaction_BA == "Positive" ~ "comensalism",
			interaction_AB == "Positive" & is.na(interaction_BA) ~ "comensalism")) %>% 
		mutate(W = case_when(
			!is.na(W_AB) & !is.na(W_BA) ~ abs(W_AB) + abs(W_BA),
			!is.na(W_AB) & is.na(W_BA) ~ abs(W_AB),
			is.na(W_AB) & !is.na(W_BA) ~ abs(W_BA)))

if(is_empty(index2remove_i)==F){net_enets <- net_enets[-index2remove_i,]}
net_enets_eco <- rbind(net_enets_eco, net_enets_i)
rm(net_enets_i, tmp_a, a, b, index2remove_i)
}
ifelse(pc_y == "no",
			 write_tsv(net_enets_eco, paste0(out_dir, "/", omic, win, "_enet_eco.tsv")),
			 write_tsv(net_enets_eco, paste0(out_dir, "/", omic, win, "_enet_eco_pc.tsv")))
}
