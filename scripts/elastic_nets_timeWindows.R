#!/bin/R

# This script performs elastic net modeling using the caret package

# define arguments
args <- commandArgs(TRUE)
data_in <- args[1]
out_dir <- args[2]
src_dir <- args[3]
pcs <- args[4]
pc_in <- args[5]
mod_all <- args[6]
mod_this <- args[7]


packs <- c("tidyverse", "caret", "glmnet", "imputeTS")
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}

#  data_in="~/Documents/LAO/Tables/mg_avg_depth_filt_relab_clr.tsv"
#  out_dir="~/Documents/LAO/LAO_DYNAMICS/elas_nets_MG_timeWindows"
#  src_dir="~/git/dynamics-characterization-lao/src"
#  pcs="no"
#  pc_in="~/Documents/LAO/Tables/pc_params.processed.tsv"
#  mod_all="no"
#  mod_this="D51_G1.1.2"

dir.create(out_dir, recursive = T, showWarnings = F)

funcs_doc <- paste0(src_dir, "/functions.R")
source(funcs_doc)

# Read data table
mydata_all <- read_tsv(data_in)

if(pcs == "yes"){
	pc_tb <- read_tsv(pc_in) %>% select(-Type) %>% 
		spread(variable, value) %>% 
		imputeTS::na_interpolation(option = "stine")
	mydata_all <- left_join(mydata_all, pc_tb)
}

# Set training control
myTimeControl <- trainControl(method = "timeslice",
	initialWindow = 10, # initial number of consecutive values in each training set sample
	horizon = 5, # number of consecutive values in test set sample
	fixedWindow = TRUE, # if FALSE, the training set always start at the first sample and the training set size will vary over data splits
	verboseIter = TRUE,
	number = 10,
	search = "grid")
# Set grid of alpha and lambda (which are the parameters chosen by cross-validation)
# by seleting the pair tht minimises the cv error
lambda_grid <- 10^seq(2,-2,length=100)
alpha_grid <- seq(0,1, length=10)
srchGrd <- expand.grid(.alpha=alpha_grid, .lambda= lambda_grid)

# if I only want to model one MAG, loop is not necesary
# Variables vector (mydata has first column with dates, and the rest with MAGs)
variables_mags <- colnames(mydata_all)[-1]
if(mod_all == "no"){variables_mags = mod_this}

# loop through windows (3 overlapping windows of 20 timepoints each)
windows <- c("BF", "DR", "AF")

for(win in windows){

if(win == "BF"){mydata <- mydata_all[1:20,]}
if(win == "DR"){mydata <- mydata_all[17:36,]}
if(win == "AF"){mydata <- mydata_all[32:51,]}

# Train the model (per MAG, so here the loop starts)
for (res_var in variables_mags){
	formula_var <- formula(paste0(res_var, " ~ ."))
	set.seed(666)
	plsFitTime <- train(formula_var,
			data = mydata[,-1], # it keeps all the MAGs to make feature selection through the lasso of the elastic net
			method = "glmnet",
			tuneGrid = srchGrd,
			tuneLength = 10,
			trControl = myTimeControl,
			verbose = FALSE)

# multiple R-squared
## calculate predicted values with the model
pred_enet <- predict(plsFitTime, mydata)
## calculate the multiple R squared
rsq <- cor(mydata[,which(colnames(mydata) == res_var)], pred_enet)^2
# save coefficients of the best model
coef_mat <- as.matrix(coef(plsFitTime$finalModel, plsFitTime$bestTune$lambda))
coef_tb <- tibble(response_var = res_var,
									predictor_var = rownames(coef_mat), 
									coefficient = coef_mat[,1], 
									rsq_model = rsq) %>%
						filter(coefficient !=0)
ifelse(pcs == "yes", 
			 write_tsv(coef_tb, paste0(out_dir, "/", win, "_", res_var, "_elasNet_coefficients_pc.tsv")),
			 write_tsv(coef_tb, paste0(out_dir, "/", win, "_", res_var, "_elasNet_coefficients.tsv")))

# save residuals of best model
resid_mode <- tibble(response_var = res_var, date = mydata$Date, residuals = residuals(plsFitTime))
ifelse(pcs == "yes", 
			 write_tsv(resid_mode, paste0(out_dir, "/", win, "_", res_var, "_elasNet_residuals_pc.tsv")), 
			 write_tsv(resid_mode, paste0(out_dir, "/", win, "_", res_var, "_elasNet_residuals.tsv")))

rm(formula_var, plsFitTime, pred_enet, rsq, coef_mat, coef_tb, resid_mode)
}
}
