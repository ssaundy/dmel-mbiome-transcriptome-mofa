rm(list=ls())

library(pheatmap)
library(data.table)
library(MOFA2)

	#threshold for quantile of coefficient of variation (transcriptome)
quantThresh <- 0.8

	#percentage (microbiome) samples in which you require data
percThresh <- 80	#80 good

	#samples to exclude
		#(for one reason and another...)
#excl <- c("DS37", "DS2")

##################
##	read data	##
##################

	#microbiota
#micro <- read.table("noWD_featuretable_L6_1_10.txt", sep="\t", header=T)
#micro <- read.table("feature_table_full.txt", sep="\t", header=T)
micro <- read.table("noWD_featuretable_L5_relativeLog2Abundance.txt", sep="\t", header=T)

#colnames(micro)
#nrow(micro)
#rownames(micro) <- micro$OTU.ID
#micro <- micro[,grepl("DS", colnames(micro))]
#micro <- log2(micro+1)

	#bring in fly transcripts
fly <- read.table("em_symbols_vst.txt", quote="")
colnames(fly)
#rownames(fly) <- fly$X
#fly <- fly[,2:ncol(fly)]
complete.cases(fly)

fly <- fly[,!colnames(fly) %in% excl]
micro <- micro[,!colnames(micro) %in% excl]

	#ensure you have the same samples in both matrices
all(colnames(fly) %in% colnames(micro))
table(colnames(fly) %in% colnames(micro))
table(colnames(micro) %in% colnames(fly))
inBoth <- intersect(colnames(micro), colnames(fly))
length(inBoth)

fly <- fly[,match(inBoth, colnames(fly))]
micro <- micro[,match(inBoth, colnames(micro))]
all(colnames(fly) %in% colnames(micro))
colnames(micro) == inBoth
colnames(fly) == inBoth
colnames(fly) == colnames(micro)
all(colnames(fly) == colnames(micro))

##################
##	filtering	##
##################

#microBinary <- ifelse(micro==0,0,1)
#microPercFilter <- (rowSums(microBinary)/ncol(microBinary)*100) >= percThresh
#table(microPercFilter)
#barplot(colSums(micro), las=2)
#micro <- micro[microPercFilter,]
#dim(micro)

	#express microbiota as proportion
#micro <- t(apply(micro, 1, function(x){x/colSums(micro)}))
#barplot(colSums(micro), las=2)

#any(fly==0)	#nup.

flyCv <- t(apply(fly, 1, function(x){sd(x)/mean(x)}))
#microCv <- t(apply(micro, 1, function(x){sd(x)/mean(x)}))
plot(density(flyCv, na.rm=T))
#plot(density(microCv, na.rm=T))

fly <- fly[flyCv >= quantile(flyCv, quantThresh, na.rm=T),]
#micro <- micro[microCv >= quantile(microCv, quantThresh, na.rm=T),]

nrow(fly)
nrow(micro)

fly <- t(apply(fly, 1, scale, center=T, scale=T))
micro <- t(apply(micro, 1, scale, center=T, scale=T))

table(complete.cases(fly))
table(complete.cases(micro))

fly <- fly[complete.cases(fly),]
micro <- micro[complete.cases(micro),]

colnames(fly) <- colnames(micro) <- inBoth

##############
##	plot	##
##############

#pheatmap(cor(t(micro)),
#	show_rownames=F,
#	show_colnames=F
#	)

##############
##	MOFA	##
##############

reticulate::use_python("/Users/dob/opt/anaconda3/bin/python")

data <- list(fly=as.matrix(fly), micro=as.matrix(micro))
rownames(data$fly) <- paste(rownames(data$fly), "fly", sep="_")
rownames(data$micro) <- paste(rownames(data$micro), "micro", sep="_")

lapply(data,dim)

MOFAobject <- create_mofa(data)

print(MOFAobject)

plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
#model_opts$num_factors <- 3

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


outfile = file.path(getwd(), "model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk=F)