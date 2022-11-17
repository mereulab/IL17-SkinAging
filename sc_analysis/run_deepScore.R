#' This script script runs DeepSCore (A multi-language multi-omics deep learning
#'  model for automatic single-cell label transfer).

# Load packages
library("devtools")
install_github('elimereu/matchSCore2',ref="vignette-workout") 
install_github("rstudio/keras")
library(keras)
library(matchSCore2)
library(Seurat)
library(tidyverse)
library(here)

# Load data
data <- readRDS(here::here("data", "adult_aged", "lymphoid.RDS"))


Idents(data) <-  data$condition
data.list = SplitObject(data)

adult_scaled <- ScaleData(data.list$Adult,features = rownames(data.list$Adult)) 
adult_scaled <- adult_scaled@assays$RNA@scale.data

aged_scaled <- ScaleData(data.list$Aged,features = rownames(data.list$Aged))
aged_scaled <- aged_scaled@assays$RNA@scale.data

features <- rownames(adult_scaled)
enc <- ds_encoder(data = adult_scaled, 
                  genes = features, 
                  dims = 2000,
                  hnodes = c(5000),
                  verbose = T,
                  name="encoder")

adult_feat <- ds_get_features(enc = enc,data = adult_scaled,genes = features)
cluster <- factor(data.list$Adult$cell.type) # cell types

out <- ds_split_data_encoder(features = adult_feat,
                             clus = cluster,
                             prop = 0.8,
                             verbose = T)
mod <- ds_dnn_model(out = out,hnodes = c(100),
                    epochs = 10,
                    batch_size = 32,
                    verbose = T,
                    add_dropout = T,
                    pct_dropout = 0.3)

probs <- mod %>% predict(adult_feat)
head(probs)
probs <- data.frame(probs,check.names = F)
names(probs) <- c("unclassified",out$classes)
rownames(probs) <- rownames(adult_feat)
young_pred <- apply(probs,1,function(x) ifelse(max(x)>0.5,names(probs)[which(x==max(x))],"unclassified"))
t <-table(young_pred,data.list$Adult$cell.type)
t

predy <- vector()
celltype <- data.list$Adult$cell.type
names(celltype) <- rownames(data.list$Adult@meta.data)
for(i in c(1:nrow(probs))){
  cell <- rownames(probs)[i]
  ct <- as.character(celltype[cell])
  p <- probs[i,ct]
  if(i>1){
    predy <- append(x = predy,values = p,after = length(predy))
  }else{
    predy <- p
  }
}
length(predy)
length(celltype)
nrow(probs)

summary(predy)

aged_feat <- ds_get_features(enc = enc,data = aged_scaled,genes = features)

probs <- mod %>% predict(aged_feat)
head(probs)
probs <- data.frame(probs,check.names = F)
names(probs) <- c("unclassified",out$classes)
rownames(probs) <- rownames(aged_feat)
head(probs)
old_pred <- apply(probs,1,function(x) ifelse(max(x)>0.5,names(probs)[which(x==max(x))],"unclassified"))
t <-table(old_pred,data.list$Aged$cell.type)
t




predh <- vector()
celltype <- data.list$Aged$cell.type
names(celltype) <- rownames(data.list$Aged@meta.data)
for(i in c(1:nrow(probs))){
  cell <- rownames(probs)[i]
  ct <- as.character(celltype[cell]) 
  p <- probs[i,ct]
  if(i>1){
    predh <- append(x = predh,values = p,after = length(predh))
  }else{
    predh <- p
  }
}
length(predh)
length(celltype)
nrow(probs)

summary(predh)

pred <- c(predy,predh)
summary(pred)
names(pred) <- c(rownames(adult_feat),rownames(aged_feat))

pred.ord <- pred[colnames(data@assays$RNA@data)]
data$aging.score <- 0.5-pred.ord


# Load data
saveRDS(data, here::here("data", "adult_aged", "lymphoid_deepscore.RDS"))

