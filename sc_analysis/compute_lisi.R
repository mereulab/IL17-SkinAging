#' This script computes the local inverse Simpson’s Index.

# Load packages
library(lisi)
library(tidyverse)
library(here)

# Load data
se_obj1 <- readRDS(here::here("data", "adult_aged", "lymphoid.RDS"))
se_obj2 <- readRDS(here::here("data", "adult_aged", "myeloid.RDS"))
se_obj3 <- readRDS(here::here("data", "adult_aged", "non_immune.RDS"))
se_obj4 <- readRDS(here::here("data", "adult_aged", "lymphoid_batch_effect.RDS"))
se_obj5 <- readRDS(here::here("data", "adult_aged", "myeloid_batch_effect.RDS"))
se_obj6 <- readRDS(here::here("data", "adult_aged", "non_immune_batch_effect.RDS"))
se_obj7 <- readRDS(here::here("data", "control_treated", "lymphoid.RDS"))
se_obj8 <- readRDS(here::here("data", "control_treated", "myeloid.RDS"))
se_obj9 <- readRDS(here::here("data", "control_treated", "non_immune.RDS"))

#' Create input lists of Seurat objects corresponding to lymphoid, lymphoid subclustered, myeloid and non-immune 
#' compartments from adult vs. aged dermal cells

seurat_objs<- list("lymphoid_adult_vs_aged" = se_obj1,
                   "myeloid_adult_vs_aged" = se_obj2 ,
                   "non_immune_adult_vs_aged" = se_obj3,
                   "lymphoid_batch_effect_adult_vs_aged" = se_obj4, 
                   "myeloid_batch_effect_adult_vs_aged" = se_obj5,
                   "non_immune_batch_effect_adult_vs_aged" = se_obj6,
                   "lymphoid_control_vs_treated" = se_obj7,
                   "myeloid_control_vs_treated" = se_obj8,
                   "non_immune_control_vs_treated" = se_obj9)


# Get cell embeddings and metadata (variable of interest)
cell.embeddings <- seurat_objs %>% 
  purrr::imap(~.x@reductions$umap@cell.embeddings)
meta_data <- seurat_objs %>% 
  purrr::imap(~.x@meta.data %>% select('annotation', 'replicate'))

# Compute LISI score
lisi.scores <- purrr::map2(cell.embeddings, meta_data, ~compute_lisi(
  X = .x,
  meta_data = .y,
  label_colnames = "replicate"
) %>%
  cbind(.x))

# save
lisi.scores %>%
  purrr::imap(~saveRDS(.x, here::here("results", "lisi.scores", paste0("lisi_value_", .y, ".RDS"))))

saveRDS(lisi.scores, here::here("results", "lisi.scores", paste0("lisi_scores.RDS")))

