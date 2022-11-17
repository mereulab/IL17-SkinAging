#' This script performs differential abundance with Sccomp.

# Load packages
library(sccomp)
library(tidyverse)
library(here)

# Load data
lymphoid_adult_vs_aged <- readRDS(here::here("data", "adult_aged", "lymphoid.RDS"))
myeloid_adult_vs_aged <- readRDS(here::here("data", "adult_aged", "myeloid.RDS"))
non_immune_adult_vs_aged <- readRDS(here::here("data", "adult_aged", "non_immune.RDS"))
lymphoid_subclusters_adult_vs_aged <- readRDS(here::here("data", "adult_aged", "lymphoid_subclustered.RDS"))
lymphoid_control_vs_treated <- readRDS(here::here("data", "control_treated", "lymphoid.RDS"))
myeloid_control_vs_treated <- readRDS(here::here("data", "control_treated", "myeloid.RDS"))
non_immune_control_vs_treated <- readRDS(here::here("data", "control_treated", "non_immune.RDS"))

#' Create input lists of Seurat objects corresponding to lymphoid, lymphoid subclustered, myeloid and non-immune 
#' compartments from adult vs. aged dermal cells
seurat_objs_adult_vs_aged <- list(lymphoid_adult_vs_aged,
                                  myeloid_adult_vs_aged,
                                  non_immune_adult_vs_aged,
                                  lymphoid_subclusters_adult_vs_aged)

names(seurat_objs_adult_vs_aged) <- c("lymphoid",
                                      "myeloid",
                                      "non_immune",
                                      "lymphoid_subclusters")

#' Create input lists of Seurat objects corresponding to lymphoid, myeloid and non-immune compartments from IgG control
#' vs. anti-IL-17A/F dermal cells. 
seurat_objs_control_vs_treated <- list(lymphoid_control_vs_treated,
                                       myeloid_control_vs_treated,
                                       non_immune_control_vs_treated)

names(seurat_objs_control_vs_treated) <- c("lymphoid",
                                           "myeloid",
                                           "non_immune")

# Run Sccomp and save Sccomp objects
sccomp <- function(seurat_obj) {
  seurat_obj$cell_group = seurat_obj$annotation
  seurat_obj$sample = seurat_obj$replicate
  seurat_obj$type = seurat_obj$condition
  
  sccomp_obj = seurat_obj |> 
    sccomp::sccomp_glm( formula_composition = ~ type, 
                formula_variability = ~ 1, 
                percent_false_positive = 5, 
                .sample = sample, 
                .cell_group = cell_group)
  
  return(sccomp_obj)
}

sccomp_objs_adult_vs_aged <- seurat_objs_adult_vs_aged %>%
  purrr::imap(~sccomp(.x) %>% 
                saveRDS(., here::here("data", "adult_aged", paste0("sccomp_", .y, ".RDS"))))

sccomp_objs_control_vs_treated <- seurat_objs_control_vs_treated %>%
  purrr::imap(~sccomp(.x) %>% 
                saveRDS(., here::here("data", "control_treated", paste0("sccomp_", .y, ".RDS"))))

