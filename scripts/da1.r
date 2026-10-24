library(tidyverse)
library(tximport)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(apeglm)

samples <- read_csv("data/airway_samples.csv")
samples

samples <- samples |>
  select(run, cell, dex) |>
  mutate(cell = factor(cell),
         dex = factor(dex, levels = c("untrt","trt")))
samples

tx2gene <- read_tsv("data/tx2gene.gencode.v50.tsv.gz")
genes <- read_tsv("data/genes.gencode.v50.tsv.gz")

samples <- samples |>
  select(run, cell, dex) |>
  mutate(cell = factor(cell),
         dex = factor(dex, levels = c("untrt", "trt")))
samples
tx2gene <- read_tsv("data/tx2gene.gencode.v50.tsv.gz")
genes   <- read_tsv("data/genes.gencode.v50.tsv.gz")
own <- c("SRR1039508", "SRR1039509")


files <- if_else(samples$run %in% own,
                 file.path("results/salmon", samples$run, "quant.sf"),
                 file.path("data/quant", samples$run, "quant.sf.gz"))
names(files) <- samples$run
files
all(file.exists(files))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

names(txi)
dim(txi$counts)
round(head(txi$counts), 1)