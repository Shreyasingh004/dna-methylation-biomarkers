if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("minfi", "limma", "GEOquery", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
install.packages("ggplot2")

library(minfi)
library(limma)
library(GEOquery)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

geo_accession <- "GSE53051"
baseDir <- tempdir()
getGEOSuppFiles(geo_accession, baseDir = baseDir, makeDirectory = TRUE)
untar(file.path(baseDir, geo_accession, paste0(geo_accession, "_RAW.tar")), 
      exdir = file.path(baseDir, geo_accession))

idat_dir <- file.path(baseDir, geo_accession)
rgset <- read.metharray.exp(base = idat_dir, extended = TRUE, verbose = TRUE)

mset <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls", reference = 1)
beta_values <- getBeta(mset)

pheno_data <- data.frame(
  sample = colnames(beta_values),
  group = c(rep("healthy", 10), rep("diseased", 10))
)
rownames(pheno_data) <- pheno_data$sample

design <- model.matrix(~group, data = pheno_data)
fit <- lmFit(beta_values, design)
fit <- eBayes(fit)
diff_methyl <- topTable(fit, coef = 2, number = Inf)

heatmap_data <- beta_values[rownames(diff_methyl)[1:100], ]
heatmap(heatmap_data, ColSideColors = as.character(pheno_data$group))

volcanoplot(fit, coef = 2, highlight = 10, names = rownames(diff_methyl))

significant_cpgs <- subset(diff_methyl, adj.P.Val < 0.05)
write.csv(significant_cpgs, "significant_cpgs.csv", row.names = TRUE)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotated_cpgs <- merge(significant_cpgs, ann450k, by.x = "row.names", by.y = "Name")
write.csv(annotated_cpgs, "annotated_cpgs.csv", row.names = FALSE)
