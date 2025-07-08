
# 🧬 DNA Methylation Biomarker Discovery

**Goal:**  
Identify DNA methylation markers that differentiate between healthy and diseased samples.

---

## 🧠 Skills Used

- R (Bioconductor)
- DNA methylation analysis
- Differential methylation statistics
- Data visualization (heatmaps, volcano plots)
- Genomic annotation

---

## 📁 Dataset

- GEO Accession: [GSE53051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53051)

---

## 🧪 Methodology

1. **Download** methylation beta values from GEO.
2. **Preprocess** IDAT files using Illumina pipeline.
3. **Perform** differential methylation analysis using `limma`.
4. **Visualize** top CpG sites using heatmaps and volcano plots.
5. **Extract** significant CpG markers (adjusted p-value < 0.05).
6. **Annotate** CpG sites with gene information using Illumina 450k annotation.

---

## 📂 Files

| File                      | Description                                           |
|---------------------------|-------------------------------------------------------|
| `methylation_analysis.R`  | Full R script for the analysis pipeline               |
| `significant_cpgs.csv`    | List of significantly differentially methylated CpGs |
| `annotated_cpgs.csv`      | Annotated CpGs with gene information                 |

---

## 📷 Output Visualizations

- **Heatmap:** Top 100 CpGs across all samples  
- **Volcano Plot:** p-value vs logFC for all probes
### 🧪 Output Plots

#### 🔥 Volcano Plot
![Volcano Plot](volcano.png)

#### 🧊 Heatmap of Top 100 CpGs
![Heatmap](heatmap.png)

---

## ✅ Requirements

Install these R packages before running:

```r
install.packages("BiocManager")
BiocManager::install(c("minfi", "limma", "GEOquery", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
install.packages("ggplot2")
