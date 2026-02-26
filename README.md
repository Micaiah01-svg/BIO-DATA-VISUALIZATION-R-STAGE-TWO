# BIO-DATA-VISUALIZATION-R-STAGE-TWO
Gene Expression Analysis, Breast Cancer Exploration, and Immune Kinetics 
BY: MICAIAH ADEOLUWA ADEDEJI
This project brings together several small but important analytical workflows commonly used in bioinformatics and computational biology. It includes:
•	Working with normalized RNA-seq counts
•	Performing differential expression visualization
•	Exploring diagnostic patterns in breast cancer data
•	Reconstructing multi-panel immune kinetic figures (HackBio Stage 2)
Every part of the analysis is reproducible, and each section builds on core ideas used in modern genomics data science.
________________________________________
 Setup and Required Packages
The analysis relies on the following R packages:
•	ggplot2 — general plotting
•	pheatmap — heatmap visualization
•	readxl — reading Excel sheets
•	igraph — network reconstruction
•	reshape2 — matrix melting for heatmaps
•	gridExtra — arranging multi-panel figures
•	dplyr — data manipulation
•	RColorBrewer — color palettes
Install if necessary:
install.packages(c("ggplot2", "pheatmap", "readxl", "igraph", "reshape2",
                   "gridExtra", "dplyr", "RColorBrewer"))
________________________________________ PART 1 — Gene Expression Analysis
This section focuses on gene-level visualization between two samples: HBR and UHR.
Two plots are produced:
1.	A heatmap of normalized counts
2.	A volcano plot of differentially expressed genes
________________________________________
 1A. Heatmap of Normalized RNA-seq Counts
We begin by loading the dataset:
counts_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"
counts <- read.csv(counts_url, header = TRUE, row.names = 1)
We then visualize the first six samples using pheatmap, applying a blue gradient to show expression intensity:
pheatmap(
  mat = counts[, 1:6],
  color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 5,
  fontsize_col = 5,
  main = "HBR vs UHR"
)
 Conceptual Note (Heatmap)
•	Heatmaps help reveal co-expression structure across genes and samples.
•	Clustering rows and columns groups together patterns that behave similarly.
•	A smooth gradient helps quickly distinguish high vs. low expression.
________________________________________
 1B. Volcano Plot of Differential Expression (Chromosome 22)
Load DEG dataset:
deg_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
deg <- read.csv(deg_url, header = TRUE)
Categorize genes based on fold-change and adjusted p-value:
deg$significance <- "NS"
deg$significance[deg$log2FoldChange >= 1 & deg$PAdj < 0.05] <- "Up"
deg$significance[deg$log2FoldChange <= -1 & deg$PAdj < 0.05] <- "Down"
Plot volcano:
ggplot(deg, aes(x = log2FoldChange, y = -log10(PAdj), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "green", "Down" = "orange", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal()
 Conceptual Note (Volcano Plot)
•	log2 fold-change gives symmetry between up- and down-regulation.
•	–log10(Padj) emphasizes low p-values (strong significance).
•	Volcano plots show a clean contrast between effect size and statistical confidence.
________________________________________

 PART 2 — Breast Cancer Diagnostic Data Exploration
This dataset comes from the Wisconsin Breast Cancer Diagnostic dataset.
We explore its structure visually through scatter plots, correlation heatmaps, and density curves.
bc_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
bc <- read.csv(bc_url)
________________________________________
 2C. Scatter Plot — Radius vs Texture
ggplot(bc, aes(radius_mean, texture_mean, color = diagnosis)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Radius vs Texture")
 Conceptual Note
These two features are among the most diagnostic.
Malignant samples tend to form distinct clusters.
________________________________________
 2D. Correlation Heatmap of Six Key Features
Compute correlation matrix:
features <- bc[, c("radius_mean", "texture_mean", "perimeter_mean",
                   "area_mean", "smoothness_mean", "compactness_mean")]
cor_matrix <- cor(features)
Plot heatmap with values:
ggplot(melt(cor_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal()
 Conceptual Note
•	Radiological features such as radius, perimeter, and area are strongly correlated.
•	Smoothness and compactness capture texture characteristics rather than shape.
________________________________________
 2E. Scatter Plot — Smoothness vs Compactness
ggplot(bc, aes(smoothness_mean, compactness_mean, color = diagnosis)) +
  geom_point() +
  theme_minimal()
 Conceptual Note
This reveals structure between two texture-based features.
Malignant tumors often show higher compactness.
________________________________________
 2F. Density Plot — Area Mean
ggplot(bc, aes(area_mean, fill = diagnosis)) +
  geom_density(alpha = 0.4) +
  theme_minimal()
 Conceptual Note
Density curves show differences in distribution rather than individual values.
Malignant tumors generally have larger area measurements.

 PART 3 — HackBio Stage 2 Tasks 
Using the file hb_stage_2.xlsx, we reconstruct seven figure panels and a final multi-panel assembly.
________________________________________
 Task 1 — Panel 2a: Cell-Type Ratio Distributions
p2a <- ggplot(df_a, aes(cell_type, new_ratio, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = hb_pal)
 Conceptual Check
•	Boxplots summarize relative abundance of immune cells.
•	Outliers show rare or burst-like populations.
•	Rotated labels improve readability when cell types are many.
________________________________________
 Task 2 — Panel 2b: Half-Life vs Alpha (Kinetic Regimes)
Data is transformed into log2-space and split into quadrants using medians:
df_b$log_half <- log2(df_b$half_life)
df_b$log_alpha <- log2(df_b$alpha)
 Conceptual Checks
•	Why log2?
Half-life and induction rates vary over orders of magnitude.
Log2 compresses the range and centers fold-changes.
•	What do the quadrants represent?
o	High alpha + high half-life → strong, stable responders
o	High alpha + low half-life → fast but unstable
o	Low alpha + high half-life → slow but persistent
o	Low alpha + low half-life → weak responders
These kinetic regimes describe temporal immune behavior.
________________________________________
 Task 3 — Panel 2c: Heatmap Across Cell Types & Time
This heatmap features:
•	Clustered genes
•	Unclustered timepoints
•	Column annotations (cell type & time)
 Conceptual Check
Time is a fixed biological sequence, so it should not be rearranged.
Clustering genes reveals shared transcriptional programs.
________________________________________
 Task 4 — Panel 2d: Pathway Enrichment Heatmap
 Conceptual Checks
•	No clustering: pathway order is biologically curated.
•	Diverging palette: values span negative (down) to positive (up) enrichment, requiring a centered color scale.
________________________________________
 Task 5 — Bubble Plot of Kinetic Regimes
Bubble area reflects gene count:
geom_point(aes(size = count, color = stage))
 Conceptual Check
The plot communicates which kinetic regimes are populated heavily or sparsely—an overview of system-level behavior.
________________________________________
 Task 6 — Stacked Barplot (B vs Plasma)
geom_bar(position = "stack")
 Conceptual Check
Stacked bars emphasize the relative fractions of two cell types, not their absolute numbers.
Side-by-side bars wouldn’t show how they divide a shared total.
________________________________________
Task 7 — Directed Cell–Cell Interaction Network
 Conceptual Checks
•	The graph is directed because ligand–receptor communication flows from sender → receiver.
•	Edge weight represents interaction strength, such as signaling magnitude.
________________________________________
Task 8 — Final Multi-Panel Assembly
Panels are arranged using grid.arrange and exported as a publication-ready PNG.
________________________________________
Final Notes
This project brings together multiple skills in exploratory analysis, graphical representation, kinetic modeling, and biological interpretation. Every plot has been constructed with clarity and biological relevance in mind.
If you intend to publish the outputs (e.g., on LinkedIn or GitHub), the multi-panel figure exported in Task 8 makes a great visual summary.

