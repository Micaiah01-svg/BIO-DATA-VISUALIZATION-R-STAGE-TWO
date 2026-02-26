# Normalized counts (HBR vs UHR)
counts_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"
counts <- read.csv(counts_url, header = TRUE, row.names = 1)

# Differential expression results (Chromosome 22)
deg_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
deg <- read.csv(deg_url, header = TRUE)

# Breast Cancer Wisconsin dataset
bc_url <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
bc <- read.csv(bc_url, header = TRUE)

# STAGE TWO TASK
#PART 1: GENE EXPRESSION ANALYSIS
# A
#  Draw the heatmap
library(pheatmap)
pheatmap::pheatmap(
  mat = counts[ , 1:6],
  color = colorRampPalette(c('white','lightblue','navyblue'))(100),
  cluster_rows = T,
  cluster_cols = T,
  fontsize_row = 5,
  fontsize_col = 5,
  labels_row = 'gene_labels',
  main = 'HBR vs UHR',
  legend = T, 
  border_color = NA
)


#B
#Volcano Plot

library(ggplot2)
#Add a significance category
deg$significance <- 'NS'
deg$significance [deg$log2FoldChange >= 1 & deg$PAdj < 0.05] <- 'Up'
deg$significance [deg$log2FoldChange <= -1 & deg$PAdj < 0.05] <- 'Down'
#Volcano plot
ggplot(deg, aes(x = log2FoldChange, y = -log10(PAdj), color = significance)) +
  geom_point() +
  scale_color_manual(values = c('Up' = 'green', 'Down' = 'orange','NS' = 'grey'))+
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed')+
  theme_minimal()

#PART 2- Breast cancer Data Exploration
#C. Scatter Plot (radius vs texture)
library(ggplot2)
ggplot(bc, aes(x = radius_mean, y = texture_mean, color = diagnosis))+
  geom_point() +
  labs(x = 'radius mean', y = 'texture mean', title = 'radius vs texture') +
  theme_minimal()


# D. Correlation Heatmap
# Select the 6 features
features <- bc[, c('radius_mean', 'texture_mean', 'perimeter_mean',
                   'area_mean', 'smoothness_mean', 'compactness_mean')]

# Compute Correlation
cor_matrix <- cor(features)

# Install reshape2 if needed
install.packages("reshape2")

# Load libraries
library(ggplot2)
library(reshape2)

# Heatmap
ggplot(melt(cor_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Correlation Heatmap", x = "", y = "") +
  theme_minimal()

#E. Scatter Plot (Smoothness vs Compactness)
ggplot(bc, aes(x = smoothness_mean, y = compactness_mean, color = diagnosis)) +
  geom_point() +
  labs(x = 'smoothness mean', y = 'compactness mean',
       title = 'smoothness vs compactness') +
  theme_minimal() +
  theme(panel.grid.major = element_line(),
        panel.grid.minor = element_line())

#F. Density Plot(area distribution)
ggplot(bc, aes(x = area_mean, fill = diagnosis)) +
  geom_density(alpha = 0.4) +
  labs(x = 'Area Mean', y = 'Density',
       title = 'Area Mean Distribution by Diagnosis') +
  theme_minimal()





library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(gridExtra)


# helper functions
transparent_color <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", "#59a14f",
            "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63", "#79706e", "#d4a6c8",
            "#e9e9e9", "#ffbe7d", "#bab0ac", "#9d7660", "#d37295", "#86bcb6",
            "#362a39", "#cd9942")

# Load Excel data
file <- "hb_stage_2.xlsx"
df_a  <- read_xlsx(file, sheet="a")
df_b  <- read_xlsx(file, sheet="b")
df_c  <- read_xlsx(file, sheet="c")
df_d1 <- read_xlsx(file, sheet="d_1")
df_e  <- read_xlsx(file, sheet="e")
df_f  <- read_xlsx(file, sheet="f")
df_g  <- read_xlsx(file, sheet="g")

#TASK 1: REPRODUCE PANEL 2A: CELL-TYPE RATIO DISTRIBUTIONS
p2a <- ggplot(df_a, aes(x = cell_type, y = new_ratio, fill = cell_type)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = hb_pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(title = "Panel 2a — Cell-type Ratio Distributions")
print(p2a)
# Conceptual check (Task 1):
# We plot cell-type ratio distributions to understand the baseline immune composition.
# This shows which immune subsets are abundant, variable, or contain outliers before deeper analysis.



# TASK 2: Panel 2b — Half-life vs Alpha with cutoffs, colored regimes, and exemplar genes
library(readxl)
library(ggplot2)
library(dplyr)

# 1️⃣ Load data
df_b <- read_xlsx("hb_stage_2.xlsx", sheet = "b") %>%
  rename(gene = cell) %>%
  mutate(
    log_half  = log2(half_life),
    log_alpha = log2(alpha)
  )

# 2️⃣ Median cutoffs
xcut <- median(df_b$log_half, na.rm = TRUE)
ycut <- median(df_b$log_alpha, na.rm = TRUE)

# 3️⃣ Assign regimes
df_b <- df_b %>%
  mutate(
    regime = case_when(
      log_half >= xcut & log_alpha >= ycut ~ "High/High",
      log_half <  xcut & log_alpha >= ycut ~ "Low Half/High Alpha",
      log_half >= xcut & log_alpha <  ycut ~ "High Half/Low Alpha",
      TRUE                                 ~ "Low/Low"
    )
  )

# 4️⃣ Define colors
regime_colors <- c(
  "High/High"            = "#e15759",  # red
  "Low Half/High Alpha"  = "#8cd17d",  # green
  "High Half/Low Alpha"  = "#4e79a7",  # blue
  "Low/Low"              = "#79706e"   # grey
)

# 5️⃣ Plot
p2b <- ggplot() +
  
  # Base layer: all points light grey
  geom_point(data = df_b, aes(x = log_half, y = log_alpha),
             color = "lightgrey", size = 1.1, alpha = 0.5) +
  
  # Colored subsets for each regime
  geom_point(data = df_b %>% filter(regime == "High/High"),
             aes(x = log_half, y = log_alpha),
             color = regime_colors["High/High"], size = 1.2) +
  geom_point(data = df_b %>% filter(regime == "Low Half/High Alpha"),
             aes(x = log_half, y = log_alpha),
             color = regime_colors["Low Half/High Alpha"], size = 1.2) +
  geom_point(data = df_b %>% filter(regime == "High Half/Low Alpha"),
             aes(x = log_half, y = log_alpha),
             color = regime_colors["High Half/Low Alpha"], size = 1.2) +
  geom_point(data = df_b %>% filter(regime == "Low/Low"),
             aes(x = log_half, y = log_alpha),
             color = regime_colors["Low/Low"], size = 1.2) +
  
  # Median cutoff lines
  geom_vline(xintercept = xcut, linetype = "dashed") +
  geom_hline(yintercept = ycut, linetype = "dashed") +
  
  # Labels for exemplar genes
  geom_text(
    data = df_b %>% filter(gene %in% c("Camp", "Ccr2")),
    aes(x = log_half, y = log_alpha, label = gene),
    color = "black",
    fontface = "bold",
    vjust = -0.5
  ) +
  
  # Theme and labels
  theme_bw() +
  labs(
    x = "log2(Half Life)",
    y = "log2(Alpha)",
    title = "Panel 2b — Kinetic Regimes with Exemplars"
  )

print(p2b)
# Conceptual check (Task 2):
# log2 transforms compress wide numeric ranges and make fold-changes additive.
# The 4 quadrants classify genes by induction rate (alpha) and stability (half-life):
# fast+stable, fast+unstable, slow+stable, slow+unstable.


#Task 3 - Heatmap across cell types and time
library(stringr)
sheet_c <- read_excel('hb_stage_2.xlsx', sheet = 'c')
#inspect dataset
str(sheet_c)
#extract numeric matrix
mat <- as.matrix(sheet_c[ , -1])   # remove genes
rownames(mat) <- sheet_c$genes

#Build column annotations
#seperate the celltypes and time
celltypes <- str_extract(colnames(mat), "^[A-Za-z]+")
celltypes <- str_replace(celltypes, "n$", "")
print(celltypes)

times <- str_extract(colnames(mat), "[0-9]+h")
print(times)

annotation_col <- data.frame(CellType = celltypes,
                             Time = times
)
rownames(annotation_col) <- colnames(mat)

#Plot Heatmap
pheatmap::pheatmap( mat,
                    annotation_col = annotation_col,
                    fontsize = 11,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    scale = "row",
                    color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
                    show_rownames = FALSE,
                    show_colnames = FALSE
)

# Conceptual check (Task 3):
# Genes are clustered to reveal co-regulation patterns.
# Timepoints are NOT clustered because time is sequential and must remain in biological order.



#Task 4 - Pathway enrichment heatmap
sheet_d <- read_excel('hb_stage_2.xlsx', sheet = 'd_1')

#using pathway names as rownames
rownames <- sheet_d$pathway
Mat_d <- as.matrix(sheet_d[ , -1])   # remove pathway column

#centering color scale
range_val <- max(abs(Mat_d))

pheatmap::pheatmap(Mat_d,
                   fontsize = 10,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(c("red", "white", "blue"))(200),
                   breaks = seq(-range_val, range_val, length.out = 200),
                   main = "Pathway Enrichment Across Timepoints",
                   annotation_legend = TRUE,
                   legend = TRUE,
                   labels_row = rownames(Mat_d),
                   annotation_names_row = TRUE
)
# Conceptual check (Task 4):
# No clustering is used because pathways have fixed biological meaning and should not be reordered.
# A diverging palette is required because enrichment values include both up- and down-regulation.



#Task 5 - Bubble plot of kinetic regimes
sheet_e <- read_excel('hb_stage_2.xlsx', sheet = 'e')

#Bubble plot
p2e <- ggplot(sheet_e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.7) +
  
  # Apply project palette for stages
  scale_color_manual(values = c("6h" = "#e15759", "72h" = "#8cd17d")) +
  # Customize size range for readability
  scale_size_continuous(range = c(2, 10)) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Functional Kinetic Regimes",
    x = "Half Life",
    y = "Alpha Life",
    size = "Gene Count",
    color = "Stage"
  )

# Display plot
print(p2e)
# Conceptual check (Task 5):
# Bubble size reflects how many genes share a kinetic signature,
# highlighting which kinetic regimes dominate the immune response.



#Task 6 - Stacked Proportions
sheet_f <- read_excel('hb_stage_2.xlsx', sheet = 'f')

#Subset the data
#Filter for start and end timepoint
sheet_f <- sheet_f %>% 
  filter(stage %in% c("s00h", "s72h"))

#Create Stacked Barplot
p2f <- ggplot(sheet_f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  # Fixed y-axis limit as per requirements
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
  # Apply colors from the project palette (B and Plasma)
  scale_fill_manual(values = c("B" = "#4e79a7", "Plasma" = "#b07aa1")) +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Figure 2f: Stacked Proportion",
    x = "Timepoint",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

#Display plot
print(p2f)

# Conceptual check (Task 6):
# Stacked bars emphasize how two populations divide a single total pool at each timepoint.
# Side-by-side bars would hide the relative composition changes.



#Task 7
sheet_g <- read_excel('hb_stage_2.xlsx', sheet = 'g')

#Convert to Adjacency Matrix
mat_g <- as.matrix(sheet_g[,-1])
rownames(mat_g) <- sheet_g[[1]]
colnames(mat_g) <- sheet_g[[1]]

#Build Directed Graph
net <- graph_from_adjacency_matrix(mat_g, mode = "directed", weighted = TRUE, diag = FALSE)

#Remove Zero-weight Edges
net <- delete.edges(net, which(E(net)$weight <= 0))

#Define Layout and Styling
set.seed(42) 
layout_fr <- layout_with_fr(net)

#Map edge width and arrow size to the weight
E(net)$width <- E(net)$weight * 5
E(net)$arrow.size <- E(net)$weight * 1.5

#Plot the Network
plot(net, 
     layout = layout_fr,
     vertex.color = node_colors,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.label.dist = 2,
     vertex.frame.color = "white",
     edge.color = transparent_color("grey", 30), # Now this function will work!
     edge.curved = 0.2,
     main = "Figure 2g: Cell-Cell Interaction Network")

# Conceptual check (Task 7):
# The network is directed because ligand signals flow from a sender cell to a receiver cell.
# Edge weight encodes communication strength (e.g., ligand–receptor activity level).



# =========================
# Task 8: Final Assembly 
# =========================

# 1. Load Libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)

# 2. File path
file_path <- "hb_stage_2.xlsx"

# 3. Load all sheets
df_a <- read_excel(file_path, sheet = "a")
df_b <- read_excel(file_path, sheet = "b")
sheet_c <- read_excel(file_path, sheet = "c")
sheet_d <- read_excel(file_path, sheet = "d_1")
sheet_e <- read_excel(file_path, sheet = "e")
sheet_f <- read_excel(file_path, sheet = "f") %>% filter(stage %in% c("s00h", "s72h"))

# 4. Clean column names
clean_names <- function(df) {
  colnames(df) <- tolower(gsub(" ", "_", colnames(df)))
  return(df)
}
df_a <- clean_names(df_a)
df_b <- clean_names(df_b)
sheet_c <- clean_names(sheet_c)
sheet_d <- clean_names(sheet_d)
sheet_e <- clean_names(sheet_e)
sheet_f <- clean_names(sheet_f)

# 5. Dynamic color palette for panel 2a
cell_types <- unique(df_a$cell_type)
hb_pal <- setNames(brewer.pal(n = max(3, length(cell_types)), name = "Set3"), cell_types)

# =========================
# 6. Panel 2a: Boxplot
# =========================
p2a <- ggplot(df_a, aes(x = reorder(cell_type, -new_ratio, FUN = median), y = new_ratio, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5) + scale_fill_manual(values = hb_pal) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(title = "a) Cell-type Ratio distribution", x = "", y = "Ratio")

# =========================
# 7. Panel 2b: Scatter
# =========================
p2b <- ggplot(df_b, aes(x = log2(alpha), y = log2(half_life))) +
  geom_point(alpha = 0.2, color = "grey50") + theme_classic() +
  geom_text(
    data = filter(df_b, cell %in% c("Camp", "Ccr2")),
    aes(label = cell),
    color = "red", fontface = "bold"
  ) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(title = "b) Kinetic Landscape (log2 scale)", x = "log2 (Alpha)", y = "log2 (Half-life)")

# =========================
# 8. Panel 2c: Gene Heatmap
# =========================
mat_c <- as.matrix(sheet_c[,-1])
rownames(mat_c) <- sheet_c[[1]]
p2c <- pheatmap(mat_c, cluster_cols = FALSE, scale = "row", silent = TRUE,
                color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
                main = "c) Gene Expression")$gtable

# =========================
# 9. Panel 2d: Pathway Heatmap
# =========================
mat_d <- as.matrix(sheet_d[,-1])
rownames(mat_d) <- sheet_d[[1]]
p2d <- pheatmap(mat_d, cluster_cols = FALSE, cluster_rows = FALSE, silent = TRUE,
                color = colorRampPalette(c("red", "white", "blue"))(200),
                main = "d) Pathway Enrichment Activity")$gtable

# =========================
# 10. Panel 2e: Bubble Plot
# =========================
p2e <- ggplot(sheet_e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.6) + scale_x_log10() + scale_y_log10() + theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size=7),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(title = "e) Functional Kinetic Regimes", x = "Half-life", y = "Alpha")

# =========================
# 11. Panel 2f: Stacked Bar
# =========================
p2f <- ggplot(sheet_f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) + theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(title = "f) Lineage Shift", x = "", y = "Prop.")

#12 Create the layout grid
# We arrange them in a 3-row layout
final_plot <- grid.arrange(
  p2a, p2b, 
  p2c, p2d, 
  p2e, p2f,
  ncol = 2,
  widths = c(1.2, 0.8),
  top = textGrob("Figure 2: Global Immune Response Kinetics", gp = gpar(fontsize=16, fontface="bold"))
)


# =========================
# 13. PRINT TO SCREEN
# =========================
grid.newpage()
grid.draw(final_plot)

# =========================
# 14. EXPORT HIGH-RES PNG
# =========================
ggsave("Figure2_Publication_Ready.png", final_plot, width = 14, height = 17, dpi = 300)


