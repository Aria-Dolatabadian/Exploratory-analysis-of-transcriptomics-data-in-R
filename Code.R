#https://tavareshugo.github.io/data-carpentry-rnaseq/

library(tidyverse)


raw_cts <- read_csv("counts_raw.csv")
trans_cts <- read_csv("counts_transformed.csv")
sample_info <- read_csv("sample_info.csv")
test_result <- read_csv("test_result.csv")


# "gather" the counts data
trans_cts_long <- trans_cts %>% 
  pivot_longer(cols = wt_0_r1:mut_180_r3, 
               names_to = "sample", 
               values_to = "cts")

trans_cts_long


trans_cts_long %>% 
  pivot_wider(names_from = "sample", values_from = "cts")


sample_info

trans_cts_long <- full_join(trans_cts_long, sample_info, by = "sample")

trans_cts_long

trans_cts_long %>%
  ggplot(aes(cts, colour = replicate)) + 
  geom_freqpoly(binwidth = 1) + 
  facet_grid(strain ~ minute)


trans_cts %>% 
  ggplot(aes(wt_0_r1, wt_0_r2)) + geom_point() +
  geom_abline(colour = "brown")


# Calculate all correlations 
trans_cts_corr <- trans_cts %>% 
  # remove the column "gene", which we do not want to calculate correlation on
  select(-gene) %>% 
  # we use Spearman's correlation, a non-parametric metric based on ranks
  cor(method = "spearman")

# Visualise the correlations between the first 5 samples
trans_cts_corr[1:5, 1:5]


library(corrr)

rplot(trans_cts_corr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




raw_cts_long <- raw_cts %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "cts") %>% 
  full_join(sample_info, by = "sample")



summary(raw_cts_long$cts)


summary(trans_cts_long$cts)


raw_cts %>% 
  ggplot(aes(wt_0_r1, wt_0_r2)) + 
  geom_point()



# We add a "pseudocount" of 1 to the count data because otherwise log(0) = -Inf
raw_cts %>% 
  ggplot(aes(wt_0_r1 + 1, wt_0_r2 + 1)) + 
  geom_point() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2")



# Summarise the data and then plot
raw_cts_long %>% 
  # for each gene
  group_by(gene) %>% 
  # get mean and variance
  summarise(mean_cts = mean(cts),
            var_cts = var(cts)) %>% 
  # plot one against the other
  ggplot(aes(mean_cts, var_cts)) +
  geom_point() +
  geom_abline(colour = "brown") +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2")



raw_cts_long <- raw_cts %>% 
  pivot_longer(wt_0_r1:mut_180_r3, names_to = "sample", values_to = "cts")


# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))

# Make the plot
raw_cts_long %>%
  # add pseudo-count of 1 because log(0) = -Inf
  ggplot(aes(log10(cts + 1), colour = replicate)) + 
  geom_freqpoly(binwidth = 1) +
  facet_grid(rows = vars(strain), cols = vars(minute))


# Make a boxplot
raw_cts_long %>%
  # make sure minute is specified as a factor
  ggplot(aes(factor(minute), log10(cts + 1), fill = strain)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(replicate))


# Scatterplot between T0 and T30
# the correlation is lower than between replicates at T0, for example
trans_cts %>% 
  ggplot(aes(wt_0_r1, wt_30_r1)) + geom_point() +
  geom_abline(colour = "brown")


# 1. class of the object
class(sample_pca)

# 2. structure of the object
str(sample_pca)

# 3. checking the help ?prcomp, under the section "Value" is says:
# "sdev" contains the standard deviation explained by each PC, so if we square it we get the eigenvalues (or explained variance)
# "rotation" contains the variable loadings for each PC, which define the eigenvectors
# "x" contains the PC scores, i.e. the data projected on the new PC axis
# "center" in this case contains the mean of each gene, which was subtracted from each value
# "scale" contains the value FALSE because we did not scale the data by the standard deviation

# 4. we can use the 'dollar sign' to access these elements
pc_scores <- sample_pca$x              # PC scores (a matrix)
pc_eigenvalues <- sample_pca$sdev^2    # eigenvalues (a vector) - notice we square the values
pc_loadings <- sample_pca$rotation     # variable loadings (a matrix)

# 5. here's three ways to check this
ncol(pc_scores)
length(pc_eigenvalues)
ncol(pc_loadings)


pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(minute), shape = strain)) +
  geom_point()

# print the result (in this case a ggplot)
pca_plot


loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot


# 1. making the MA plot
test_result %>% 
  ggplot(aes(log10(baseMean), log2FoldChange)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(comparison))

test_result %>% 
  # add column which contains value only if padj < 0.01
  mutate(sig = ifelse(padj < 0.01, log2FoldChange, NA)) %>% 
  # make the plot
  ggplot(aes(baseMean, log2FoldChange)) +
  geom_point(alpha = 0.1) +
  geom_point(aes(y = sig), colour = "brown", size = 1) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(vars(comparison))



##### setup ####

# load packages
library(tidyverse)

# read the data
trans_cts <- read_csv("counts_transformed.csv")
sample_info <- read_csv("sample_info.csv")



# Create a matrix from our table of counts
pca_matrix <- trans_cts %>% 
  # make the "gene" column become the rownames of the table
  column_to_rownames("gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

# Perform the PCA
sample_pca <- prcomp(pca_matrix)


# Look at the first 10 rows and first 5 columns of the matrix
pca_matrix[1:10, 1:5]

# Convert matrix to tibble
as_tibble(pca_matrix)


# Convert matrix to tibble - add colnames to a new column called "gene"
as_tibble(pca_matrix, rownames = "sample")


pc_eigenvalues <- sample_pca$sdev^2


# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues


pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x


pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores


pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

pc_loadings <- sample_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings

top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes


top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)


library(patchwork)

# Adjust some aspects of each plot
pca_plot <- pca_plot + 
  coord_fixed(ratio = 0.4) + 
  labs(title = "PC scores") +
  theme(legend.position = "none")

loadings_plot <- loadings_plot + 
  coord_fixed(ratio = 0.4) + 
  labs(title = "PC loadings")

# Put them together
(pca_plot | loadings_plot) + plot_annotation(tag_levels = "A")

library(ggfortify)
autoplot(sample_pca)

autoplot(sample_pca, data = sample_info, colour = "minute", shape = "strain")

# Example using iris dataset
autoplot(prcomp(iris[, -5]), data = iris, colour = "Species",
         loadings = TRUE, loadings.label = TRUE)

library(broom)

# PC variances (eigen values)
tidy(sample_pca, matrix = "eigenvalues")

# variable loadings
tidy(sample_pca, matrix = "loadings")




##### setup ####

# load packages
library(tidyverse)

# read the data
trans_cts <- read_csv("counts_transformed.csv")
sample_info <- read_csv("sample_info.csv")
test_result <- read_csv("test_result.csv")


candidate_genes <- test_result %>% 
  filter(padj < 0.01) %>%    # filter table
  pull(gene) %>%             # extract the gene column as a vector
  unique()                   # retain only unique values


trans_cts_long <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = wt_0_r1:mut_180_r3, names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample"))


trans_cts_mean <- trans_cts_long %>% 
  # filter genes of interest
  filter(gene %in% candidate_genes) %>% 
  # for each gene, strain and minute
  group_by(gene, strain, minute) %>% 
  # calculate mean and number of replicates
  summarise(mean_cts = mean(cts),
            nrep = n()) %>% 
  # remove grouping from downstream analysis
  ungroup()

head(trans_cts_mean)


trans_cts_mean %>% 
  ggplot(aes(minute, mean_cts)) +
  geom_line(aes(group = gene), alpha = 0.3) +
  facet_grid(rows = vars(strain))


trans_cts_mean <- trans_cts_long %>% 
  # filter to retain only genes of interest
  filter(gene %in% candidate_genes) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene, strain, minute) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()


trans_cts_mean %>%
  ggplot(aes(minute, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.2) + 
  geom_hline(yintercept = 0, colour = "brown", linetype = "dashed") +
  facet_grid(rows = vars(strain))




##### setup ####

# load packages
library(tidyverse)

# read the data
trans_cts <- read_csv("counts_transformed.csv")
sample_info <- read_csv("sample_info.csv")
test_result <- read_csv("test_result.csv")


##### get counts for candidate genes ####

# set of candidate genes for clustering
candidate_genes <- test_result %>% 
  filter(padj < 0.01) %>%    # filter table
  pull(gene) %>%             # extract the gene column as a vector
  unique()                   # retain only unique values

# Summarise counts 
trans_cts_mean <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = wt_0_r1:mut_180_r3, names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # filter to retain only genes of interest
  filter(gene %in% candidate_genes) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene, strain, minute) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()



# Create a matrix
hclust_matrix <- trans_cts %>% 
  select(-gene) %>% 
  as.matrix()

# assign rownames
rownames(hclust_matrix) <- trans_cts$gene


hclust_matrix <- hclust_matrix[candidate_genes, ]

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


gene_dist <- dist(hclust_matrix)

gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram


cutree(gene_hclust, k = 5)


gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

head(gene_cluster)


trans_cts_cluster <- trans_cts_mean %>% 
  inner_join(gene_cluster, by = "gene")

head(trans_cts_cluster)


trans_cts_cluster %>% 
  ggplot(aes(minute, mean_cts_scaled)) +
  geom_line(aes(group = gene)) +
  facet_grid(rows = vars(strain), cols = vars(cluster))

trans_cts_cluster %>% 
  ggplot(aes(minute, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(strain), cols = vars(cluster))



#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)
Heatmap(hclust_matrix, show_row_names = FALSE)





