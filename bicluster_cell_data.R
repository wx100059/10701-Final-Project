"This program compares the performance between the k means, spectral clustering,
convex biclustering and the iterative spectral biclustering in different decomposition level."

# Initialization
setwd("/Users/xiangwang/Library/CloudStorage/Box-Box/courses/10701 Machine Learning/Homework/Final Project/simulation")

library(Matrix)
library(igraph)
library(RMThreshold)
library(lpSolve)
library(mcclust)
library(ggplot2)
library(sp)
library(maps)
library(foreign)
library(shapefiles)
library(CVXR)
library(anocva)
library(fpc)
library(kernlab)
library(powerplus)
library(stats)
library(cccd)
library(pracma)
library(usethis)
library(devtools)
library(roxygen2)
library(proxy)
library(devtools)
library(hdf5r)
document(pkg='./mulbicluster')
library(mulbicluster)
document(pkg='cvxbiclustr')
library(cvxbiclustr)
start_time <- Sys.time()

## control the number of iterations
iterate <- 1;  
start_time<-Sys.time()

## read lung cancer data
lung_cancer<-read.table("D:\\OneDrive\\Research\\Master thesis\\simulation\\data.txt")
lung_cancer<-data.matrix(lung_cancer)
storage.mode(lung_cancer)<-"double"

## choose 500 genes with biggest variances and show the heatmap after reorder
num_of_gene <- 200 # how many genes with biggest variation are selected
var_lung_cancer<-apply(lung_cancer,1,var)
threshold<-sort(var_lung_cancer,decreasing = TRUE)[1:num_of_gene]
index<-match(threshold,var_lung_cancer)
X<-lung_cancer[index,]

# # set different level of noise
n <- 50
noise_low <- 0
noise_high <- 5.0
noise_seq <- seq(from = noise_low, to = noise_high, by = (noise_high - noise_low)/(n-1))

## store clustering accuracy of different algorithm under different noise level
## RI, ARI and VI are three metrics for clustering accuracy
## multiple is the multiple sample based clustering algorithm
## multiple_mean is the multiple sample based clustering algorithm only use the mean vector
## svd_cluster is the svd separate clustering algorithm
## cvx_cluster is the convex clustering algorithm
## kmeans is the k means algorithm
## spectral is the spectral clustering algorithm

# VI_kmeans_multiple <- rep(x = 0, times = n)
# VI_kmeans_multiple_mean <- rep(x = 0, times = n)
VI_spectral_multiple <- rep(x = 0, times = n)
# VI_svd_cluster <- rep(x = 0, times = n)
# VI_kmeans <- rep(x = 0, times = n)
VI_spectral <- rep(x = 0, times = n)
VI_cvx_bicluster <- rep(x = 0, times = n)
VI_cvx_cluster <- rep(x = 0, times = n)

# RI_kmeans_multiple <- rep(x = 0, times = n)
# RI_kmeans_multiple_mean <- rep(x = 0, times = n)
RI_spectral_multiple <- rep(x = 0, times = n)
# RI_svd_cluster <- rep(x = 0, times = n)
# RI_kmeans <- rep(x = 0, times = n)
RI_spectral <- rep(x = 0, times = n)
RI_cvx_bicluster <- rep(x = 0, times = n)
RI_cvx_cluster <- rep(x = 0, times = n)

# ARI_kmeans_multiple <- rep(x = 0, times = n)
# ARI_kmeans_multiple_mean <- rep(x = 0, times = n)
ARI_spectral_multiple <- rep(x = 0, times = n)
# ARI_svd_cluster <- rep(x = 0, times = n)
# ARI_kmeans <- rep(x = 0, times = n)
ARI_spectral <- rep(x = 0, times = n)
ARI_cvx_bicluster <- rep(x = 0, times = n)
ARI_cvx_cluster <- rep(x = 0, times = n)

# rank of the approximated matrix
r <- 5
X_noised <- list()

## compute the baseline of each clustering algorithm
## column_cluster is given by the previous knowledge
## row_cluster is given by applying each clustering algortihm on the gene data
column_cluster <- c(rep(x = 1, times = 20), rep(x = 2, times = 13), rep(x = 3, times = 17), rep(x = 4, times = 6))
row_cluster <- c(rep(x = 1, times = as.integer(num_of_gene/2)), rep(x = 2, times = num_of_gene - as.integer(num_of_gene/2)))
k_col <- length(unique(column_cluster))
k_row <- length(unique(row_cluster))

X_baseline <- matrix(data = 0, nrow= size(X)[1], ncol = size(X)[2])
for(i in 1:num_of_gene)
  for(j in 1:ncol(X))
  {
    X_baseline[i,j] = (row_cluster[i] - 1) * k_col + column_cluster[j]
  }

## Set baselines for each biclustering aglorithm
# set the initial store space
X_baseline_spec_mul <- matrix(data = 0, nrow = size(X)[1], ncol = size(X)[2])
X_baseline_spec <- matrix(data = 0, nrow = size(X)[1], ncol = size(X)[2])
X_baseline_cvxbi <- matrix(data = 0, nrow = size(X)[1], ncol = size(X)[2])
X_baseline_cvx <- matrix(data = 0, nrow = size(X)[1], ncol = size(X)[2])

# Compute the initial clustering results with four algorithms
spectral_multiple_accuracy <- spectral_multiple_bicluster(X = X, X_baseline = X_baseline, type = "MEAN")
# svd_accuracy <- svd_separated_bicluster(X = X_noised, X_baseline = X_baseline, r = 3)
# kmeans_accuracy <- kmeans_separate_bicluster(X = X_noised[[it]], X_baseline = X_baseline)
spectral_accuracy <- spectral_separate_bicluster(X = X, X_baseline = X_baseline)
cvx_bicluster_accuracy <- cvx_bicluster(X = X, X_baseline = X_baseline, phi = 5)
cvx_cluster_accuracy <- cvx_separate_cluster(X = X, X_baseline = X_baseline, col_k = 4, row_gamma = 50, col_gamma = 50)

k_col_spec_mul <- length(unique(spectral_multiple_accuracy$col_cluster))
k_row_spec_mul <- length(unique(spectral_multiple_accuracy$row_cluster))
k_col_spec <- length(unique(spectral_accuracy$col_cluster))
k_row_spec <- length(unique(spectral_accuracy$row_cluster))
k_col_cvxbi <- length(unique(cvx_bicluster_accuracy$clu_col))
k_row_cvxbi <- length(unique(cvx_bicluster_accuracy$clu_row))
k_col_cvx <- length(unique(cvx_bicluster_accuracy$clu_row))
k_row_cvx <- length(unique(cvx_bicluster_accuracy$clu_col))


for(i in 1:num_of_gene)
  for(j in 1:ncol(X))
  {
    X_baseline_spec_mul[i,j] = (spectral_multiple_accuracy$row_cluster[i] - 1) * k_col_spec_mul + spectral_multiple_accuracy$col_cluster[j]
    X_baseline_spec[i,j] = (spectral_accuracy$row_cluster[i] - 1) * k_col_spec + spectral_accuracy$col_cluster[j]
    X_baseline_cvxbi[i,j] = (cvx_bicluster_accuracy$clu_row[i] - 1) * k_col_cvxbi + cvx_bicluster_accuracy$clu_col[j]
    X_baseline_cvx[i,j] = (cvx_cluster_accuracy$clu_row[i] - 1) * k_col_cvx + cvx_cluster_accuracy$clu_col[j]
  }

for (i in 1:length(noise_seq))
{
  ## Add noise to the ground truth matrix and build low rank approximation matrix.
  for (it in 1:iterate)
  {
    X_noised[[it]] <- RMThreshold::add.Gaussian.noise(mat = X, mean = 0, stddev = noise_seq[i], symm = FALSE)
    
    ## SVD pre-processing and get the low-rank approximate
    # s <- svd(X_noised[[it]])
    # Sigma_matrix <- diag(x = s$d[1:r], nrow = r, ncol = r)
    # X_noised[[it]] <- s$u[ ,1:r] %*% Sigma_matrix %*% t(s$v[ ,1:r])
  }
  
  for (it in 1:iterate)
  {
    
    ## Use three algorithms to find row clusters
    # spectral_multiple_accuracy <- spectral_multiple_bicluster(X = X_noised, X_baseline = X_baseline)
    # kmeans_multiple_accuracy <- kmeans_multiple_bicluster(X = X_noised[[it]], X_baseline = X_baseline, type = 'KL')
    # kmeans_multiple_accuracy <- kmeans_multiple_bicluster(X = X_noised[[it]], X_baseline = X_baseline, type = 'mean')
    spectral_multiple_accuracy <- spectral_multiple_bicluster(X = X_noised[[it]], X_baseline = X_baseline_spec_mul, type = "MEAN")
    # svd_accuracy <- svd_separated_bicluster(X = X_noised, X_baseline = X_baseline, r = 3)
    # kmeans_accuracy <- kmeans_separate_bicluster(X = X_noised[[it]], X_baseline = X_baseline)
    spectral_accuracy <- spectral_separate_bicluster(X = X_noised[[it]], X_baseline = X_baseline_spec)
    cvx_bicluster_accuracy <- cvx_bicluster(X = X_noised[[it]], X_baseline = X_baseline_cvxbi)
    cvx_cluster_accuracy <- cvx_separate_cluster(X = X_noised[[it]], X_baseline = X_baseline_cvx)
    
    ## Store the clustering accuracy
    # VI_kmeans_multiple[i] <- VI_kmeans_multiple[i] + kmeans_multiple_accuracy$VI
    # VI_kmeans_multiple_mean[i] <- VI_kmeans_multiple_mean[i] + kmeans_multiple_accuracy$VI
    VI_spectral_multiple[i] <- VI_spectral_multiple[i] + spectral_multiple_accuracy$VI
    # VI_svd_cluster[i] <- VI_svd_cluster[i] + svd_accuracy$VI
    # VI_kmeans[i] <- VI_kmeans[i] + kmeans_accuracy$VI
    VI_spectral[i] <- VI_spectral[i] + spectral_accuracy$VI
    VI_cvx_bicluster[i] <- VI_cvx_bicluster[i] + cvx_bicluster_accuracy$VI
    VI_cvx_cluster[i] <- VI_cvx_cluster[i] + cvx_cluster_accuracy$VI
    
    # RI_kmeans_multiple[i] <- RI_kmeans_multiple[i] + kmeans_multiple_accuracy$RI
    # RI_kmeans_multiple_mean[i] <- RI_kmeans_multiple_mean[i] + kmeans_multiple_accuracy$RI
    RI_spectral_multiple[i] <- RI_spectral_multiple[i] + spectral_multiple_accuracy$RI
    # RI_svd_cluster[i] <- RI_svd_cluster[i] + svd_accuracy$RI
    # RI_kmeans[i] <- RI_kmeans[i] + kmeans_accuracy$RI
    RI_spectral[i] <- RI_spectral[i] + spectral_accuracy$RI
    RI_cvx_bicluster[i] <- RI_cvx_bicluster[i] + cvx_bicluster_accuracy$RI
    RI_cvx_cluster[i] <- RI_cvx_cluster[i] + cvx_cluster_accuracy$RI
    # 
    # ARI_kmeans_multiple[i] <- ARI_kmeans_multiple[i] + kmeans_multiple_accuracy$ARI
    # ARI_kmeans_multiple_mean[i] <- ARI_kmeans_multiple_mean[i] + kmeans_multiple_accuracy$ARI
    ARI_spectral_multiple[i] <- ARI_spectral_multiple[i] + spectral_multiple_accuracy$ARI
    # ARI_svd_cluster[i] <- ARI_svd_cluster[i] + svd_accuracy$ARI
    # ARI_kmeans[i] <- ARI_kmeans[i] + kmeans_accuracy$ARI
    ARI_spectral[i] <- ARI_spectral[i] + spectral_accuracy$ARI
    ARI_cvx_bicluster[i] <- ARI_cvx_bicluster[i] + cvx_bicluster_accuracy$ARI
    ARI_cvx_cluster[i] <- ARI_cvx_cluster[i] + cvx_cluster_accuracy$ARI
    
    # accuracy <- c(multiple_accuracy1$RI, multiple_accuracy1$ARI,multiple_accuracy1$VI,
    #               cvx_bicluster_accuracy1$RI, cvx_bicluster_accuracy1$ARI, cvx_bicluster_accuracy1$VI,
    #               svd_accuracy1$RI, svd_accuracy1$ARI, svd_accuracy1$VI)
    # accuracy_sum1 <- accuracy_sum1 + accuracy
  }
}


# VI_kmeans_multiple <- VI_kmeans_multiple / iterate
# VI_kmeans_multiple_mean <- VI_kmeans_multiple_mean / iterate
VI_spectral_multiple <- VI_spectral_multiple / iterate
# VI_svd_cluster <- VI_svd_cluster / iterate
# VI_kmeans <- VI_kmeans / iterate
VI_spectral <- VI_spectral / iterate
VI_cvx_bicluster <- VI_cvx_bicluster / iterate
VI_cvx_cluster <- VI_cvx_cluster / iterate

# RI_kmeans_multiple <- RI_kmeans_multiple / iterate
# RI_kmeans_multiple_mean <- RI_kmeans_multiple_mean / iterate
# RI_spectral_multiple <- RI_spectral_multiple / iterate
RI_spectral_multiple <- RI_spectral_multiple / iterate
# RI_svd_cluster <- RI_svd_cluster / iterate
# RI_kmeans <- RI_kmeans / iterate
RI_spectral <- RI_spectral / iterate
RI_cvx_bicluster <- RI_cvx_bicluster / iterate
RI_cvx_cluster <- RI_cvx_cluster / iterate

# ARI_kmeans_multiple <- ARI_kmeans_multiple / iterate
# ARI_kmeans_multiple_mean <- ARI_kmeans_multiple_mean / iterate
# ARI_spectral_multiple <- ARI_spectral_multiple / iterate
ARI_spectral_multiple <- ARI_spectral_multiple / iterate
# ARI_svd_cluster <- ARI_svd_cluster /iterate
# ARI_kmeans <- ARI_kmeans / iterate
ARI_spectral <- ARI_spectral / iterate
ARI_cvx_bicluster <- ARI_cvx_bicluster / iterate
ARI_cvx_cluster <- ARI_cvx_cluster / iterate

algorithm <- rep(x = c("multiple spectral mean", "spectral clustering", "cvx bicluster", "cvx cluster"), each = n)
VI_accuracy <- data.frame("index" = c(VI_spectral_multiple, VI_spectral,VI_cvx_bicluster, VI_cvx_cluster),
                          "noise" = rep(x = noise_seq, times = 4),
                          "algorithm" = algorithm)
RI_accuracy <- data.frame("index" = c(RI_spectral_multiple, RI_spectral, RI_cvx_bicluster, RI_cvx_cluster),
                          "noise" = rep(x = noise_seq, times = 4),
                          "algorithm" = algorithm)
ARI_accuracy <- data.frame("index" = c(ARI_spectral_multiple, ARI_spectral, ARI_cvx_bicluster, ARI_cvx_cluster),
                           "noise" = rep(x = noise_seq, times = 4),
                           "algorithm" = algorithm)

algorithm <- rep(x = c("multiple spectral mean", "SVD multiple spectral mean"), each = n)
VI_accuracy <- data.frame("index" = c(VI_spectral_multiple, svd_VI_spectral_multiple),
                          "noise" = rep(x = noise_seq, times = 2),
                          "algorithm" = algorithm)
RI_accuracy <- data.frame("index" = c(RI_spectral_multiple, svd_RI_spectral_multiple),
                          "noise" = rep(x = noise_seq, times = 2),
                          "algorithm" = algorithm)
ARI_accuracy <- data.frame("index" = c(ARI_spectral_multiple, svd_ARI_spectral_multiple),
                           "noise" = rep(x = noise_seq, times = 2),
                           "algorithm" = algorithm)

# algorithm = rep(x = c("kmeans_multiple", "kmeans_multiple_mean", "spectral_multiple", "spectral_multiple", "k means", "spectral", "cvx_bicluster"), each = n)
# VI_accuracy <- data.frame("index" = c(VI_kmeans_multiple, VI_kmeans_multiple_mean, VI_spectral_multiple, VI_spectral_multiple, VI_kmeans, VI_spectral, VI_cvx_bicluster),
#                           "noise" = rep(x = noise_seq, times = 7),
#                           "algorithm" = algorithm)
# RI_accuracy <- data.frame("index" = c(RI_kmeans_multiple, RI_kmeans_multiple_mean, RI_spectral_multiple, RI_spectral_multiple, RI_kmeans, RI_spectral, RI_cvx_cluster),
#                           "noise" = rep(x = noise_seq, times = 7),
#                           "algorithm" = algorithm)
# ARI_accuracy <- data.frame("index" = c(ARI_kmeans_multiple, ARI_kmeans_multiple_mean, ARI_spectral_multiple, ARI_spectral_multiple, ARI_kmeans, ARI_spectral, ARI_cvx_cluster),
#                           "noise" = rep(x = noise_seq, times = 7),
#                           "algorithm" = algorithm)

## Compare the performance of three algorithms under difference noise level
ggplot(data = VI_accuracy, aes(x= noise)) + geom_line(aes(y = index, col = algorithm), size = 1) +
  theme_bw() + ggtitle("VI for difference noise levels") + 
  theme(plot.title = element_text(size = 20, face = "bold", colour = "tomato", hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        plot.caption = element_text(size = 15),
        legend.title = element_text(size = 15, color = "firebrick"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Noise", y = "Variance of Information(VI)", caption = "Semi-biclustering")  


ggplot(data = RI_accuracy, aes(x= noise)) + geom_line(aes(y = index, col = algorithm), size = 1) +
  theme_bw() + ggtitle("RI for difference noise levels") + 
  theme(plot.title = element_text(size = 20, face = "bold", colour = "tomato", hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        plot.caption = element_text(size = 15),
        legend.title = element_text(size = 15, color = "firebrick"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Noise", y = "Rand Index(RI)", caption = "Semi-biclustering")  

ggplot(data = ARI_accuracy, aes(x= noise)) + geom_line(aes(y = index, col = algorithm), size = 1) +
  theme_bw() + ggtitle("ARI for difference noise levels") + 
  theme(plot.title = element_text(size = 20, face = "bold", colour = "tomato", hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        plot.caption = element_text(size = 15),
        legend.title = element_text(size = 15, color = "firebrick"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Noise", y = "Adjusted Rand Index(ARI)", caption = "Semi-biclustering")  

## Count the total runtime
end_time <- Sys.time()
run_time <- end_time - start_time

# ggsave("D:/OneDrive/Research/Master Thesis/simulation/noise1.eps", device = "eps")
save.image(file = "D:/OneDrive/Research/Master Thesis/simulation/RData/both_semi_bicluster_lung_cancer.RData")

