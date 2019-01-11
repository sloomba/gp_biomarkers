library(gptk)
library(ggplot2)
library(reshape2)
library(infotheo)
library(files2)
#features: (sample, label, time, gene1, gene2, ...)
normalize <- FALSE
num_sigmas <- 2
support_length <- 100
num_bins <- sqrt(support_length)
top_k <- 500
cutoff <- 9

num_samples <- nrow(features)
num_feats <- ncol(features)-3 #first three are sample ID, labels and time
sample_names <- sort(unique(features[,1]))
label_names <- sort(unique(features[,2]))
time_points <- sort(unique(features[,3]))
feat_names <- row.names(features)
num_labels <-length(label_names)

options <- gpOptions(approx='ftc')
options$kern$comp = list('rbf','white') 
options$scale2var1 = TRUE
#options$optimiseBeta = TRUE
#options$beta = 1e+6


### univariate GP model

### 2D analysis
support <- as.matrix(seq(min(time_points),max(time_points),length.out=support_length))
mu <- array(dim=c(num_labels, num_feats, support_length))
S <- array(dim=c(num_labels, num_feats, support_length))

for (label_idx in 1:num_labels){
  label <- label_names[label_idx]
  data_of_interest <- features[features$label==label,]
  if (normalize){
    for (sample in unique(data_of_interest$sample)){
      min_t <- min(data_of_interest[data_of_interest$sample==sample,3])
      data_of_interest[data_of_interest$sample==sample, 4:num_feats+3] <- data_of_interest[data_of_interest$sample==sample, 4:num_feats+3]/data_of_interest[data_of_interest$sample==sample & data_of_interest$time==min_t, 4:num_feats+3] - 1
    }
  }
  for (feat in 1:num_feats){
    gp_model <- gpCreate(1,1, as.matrix(data_of_interest[,3]), as.matrix(data_of_interest[,feat+3]), options)
    gp_model <- gpOptimise(gp_model)
    post_mean_var<- gpPosteriorMeanVar(gp_model, support, varsigma.return=TRUE)
    mu[label_idx,feat,] <- post_mean_var$mu
    S[label_idx,feat,] <- post_mean_var$varsigma
  }
}
#par(ask=TRUE)
for (i in 1:num_feats){
  data_to_plot <- as.data.frame(cbind(support,t(mu[,i,])))
  data_err_min <- c(t(mu[,i,]-num_sigmas*sqrt(abs(S[,i,]))))
  data_err_max <- c(t(mu[,i,]+num_sigmas*sqrt(abs(S[,i,]))))
  colnames(data_to_plot) <- c("time",label_names)
  ggg2 <- ggplot((melt(data_to_plot, id='time')), aes(x=time, y=value, colour=variable)) + geom_line(lwd=2) + geom_ribbon(aes(ymin=data_err_min,ymax=data_err_max,fill=variable),alpha=0.10,linetype=0) + theme_classic() + labs(title=paste("gene",i,"|",feat_names[i+3]), x="time", y="expression")
  ggsave(paste("figures/all/gene_",i,"_",gsub("/", ".", feat_names[i+3]),".png", sep=""), ggg2)
  #print(ggg2)
}

####ANALYSIS STUFF

mi_inter_matrix <- array(dim=c(num_labels, num_feats, num_feats))
co_inter_matrix <- array(dim=c(num_labels, num_feats, num_feats))
eu_inter_matrix <- array(dim=c(num_labels, num_feats, num_feats))
mi_intra_matrix <- array(dim=c(num_feats, num_labels, num_labels))
co_intra_matrix <- array(dim=c(num_feats, num_labels, num_labels))
eu_intra_matrix <- array(dim=c(num_feats, num_labels, num_labels))

for (i in 1:num_feats){
  print(i)
  mu_curr <- t(mu[,i,])
  mi_intra_matrix[i,,] <- mutinformation(discretize(mu_curr,disc="globalequalwidth",nbins=num_bins), method="mm")
  co_intra_matrix[i,,] <- cor(mu_curr)
  eu_intra_matrix[i,,] <- as.matrix(dist(t(mu_curr)))
}

for (label_idx in 1:num_labels){
  print(label_idx)
  mu_curr <- t(mu[label_idx,,])
  mi_inter_matrix[label_idx,,] <- mutinformation(discretize(mu_curr,disc="globalequalwidth",nbins=num_bins), method="mm")
  co_inter_matrix[label_idx,,] <- cor(mu_curr)
  eu_inter_matrix[label_idx,,] <- as.matrix(dist(t(mu_curr)))
}

feat_ranks <- vector(length=num_feats)
for (i in 1:num_feats){
  print(i)
  combo_score <- eu_intra_matrix[i,,]#*(1+mi_intra_matrix[i,,])#eu_intra_matrix[i,,]#co_intra_matrix[i,,]#(1-co_intra_matrix[i,,])*mi_intra_matrix[i,,]*log1p(eu_intra_matrix[i,,])
  combo_score <- combo_score[lower.tri(combo_score)]
  feat_ranks[i] <- sum(combo_score)
}

feat_ranks_euclid_file <- "docs/gene_ranks_euclidean.csv"
feat_ranks_combo_file <- "docs/gene_ranks_combo.csv"
write.csv(feat_ranks, feat_ranks_euclid_file)

feat_ranks_sorted<-sort(feat_ranks,decreasing=TRUE,index.return=TRUE)

top_gene_euclid_dir <- "figures/top_genes_euclidean"
top_gene_combo_dir <- "figures/top_genes_combo"
if (dir.exists(top_gene_euclid_dir)){
  unlink(top_gene_euclid_dir, recursive=TRUE)
}
dir.create(top_gene_euclid_dir)
if (dir.exists(top_gene_combo_dir)){
  unlink(top_gene_combo_dir, recursive=TRUE)
}
dir.create(top_gene_combo_dir)

#par(ask=TRUE)
rank <- 1
for (i in feat_ranks_sorted$ix[1:min(c(top_k,length(feat_ranks_sorted$ix)))]){
  data_to_plot <- as.data.frame(cbind(support,t(mu[,i,])))
  data_err_min <- c(t(mu[,i,]-num_sigmas*sqrt(abs(S[,i,]))))
  data_err_max <- c(t(mu[,i,]+num_sigmas*sqrt(abs(S[,i,]))))
  colnames(data_to_plot) <- c("time",label_names)
  ggg2 <- ggplot((melt(data_to_plot, id='time')), aes(x=time, y=value, colour=variable)) + geom_line(lwd=2) + geom_ribbon(aes(ymin=data_err_min,ymax=data_err_max,fill=variable),alpha=0.10,linetype=0) + theme_classic() + labs(title=paste("gene",i,"|",feat_names[i+3]), x="time", y="expression")
  ggsave(paste(top_gene_combo_dir,"/rank_",rank,"_gene_",i,"_",gsub("/", ".", feat_names[i+3]),".png", sep=""), ggg2)
  rank <- rank+1
  #print(ggg2)
}

new_feat_matrix <- matrix(nrow=num_feats, ncol=num_labels*(num_labels-1))
for (i in 1:num_feats){
  temp1 <- mi_intra_matrix[i,,]
  temp1 <- temp1[lower.tri(temp1)]
  temp2 <- eu_intra_matrix[i,,]
  temp2 <- temp2[lower.tri(temp2)]
  if (max(temp2)>0){
    temp2 <- temp2/max(temp2)
  }
  new_feat_matrix[i,] <- c(temp1,temp2)
}
dist_new_feat_matrix_file <- "docs/gene_distances_2d.csv"
dist_new_feat_matrix <- dist(new_feat_matrix)
write.csv(as.matrix(dist_new_feat_matrix), dist_new_feat_matrix_file)
gene_clusters <- hclust(dist_new_feat_matrix)
gene_clusters_cut <- cutree(gene_clusters, k=2)#h=cutoff)
cluster_ids <- sort(unique(gene_clusters_cut))
gene_clusters_file <- paste("docs/gene_clusters_clusters_2.csv",sep="")
write.csv(gene_clusters_cut, gene_clusters_file)

gene_clusters_dir <- paste("figures/gene_clusters_cutoff_",cutoff,sep="")
if (dir.exists(gene_clusters_dir)){
  unlink(gene_clusters_dir, recursive=TRUE)
}
dir.create(gene_clusters_dir)
for (id in cluster_ids){
  feats_to_plot <- (1:num_feats)[gene_clusters_cut==id]
  dir.create(paste(gene_clusters_dir,"/cluster_",id,sep=""))
  for (i in feats_to_plot){
    data_to_plot <- as.data.frame(cbind(support,t(mu[,i,])))
    data_err_min <- c(t(mu[,i,]-num_sigmas*sqrt(abs(S[,i,]))))
    data_err_max <- c(t(mu[,i,]+num_sigmas*sqrt(abs(S[,i,]))))
    colnames(data_to_plot) <- c("time",label_names)
    ggg2 <- ggplot((melt(data_to_plot, id='time')), aes(x=time, y=value, colour=variable)) + geom_line(lwd=2) + geom_ribbon(aes(ymin=data_err_min,ymax=data_err_max,fill=variable),alpha=0.10,linetype=0) + theme_classic() + labs(title=paste("gene",i,"|",feat_names[i+3]), x="time", y="expression")
    ggsave(paste(gene_clusters_dir,"/cluster_",id,"/gene_",i,"_",gsub("/", ".", feat_names[i+3]),".png", sep=""), ggg2)
    #print(ggg2)
  }
}

### bivariate GP model



### 3D analysis
support_3d_1 <- seq(min(label_names),max(label_names),length.out=support_length)
support_3d_2 <- seq(min(time_points),max(time_points),length.out=support_length)
support_3d <- matrix(nrow=support_length^2, ncol=2)
dummy <- 1
for (i in support_3d_1){
  support_3d[seq((dummy-1)*support_length+1,dummy*support_length),] <- t(rbind(seq(i,i,length.out=support_length),support_3d_2))
  dummy <- dummy+1
}

mu_3d <- array(dim=c(num_feats, support_length^2))
S_3d <- array(dim=c(num_feats, support_length^2))
gpmodels_3d <- list()

data_of_interest <- features
if (normalize){
  for (sample in unique(data_of_interest$sample)){
    min_t <- min(data_of_interest[data_of_interest$sample==sample,3])
    data_of_interest[data_of_interest$sample==sample, 4:num_feats+3] <- data_of_interest[data_of_interest$sample==sample, 4:num_feats+3]/data_of_interest[data_of_interest$sample==sample & data_of_interest$time==min_t, 4:num_feats+3] - 1
  }
}

for (feat in 1:num_feats){
  print(feat)
  gp_model <- gpCreate(2,1, as.matrix(data_of_interest[,2:3]), as.matrix(data_of_interest[,feat+3]), options)
  gp_model <- gpOptimise(gp_model)
  post_mean_var_3d<- gpPosteriorMeanVar(gp_model, support_3d, varsigma.return=TRUE)
  gpmodels_3d <- c(gpmodels_3d, gp_model)
  mu_3d[feat,] <- post_mean_var_3d$mu
  S_3d[feat,] <- post_mean_var_3d$varsigma
}

#par(ask=TRUE)
for (i in 1:num_feats){
  data_to_plot <- as.data.frame(cbind(support_3d,mu_3d[i,]))
  colnames(data_to_plot) <- c("level","time","expression")
  ggg2 <- ggplot(data_to_plot, aes(x=time, y=level, z=expression)) +  geom_raster(aes(fill = expression)) + scale_fill_distiller(palette = "Spectral") + theme_classic() + labs(title=paste("gene",i,"|",feat_names[i+3]), x="time", y="pathogen level")
  ggsave(paste("figures/all_genes_3d/gene_",i,"_",gsub("/", ".", feat_names[i+3]),".png", sep=""), ggg2)
  #print(ggg2)
}

param_matrix <- array(dim=c(num_feats, num_samples+3))

for (i in 1:num_feats){
  print(i)
  param_matrix[i,1:num_samples] <- gpmodels_3d[[28*(i-1)+18]]
  param_matrix[i,num_samples+1] <- gpmodels_3d[[28*(i-1)+20]]$comp[[1]]$inverseWidth
  param_matrix[i,num_samples+2] <- gpmodels_3d[[28*(i-1)+20]]$comp[[1]]$variance
  param_matrix[i,num_samples+3] <- gpmodels_3d[[28*(i-1)+20]]$comp[[2]]$variance
}

dist_param_matrix_file <- "docs/gene_distances_3d.csv"
dist_param_matrix <- dist(param_matrix)
write.csv(as.matrix(dist_param_matrix), dist_param_matrix_file)
gene_clusters_3d <- hclust(dist_param_matrix)
gene_clusters_3d_cut <- cutree(gene_clusters_3d, h=cutoff)
cluster_3d_ids <- sort(unique(gene_clusters_3d_cut))
gene_clusters_3d_file <- paste("docs/gene_clusters_3d_cutoff_",cutoff,".csv",sep="")
write.csv(gene_clusters_3d_cut, gene_clusters_3d_file)

gene_clusters_3d_dir <- paste("figures/gene_clusters_3d_cutoff_",cutoff,sep="")
if (dir.exists(gene_clusters_3d_dir)){
  unlink(gene_clusters_3d_dir, recursive=TRUE)
}
dir.create(gene_clusters_3d_dir)
for (id in cluster_3d_ids){
  feats_to_plot <- (1:num_feats)[gene_clusters_3d_cut==id]
  dir.create(paste(gene_clusters_3d_dir,"/cluster_",id,sep=""))
  for (i in feats_to_plot){
    data_to_plot <- as.data.frame(cbind(support_3d,mu_3d[i,]))
    colnames(data_to_plot) <- c("level","time","expression")
    ggg2 <- ggplot(data_to_plot, aes(x=time, y=level, z=expression)) +  geom_raster(aes(fill = expression)) + scale_fill_distiller(palette = "Spectral") + theme_classic() + labs(title=paste("gene",i,"|",feat_names[i+3]), x="time", y="pathogen level")
    ggsave(paste(gene_clusters_3d_dir,"/cluster_",id,"/gene_",i,"_",gsub("/", ".", feat_names[i+3]),".png", sep=""), ggg2)
    #print(ggg2)
  }
}