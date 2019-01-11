library(ggplot2)
library(reshape2)
library(gptk)

datdir<-"C:\\Users\\Wyss User\\Documents\\wyss\\dlt\\dat\\select_KI\\select\\"
#groups<-c("E. Coli", "E. Coli + Antibiotics", "E. Coli + Antibiotics + DLT", "E. Coli + DLT", "Saline Dextrose")
groups<-c("E. Coli", "E. Coli + DLT(1)", "E. Coli + DLT(2)", "E. Coli + DLT(3)", "Saline Dextrose")

all_data<-vector("list", length(groups))

for (group in 1:length(groups)){
  setwd(paste(datdir, groups[[group]], sep=""))
  list_of_dat<-list()
  list_of_bm<-list.files()
  for (i in 1:length(list_of_bm)){
    bm<-as.matrix(read.table(list_of_bm[i], header=FALSE, row.names=1, sep=","))
    # if (any(bm[,1]==0)){flag=FALSE}
    # else {flag=TRUE}
    # for (j in 2:dim(bm)[1]){
    #   div<-bm[j,1]
    #   if (flag){
    #     for (k in 1:((dim(bm)[2]))){
    #       bm[j,k]<-((bm[j,k]-div)/div)*100
    #     }
    #   }
    # }
    value<-as.vector(bm[2:dim(bm)[1],])
    time<-as.vector(t(matrix(t(bm[1,]),nrow=(dim(bm)[2]),ncol=dim(bm)[1]-1)))
    new_bm<-cbind(time, value)
    colnames(new_bm)<-c("Time", strsplit(list_of_bm[[i]],"_")[[1]][1])
    list_of_dat[[i]]<-new_bm
  }
  all_data[[group]]<-list_of_dat
}

mus<-list()
Ss<-list()

support <- as.matrix(seq(-4,8,0.1))
support_length <- length(support)

for (bm in 1:length(all_data[[1]])){
  mu <- array(dim=c(length(groups), support_length))
  S <- array(dim=c(length(groups), support_length))
  for (group in 1:length(groups)){
    x<-all_data[[group]][[bm]]
    options <- gpOptions(approx='ftc')
    options$kern$comp = list('rbf','white') 
    options$scale2var1 = TRUE
    gp_data<-na.omit(x)
    gp_model <- gpCreate(1, 1, as.matrix(gp_data[,1]), as.matrix(gp_data[,2]), options)
    gp_model <- gpOptimise(gp_model)
    post_mean_var<- gpPosteriorMeanVar(gp_model, support, varsigma.return=TRUE)
    mu[group,] <- post_mean_var$mu
    S[group,] <-post_mean_var$varsigma
  }
  print(bm)
  mus[[bm]]<-mu
  Ss[[bm]]<-S
}

for (bm in 1:length(all_data[[1]])){
  bm_name <- colnames(list_of_dat[[bm]])[[2]]
  mu <- mus[[bm]]
  S <- Ss[[bm]]
  data_to_plot <- as.data.frame(cbind(support,t(mu)))
  data_err_min <- c(t(mu-sqrt(abs(S))))
  data_err_max <- c(t(mu+sqrt(abs(S))))
  #colnames(data_to_plot) <- c("time", "E. coli", "E. coli + Antibiotics", "E. coli + Antibiotics + DLT", "E. Coli + DLT", "Saline Dextrose")
  colnames(data_to_plot) <- c("time", "E. coli", "E. Coli + DLT(1)", "E. Coli + DLT(2)", "E. Coli + DLT(3)", "Saline Dextrose")
  gp_plot <- ggplot((melt(data_to_plot, id='time')), aes(x=time, y=value, colour=variable)) + geom_line(lwd=2) + geom_ribbon(aes(ymin=data_err_min,ymax=data_err_max,fill=variable),alpha=0.10,linetype=0) + theme_classic()
  gp_plot <- gp_plot + labs(title=bm_name, x="Experimental Time (h)", y="Absolute Level of Biomarker", fill="Treatment Group", color="Treatment Group")
  gp_plot <- gp_plot + scale_x_continuous(breaks=seq(-4, 8, 1))
  ggsave(paste(datdir, "univariate_absolute_numfilters\\", bm_name, ".png", sep=""), gp_plot)
  #ggsave(paste(datdir, "univariate_absolute\\", bm_name, ".png", sep=""), gp_plot)
}

ellipse_angles <- as.matrix(seq(0,2*pi,(pi/64)))
colours<- c("#00BFC4", "#3399FF", "#7CAE00","#c77CFF", "#F8766D")
frame_times<-seq(-4,8,1)

groups_max<-list(c(3,2),c(3,4))
groups_min<-list()
metric <- lapply(frame_times, function(i) {matrix(0, nrow=length(all_data[[1]]), ncol=length(all_data[[1]]))})

for (bm1 in 1:(length(all_data[[1]])-1)){
  bm1_name <- colnames(list_of_dat[[bm1]])[[2]]
  mu1 <- mus[[bm1]]
  S1 <- sqrt(abs(Ss[[bm1]]))
  max1 = max(mu1+S1)
  min1 = min(mu1-S1)
  for (bm2 in (bm1+1):length(all_data[[1]])){
    print(bm1)
    print(bm2)
    bm2_name <- colnames(list_of_dat[[bm2]])[[2]]
    mu2 <- mus[[bm2]]
    S2 <- sqrt(abs(Ss[[bm2]]))
    max2 = max(mu2+S2)
    min2 = min(mu2-S2)
    for (t in 1:support_length){
      if (any(frame_times==support[t])){
        # metric[[which(frame_times==support[t])]][bm1,bm2] = Reduce('*', lapply(groups_max, function(g){((mu1[g[1],t]-mu1[g[2],t])**2/(S1[g[1],t]+S1[g[2],t])**2) * ((mu2[g[1],t]-mu2[g[2],t])**2/(S2[g[1],t]+S2[g[2],t])**2)})) / Reduce('*', lapply(groups_min, function(g){((mu1[g[1],t]-mu1[g[2],t])**2/(S1[g[1],t]+S1[g[2],t])**2) * ((mu2[g[1],t]-mu2[g[2],t])**2/(S2[g[1],t]+S2[g[2],t])**2)}))
        gp_plot <- ggplot() + geom_polygon(aes(x=S1[1,t]*cos(ellipse_angles)+mu1[1,t], y=S2[1,t]*sin(ellipse_angles)+mu2[1,t], colour=colours[[1]], fill=colours[[1]]), alpha=0.10, show.legend=FALSE)
        gp_plot <- gp_plot + geom_polygon(aes(x=S1[2,t]*cos(ellipse_angles)+mu1[2,t], y=S2[2,t]*sin(ellipse_angles)+mu2[2,t], colour=colours[[2]], fill=colours[[2]]), alpha=0.10, show.legend=FALSE)
        gp_plot <- gp_plot + geom_polygon(aes(x=S1[3,t]*cos(ellipse_angles)+mu1[3,t], y=S2[3,t]*sin(ellipse_angles)+mu2[3,t], colour=colours[[3]], fill=colours[[3]]), alpha=0.10, show.legend=FALSE)
        gp_plot <- gp_plot + geom_polygon(aes(x=S1[4,t]*cos(ellipse_angles)+mu1[4,t], y=S2[4,t]*sin(ellipse_angles)+mu2[4,t], colour=colours[[4]], fill=colours[[4]]), alpha=0.10, show.legend=FALSE)
        gp_plot <- gp_plot + geom_text(aes(x=mu1[1,t], y=mu2[1,t], label="EC", colour=colours[[1]]), show.legend=FALSE)
        #gp_plot <- gp_plot + geom_text(aes(x=mu1[2,t], y=mu2[2,t], label="EC+ABX", colour=colours[[2]]), show.legend=FALSE)
        #gp_plot <- gp_plot + geom_text(aes(x=mu1[3,t], y=mu2[3,t], label="EC+ABX+DLT", colour=colours[[3]]), show.legend=FALSE)
        #gp_plot <- gp_plot + geom_text(aes(x=mu1[4,t], y=mu2[4,t], label="EC+DLT", colour=colours[[4]]), show.legend=FALSE)
        gp_plot <- gp_plot + geom_text(aes(x=mu1[2,t], y=mu2[2,t], label="EC+DLT(1)", colour=colours[[2]]), show.legend=FALSE)
        gp_plot <- gp_plot + geom_text(aes(x=mu1[3,t], y=mu2[3,t], label="EC+DLT(2)", colour=colours[[3]]), show.legend=FALSE)
        gp_plot <- gp_plot + geom_text(aes(x=mu1[4,t], y=mu2[4,t], label="EC+DLT(3)", colour=colours[[4]]), show.legend=FALSE)
        gp_plot <- gp_plot + geom_text(aes(x=mu1[5,t], y=mu2[5,t], label="SD", colour=colours[[5]]), show.legend=FALSE)
        gp_plot <- gp_plot + coord_cartesian(xlim=c(min1,max1), ylim=c(min2,max2))
        gp_plot <- gp_plot + labs(title=paste("Experimental Time", support[t], "(h)"), x=bm1_name, y=bm2_name) +theme_classic()
        #ggsave(paste(datdir, "bivariate_absolute\\", bm1_name, "_", bm2_name, "_t", t, ".png", sep=""), gp_plot)
        ggsave(paste(datdir, "bivariate_absolute_numfilters\\", bm1_name, "_", bm2_name, "_t", t, ".png", sep=""), gp_plot)
      }
    }
  }
}

for (t_ in 1:length(frame_times)){
  idx <- which(metric[[t_]]==max(as.vector(metric[[t_]])), arr.ind=TRUE)
  t <- which(support==frame_times[t_])
  print(idx)
  bm1 <- idx[1,1]
  bm2 <- idx[1,2]
  bm1_name <- colnames(list_of_dat[[bm1]])[[2]]
  mu1 <- mus[[bm1]]
  S1 <- sqrt(abs(Ss[[bm1]]))
  max1 = max(mu1+S1)
  min1 = min(mu1-S1)
  bm2_name <- colnames(list_of_dat[[bm2]])[[2]]
  mu2 <- mus[[bm2]]
  S2 <- sqrt(abs(Ss[[bm2]]))
  max2 = max(mu2+S2)
  min2 = min(mu2-S2)
  gp_plot <- ggplot() + geom_polygon(aes(x=S1[1,t]*cos(ellipse_angles)+mu1[1,t], y=S2[1,t]*sin(ellipse_angles)+mu2[1,t], colour=colours[[1]], fill=colours[[1]]), alpha=0.10, show.legend=FALSE)
  gp_plot <- gp_plot + geom_polygon(aes(x=S1[2,t]*cos(ellipse_angles)+mu1[2,t], y=S2[2,t]*sin(ellipse_angles)+mu2[2,t], colour=colours[[2]], fill=colours[[2]]), alpha=0.10, show.legend=FALSE)
  gp_plot <- gp_plot + geom_polygon(aes(x=S1[3,t]*cos(ellipse_angles)+mu1[3,t], y=S2[3,t]*sin(ellipse_angles)+mu2[3,t], colour=colours[[3]], fill=colours[[3]]), alpha=0.10, show.legend=FALSE)
  gp_plot <- gp_plot + geom_polygon(aes(x=S1[4,t]*cos(ellipse_angles)+mu1[4,t], y=S2[4,t]*sin(ellipse_angles)+mu2[4,t], colour=colours[[4]], fill=colours[[4]]), alpha=0.10, show.legend=FALSE)
  gp_plot <- gp_plot + geom_text(aes(x=mu1[1,t], y=mu2[1,t], label="EC", colour=colours[[1]]), show.legend=FALSE)
  gp_plot <- gp_plot + geom_text(aes(x=mu1[2,t], y=mu2[2,t], label="EC+ABX", colour=colours[[2]]), show.legend=FALSE)
  gp_plot <- gp_plot + geom_text(aes(x=mu1[3,t], y=mu2[3,t], label="EC+ABX+DLT", colour=colours[[3]]), show.legend=FALSE)
  gp_plot <- gp_plot + geom_text(aes(x=mu1[4,t], y=mu2[4,t], label="SD", colour=colours[[4]]), show.legend=FALSE)
  gp_plot <- gp_plot + coord_cartesian(xlim=c(min1,max1), ylim=c(min2,max2))
  gp_plot <- gp_plot + labs(title=paste("Experimental Time", support[t], "(h)"), x=bm1_name, y=bm2_name) +theme_classic()
  ggsave(paste(datdir, "bivariate_top\\t", t, "_", bm1_name, "_", bm2_name, ".png", sep=""), gp_plot)
}

metric_overall <- Reduce('+', metric)
















x<-all_data[[1]][[1]]

options <- gpOptions(approx='ftc')
options$kern$comp = list('rbf','white') 
options$scale2var1 = TRUE

support <- as.matrix(seq(-4,8,0.1))
support_length <- length(support)
mu <- array(dim=c(4, support_length))
S <- array(dim=c(4, support_length))

label_idx <- 4

gp_data<-na.omit(x)
gp_model <- gpCreate(1, 1, as.matrix(gp_data[,1]), as.matrix(gp_data[,2]), options)
gp_model <- gpOptimise(gp_model)
post_mean_var<- gpPosteriorMeanVar(gp_model, support, varsigma.return=TRUE)
mu[label_idx,] <- post_mean_var$mu
S[label_idx,] <-post_mean_var$varsigma

num_sigmas <- 1
data_to_plot <- as.data.frame(cbind(support,t(mu)))
data_err_min <- c(t(mu-num_sigmas*sqrt(abs(S))))
data_err_max <- c(t(mu+num_sigmas*sqrt(abs(S))))
colnames(data_to_plot) <- c("time", "E. coli", "E. coli + Antibiotics", "E. coli + Antibiotics + DLT", "Saline Dextrose")

gp_plot <- ggplot((melt(data_to_plot, id='time')), aes(x=time, y=value, colour=variable)) + geom_line(lwd=2) + geom_ribbon(aes(ymin=data_err_min,ymax=data_err_max,fill=variable),alpha=0.10,linetype=0) + theme_classic()
gp_plot <- gp_plot + labs(title="Combination of Antibiotic treatment with DLT reduces SOFA score", x="Experimental Time (h)", y="SOFA Score", fill="Treatment Group", color="Treatment Group")
gp_plot <- gp_plot + scale_x_continuous(breaks=seq(-4, 8, 1)) + scale_y_continuous(breaks=seq(0, 8, 1)) 
print(gp_plot)
ggsave(paste("dlt_plot_gpmodel.png", sep=""), gp_plot)

means <- array(dim=c(4, 13))
vars <- array(dim=c(4, 13))

label_idx <- 4
x <- dlt_control

means[label_idx,] <- rowMeans(as.matrix(x[,2:dim(x)[2]]), na.rm=TRUE)
vars[label_idx,] <- apply(as.matrix(x[,2:dim(x)[2]]), 1, var, na.rm=TRUE)

num_sigmas <- 1
data_to_plot <- as.data.frame(cbind(seq(-4,8,1),t(means)))
data_err_min <- c(t(means-num_sigmas*sqrt(abs(vars))))
data_err_max <- c(t(means+num_sigmas*sqrt(abs(vars))))
colnames(data_to_plot) <- c("time","E. coli", "E. coli + Antibiotics", "E. coli + Antibiotics + DLT", "Saline Dextrose")

gp_plot <- ggplot((melt(data_to_plot, id='time')), aes(x=time, y=value, colour=variable)) + geom_line(lwd=2) + geom_ribbon(aes(ymin=data_err_min,ymax=data_err_max,fill=variable),alpha=0.10,linetype=0) + theme_classic()
gp_plot <- gp_plot + labs(title="Combination of Antibiotic treatment with DLT reduces SOFA score", x="Experimental Time (h)", y="SOFA Score", fill="Treatment Group", color="Treatment Group")
gp_plot <- gp_plot + scale_x_continuous(breaks=seq(-4, 8, 1)) + scale_y_continuous(breaks=seq(0, 8, 1)) 
print(gp_plot)
ggsave(paste("dlt_plot_regular.png", sep=""), gp_plot)

print(means)
print(vars)
