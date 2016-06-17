#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Clustering    +++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### performs 2-step clustering on expression data or background data
### 1. Step: clustering genes based on correlation to compare the shape
### 2. Step: split the clusters from step 1 into finer clusters based on the
###          eucledian distance to compare the magnitude

# INPUT:
#   data_table: (Numeric matrix genes x samples)
#   annotation_table: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric numbers.
#   control_timecourse: (bool) [Default FALSE] control time timecourse is 
#       part of the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in annotation_table.
#   case_name: (str) name of the case condition in annotation_table.
#   plot_clusters: (bool) [Default TRUE] plot the clusters
#   n_cluster: (scalar) 
#   n_genes: (scalar) 
# OUTPUT:
#   results: (list ["kmeans_clus","cluster_means", "n_pre_clus","fine_clus"] 
#       x number of runs [combined, case, control])
#       kmeans_clus: (vector number of genes) [1, number of clusters] 
#           Indicates assignment (value of function Kmeans) 
#       cluster_means: (matrix number of clusters x gene dimensions 
#           [number of samples per run, i.e. in case])
#           Centroids of fine (2nd) clustering.
#       n_pre_clus: (scalar) Number of clusters used in pre-clustering
#       n_fine_clusters: (scalar) Number of clusters used in fine-clustering

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

cluster_genes_for_impulse <- function(data_table, data_annotation,
                                      control_timecourse = FALSE, control_name = NULL, plot_clusters = FALSE,
                                      no_of_clusters = NULL, n_genes = NULL){
  
  #' @import amap
     corr_thres = 0.85
     eucl_thres = 1.8
     ### package for Kmeans clustering method, which allows correlation as distance
     
     ### assume only one timecourse as default and perform only 1 run
     results = list()
     runs = 1
     label = "case"
     
     ### perform 3 runs if control timecourse data is present
     if(control_timecourse == TRUE & is.null(control_name) == FALSE){
       runs = 3
       label = c("combined","case","control")
     }
     ### clustering for different runs
     # This for loops over the entire rest of the function
     for (c_runs in 1:runs){
       # Select target columns of expression matrix and retain only variable gene:
       # Cut off: sd(x)/(mean(x) > 0.025
       if(c_runs == 1){                # combined or only case if no control
         print(paste("Clustering of ",label[c_runs]," data set", sep = ""))
         dat_pre = data_table
         dat_pre = dat_pre[which(apply(dat_pre,1,function(x){100*
                                                               sd(x)/(mean(x))})>= 2.5),]
         dat = t(apply(dat_pre,1,function(y){unlist(lapply(as.list(unique(
           as.numeric(as.character(data_annotation$Time)))),function(x){
             mean(y[data_annotation$Time == x])}))}))
         colnames(dat) <- unique(as.numeric(as.character(data_annotation$Time)))
       } else if(c_runs == 2){         # case
         print(paste("Clustering of ",label[c_runs]," data set", sep = ""))
         dat_pre = data_table[,!(data_annotation$Condition %in% control_name)]
         dat_pre = dat_pre[which(apply(dat_pre,1,function(x){100*
                                                               sd(x)/(mean(x))})>= 2.5),]
         dat = t(apply(dat_pre,1,function(y){unlist(lapply(as.list(unique(
           as.numeric(as.character(data_annotation[!(data_annotation$Condition
                                                     %in% control_name),"Time"])))),function(x){mean(
                                                       y[data_annotation[!(data_annotation$Condition %in% control_name),
                                                                         "Time"] == x])}))}))
         colnames(dat) <- unique(as.numeric(as.character(data_annotation[
           !(data_annotation$Condition %in% control_name),"Time"])))
       } else {                        # control
         print(paste("Clustering of ",label[c_runs]," data set", sep = ""))
         dat_pre = data_table[,data_annotation$Condition %in% control_name]
         dat_pre = dat_pre[which(apply(dat_pre,1,function(x){100*
                                                               sd(x)/(mean(x))})>= 2.5),]
         dat = t(apply(dat_pre,1,function(y){unlist(lapply(as.list(unique(
           as.numeric(as.character(data_annotation[data_annotation$Condition
                                                   %in% control_name,"Time"])))),function(x){mean(
                                                     y[data_annotation[data_annotation$Condition %in% control_name,
                                                                       "Time"] == x])}))}))
         colnames(dat) <- unique(as.numeric(as.character(data_annotation[
           (data_annotation$Condition %in% control_name),"Time"])))
       }
       print(paste("- ", (nrow(data_table) - nrow(dat_pre)),
                   " genes were excluded due to very low variation", sep = ""))
       
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       #~~~~~~~~~~~~~   Step 1: clustering based on correlation   ~~~~~~~~~~~~~~~~~#
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       
       max_n_pre_clus = 25   # restrict maximum number of pre-clusters to 25
       ### z-transformation of data for clustering
       dat_ztrans <- t(scale(t(dat)))
       
       ### ---> clustering of expression data
       if(is.null(no_of_clusters)){
         file_part = "genes"
         n_clus = 0
         for(i in 1:max_n_pre_clus){
           
           kmeans_clus <- Kmeans(dat_ztrans,i,200,method="correlation")
           
           ### accept clustering if mean correlation between genes of the cluster and
           ### centroid of the cluster is higher than 'corr_thres' (default: 0.75)
           for(j in 1:i){
             tmp <- dat_ztrans[kmeans_clus$cluster %in% j,]
             indt1 <- NULL
             if(i > 1 & sort(summary(as.factor(kmeans_clus$cluster)))[1]  < 10){
               indt1 <- "I"
               n_clus = i -1
               kmeans_clus <- kmeans_clus_pre
               break
             }
             if(is.null(nrow(tmp))){    # if cluster contains only 1 gene
               average_pearson <- cor(tmp, kmeans_clus$centers[j,])
               #             min_pearson <- average_pearson
             } else {                   # if cluster contains > 1 genes
               average_pearson <- mean(apply(tmp,1,function(x){cor(x,
                                                                   kmeans_clus$centers[j,])}))
               #             min_pearson <- min(apply(tmp,1,function(x){cor(x,
               #                kmeans_clus$centers[j,])}))
             }
             ### if one cluster does not fufill the correlation condition leave loop,
             ### increase numbers of potential clusters and try again
             if(average_pearson < corr_thres){
               break
             }
             ### accept clustering and leave loop if all clusters fufill condition
             if(j == i) {
               n_clus = i
               break
             }
           }
           if(is.null(indt1) == FALSE){
             break
           }
           if(n_clus !=0){
             break
             ### if maximum number of allowed clusters is reached accept this number
           } else if((n_clus == 0) & (i == max_n_pre_clus)){
             n_clus = i
             break
           }
           kmeans_clus_pre = kmeans_clus
         }
         print(paste("--- Number of correlation-based pre-clusters: ",
                     n_clus,sep=""))
         
         ### clustering of background data based on clusters from expression data
       } else {
         file_part = "bg"
         n_clus <- max(round(no_of_clusters[[c_runs*2 - 1]]*nrow(data_table)/
                               n_genes),1)
         kmeans_clus <- Kmeans(dat_ztrans,n_clus,200,method="correlation")
         n_fine_clust_rand <- round(sample( no_of_clusters[[c_runs*2]] ,n_clus,
                                            replace = T)*nrow(data_table)/n_genes)
         n_fine_clust_rand[n_fine_clust_rand < 1] = 1
       }
       kmeans_pre_clus_corr <- kmeans_clus
       n_pre_clus <- n_clus
       kmeans_clus_final <- kmeans_pre_clus_corr$cluster
       #    write.table(kmeans_clus_final,paste("kmeans_pre_clus_corr_",label[1],
       #        "_",file_part,".txt",sep=""), sep="\t",quote=FALSE)
       
       print(paste("------ Number of genes in pre-clusters: ",
                   paste(paste("C",names(summary(as.factor(kmeans_clus_final))),": ",
                               summary(as.factor(kmeans_clus_final)), sep=""),collapse=", "),sep=""))
       
       
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       #~~~~~~~~~   Step 2: clustering based on eucledian distance   ~~~~~~~~~~~~~~#
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       
       
       n_fine_clusters <- NULL
       cluster_means <-  NULL
       # restrict maximum number of fine clusters to 5 per pre-cluster
       max_n_fine_clus <- 10
       kmeans_clus <- NULL
       
       tempo2 <- do.call(c,as.data.frame(dat))
       tempo3 <- (tempo2 - median(tempo2))/mad(tempo2)
       dat_median_ztrans <- matrix(tempo3, nrow(dat), ncol(dat))
       rownames(dat_median_ztrans) = rownames(dat)
       colnames(dat_median_ztrans) = colnames(dat)
       
       ### cluster each pre-cluster into fine clusters
       for(cl in 1:n_pre_clus){
         
         
         ### ---> clustering of expression data
         if(is.null(no_of_clusters)){
           if(length(kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster %in%
                                                    cl]) >  max_n_fine_clus){
             n_clus = 0
             
             for(i in 1:max_n_fine_clus){
               
               kmeans_clus <- kmeans(dat_median_ztrans[names(
                 kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                                              %in% cl]),],i,200)
               ### accept clustering if mean eucledian distance between genes of
               ### the cluster and centroid of the cluster is lower than
               ### 'eucl_thres' (default: 5)
               for(j in 1:i){
                 indt2 <- NULL
                 tmp <- dat_median_ztrans[names(kmeans_pre_clus_corr$cluster[
                   kmeans_pre_clus_corr$cluster %in% cl]),][kmeans_clus$cluster
                                                            %in% j,]
                 # if smallest cluster would have less than 5 genes take previous
                 # number of clusters
                 if(i > 1 & sort(summary(as.factor(kmeans_clus$cluster)))[1]  < 5){
                   indt2 = "I"
                   n_clus = i - 1
                   kmeans_clus <- kmeans_clus_pre
                   break
                 }
                 
                 if(is.null(nrow(tmp))){    # if cluster contains only 1 gene
                   #               average_eucledian <- dist(rbind(tmp, kmeans_clus$centers[j,]))
                   average_eucledian <- 0
                   #                   min_eucledian <- average_eucledian
                 } else {                   # if cluster contains > 1 genes
                   #                   average_eucledian <- mean(apply(tmp,1,function(x){dist(
                   #                       rbind(x,kmeans_clus$centers[j,]))}))
                   average_eucledian <- sd(apply(tmp,1,function(x){dist(rbind(x,
                                                                              kmeans_clus$centers[j,]))}))
                   #                   min_eucledian <- min(apply(tmp,1,function(x){dist(rbind(x,
                   #                      kmeans_clus$centers[j,]))}))
                 }
                 
                 ### if one cluster does not fufill the eucledian condition leave
                 ### loop,increase numbers of potential clusters and try again
                 if(average_eucledian > eucl_thres){
                   break
                 }
                 ### accept clustering and leave loop if all clusters fulfill
                 ### condition
                 if(j == i) {
                   n_clus = i
                   break
                 }
               }
               if(is.null(indt2) == FALSE){
                 break
               }
               if(n_clus !=0){
                 break
                 ### if maximum number of allowed clusters is reached accept this
                 ### number
               } else if((n_clus == 0) & (i == max_n_fine_clus)){
                 n_clus = i
                 break
               }
               kmeans_clus_pre = kmeans_clus
             }
             n_fine_clusters <- c(n_fine_clusters,n_clus)
           } else {
             n_clus = 1
             n_fine_clusters <- c(n_fine_clusters,n_clus)
             kmeans_clus$cluster <- rep(1,length(kmeans_pre_clus_corr$cluster[
               kmeans_pre_clus_corr$cluster %in% cl]))
             names(kmeans_clus$cluster) <- names(kmeans_pre_clus_corr$cluster[
               kmeans_pre_clus_corr$cluster %in% cl])
           }
           ### clustering of background data based on clusters from expression data
         } else {
           if(floor(length(kmeans_pre_clus_corr$cluster[
             kmeans_pre_clus_corr$cluster %in% cl]) / n_fine_clust_rand[cl])
             >=2){
             n_clus <- n_fine_clust_rand[cl]
             kmeans_clus <- kmeans(dat_median_ztrans[names(
               kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                                            %in% cl]),], n_clus,200)
           } else if(n_fine_clust_rand[cl] > 1 & floor(length(
             kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster %in%
                                            cl]) / (n_fine_clust_rand[cl]-1)) >=2){
             n_clus <- n_fine_clust_rand[cl] - 1
             kmeans_clus <- kmeans(dat_median_ztrans[names(
               kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                                            %in% cl]),], n_clus,200)
           } else {
             n_clus = 1
             kmeans_clus$cluster <- rep(1,length(
               kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                                            %in% cl]))
             names(kmeans_clus$cluster) <- names(kmeans_pre_clus_corr$cluster[
               kmeans_pre_clus_corr$cluster %in% cl])
           }
           n_fine_clusters <- c(n_fine_clusters,n_clus)
         }
         
         #      ### calculate mean values for each final cluster
         for (i in 1:n_clus){
           tmp <- dat_pre[names(kmeans_clus$cluster)[kmeans_clus$cluster %in% i],]
           if(is.null(dim(tmp))){       # if cluster contains only 1 gene
             cluster_means <- rbind(cluster_means,tmp)
           } else {                   # if cluster contains > 1 genes
             cluster_means <- rbind(cluster_means,colMeans(tmp))
           }
         }
         ### overwrite pre-clusters by fine clusters to have cluster information
         ### for each gene
         if(cl==1){
           kmeans_clus_final[names(kmeans_pre_clus_corr$cluster[
             kmeans_pre_clus_corr$cluster %in% cl])] <-  kmeans_clus$cluster
         } else {
           kmeans_clus_final[names(kmeans_pre_clus_corr$cluster[
             kmeans_pre_clus_corr$cluster %in% cl])] <-  kmeans_clus$cluster +
             sum(n_fine_clusters[1:(cl-1)])
         }
         #      write.table(kmeans_clus$cluster,file = paste("kmeans_clus_",cl,"_",
         #            label[1],"_",file_part, ".txt",sep=""),sep="\t",quote=FALSE)
       }
       
       if(file_part != "bg"){
         write.table(kmeans_clus_final,paste("kmeans_clus_final_",label[c_runs],
                                             "_",file_part,".txt",sep=""),sep="\t",quote = FALSE,
                     row.names = TRUE, col.names = FALSE)
       }
       rownames(cluster_means)  <- 1:sum(n_fine_clusters)
       n_clus = sum(n_fine_clusters)
       
       print(paste("--- Final number of clusters: ",n_clus,sep=""))
       print(paste("------ Number of genes in final clusters: ",
                   paste(paste("C",names(summary(as.factor(kmeans_clus_final))),": ",
                               summary(as.factor(kmeans_clus_final)), sep=""),collapse=", "),sep=""))
       
       
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       #~~~~~~~~~~~~~~~~   Plot clusters if plot_clusters == TRUE ~~~~~~~~~~~~~~~~~#
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       
       #    if(plot_clusters == TRUE & is.null(no_of_clusters)){
       if(plot_clusters == TRUE & file_part != "bg"){
         count = 0
         pdf(paste("clusters_",label[c_runs],"_",file_part,".pdf",sep=""),
             height=10.0,width=16.0)
         
         for(i in 1:n_clus){
           
           ### 1 or 2 clusters
           if((i == 1) & (n_clus <= 2)){
             par(mfrow=c(1,2))
           } else if((i == 1) & (n_clus <= 4)){
             par(mfrow=c(2,2))
             ### print all clusters on one page if number of clusters is <= 12
           } else if((i == 1) & (n_clus <= 12)){
             #          count = count + 1
             #          pdf(paste("clusters_",label[c_runs],"_",file_part,".pdf",sep=""),
             #            height=10.0,width=16.0)
             #          pdf(paste("clusters_",label[c_runs],"_",file_part,"_",count,".pdf",
             #            sep=""),height=10.0,width=16.0)
             par(mfrow=c(3,4))
             
             ### if number of clusters is > 12 plot always 20 clusters on one page
           } else if(i == 1 || ((i-1) %% 20 == 0)){
             #          count = count + 1
             #          pdf(paste("clusters_",label[c_runs],"_",file_part,"_",count,".pdf",
             #            sep=""),height=10.0,width=16.0)
             #          pdf(paste("clusters_",label[c_runs],"_",file_part,".pdf",sep=""),
             #            height=10.0,width=16.0)
             par(mfrow=c(4,5))
           }
           ### expression values of genes in cluster
           tmp <- dat_ztrans[kmeans_clus_final %in% i,]
           # if cluster contains only 1 gene
           if(length(kmeans_clus_final[kmeans_clus_final %in% i]) == 1){
             n_g = 1
           } else {   # if cluster contains > 1 genes
             n_g = nrow(tmp)
           }
           
           ### order gene expression data according to timepoints
           for(j in 1:n_g){
             # if cluster contains only 1 gene
             if(length(kmeans_clus_final[kmeans_clus_final %in% i]) == 1){
               #            tmp_plot = tmp[order(as.numeric(data_annotation[colnames(dat),
               #             "Time"]))]
               tmp_plot = tmp
             } else {       # if cluster contains > 1 genes
               #            tmp_plot = tmp[j,][order(as.numeric(data_annotation[colnames(dat),
               #            "Time"]))]
               tmp_plot = tmp[j,]
             }
             ### plot clusters
             if(j == 1){                     # if cluster contains only 1 gene
               plot(tmp_plot,type="l",xlab = "Time points", ylab="Z-score",
                    main=paste("Cluster ",i,sep=""),ylim=c(min(tmp),max(tmp)),
                    xaxt = "n")
               axis(1, at = 1:length(colnames(tmp)), labels = colnames(tmp), las = 3)
             } else {                        # if cluster contains > 1 genes
               points(colnames(tmp),tmp_plot,type="l")
             }
           }
           ### close PDFs depending on number of clusters
           #        if((i == n_clus) & (n_clus <= 12)){
           #          dev.off()
           #        } else if((i %% 20 == 0) || ((i == n_clus) & (i %% 20 !=0))){dev.off()}
         }
         #      if(length(dev.list()) > 0) { dev.off() }
         dev.off()
       }
       results[[c_runs*4 - 3]] <- kmeans_clus_final
       results[[c_runs*4 - 2]] <- cluster_means
       results[[c_runs*4 - 1]] <- n_pre_clus
       results[[c_runs*4]] <- n_fine_clusters
       names(results)[c((c_runs*4 - 3):(c_runs*4))] <- paste(c("kmeans_clus",
                                                               "cluster_means", "pre_clus","fine_clus"),label[c_runs],sep="_")
     }
     
     return(results)
}
