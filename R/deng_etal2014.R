

################   Deng et al 2014  data analysis  ####################################

library(devtools)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()

deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

voom_deng_counts <- limma::voom(deng.counts)$E

dat <- t(voom_deng_counts)

pr <- prcomp(dat);

par(mfrow=c(1,1))
plot(pr$x[,1], pr$x[,2])

cor_sample <- cov2cor(cov(voom_deng_counts))
nsamples <- 10;

library(vegan)
screeplot(pr, bstick=TRUE, npcs=20)
bstick_var <- bstick(pr)
countPC <- numPC_generate(pr, npcs=20)

fields::image.plot(cor_sample)

nsamples <- 10

ash_cor_sample <- ash_cor(cor_sample, nsamples)$ash_cor_PD;
fields::image.plot(as.matrix(ash_cor_sample))

ash_cor_sample2 <- ash_cor2(cor_sample, nsamples)$ash_cor_PD;
fields::image.plot(as.matrix(ash_cor_sample2))

ash_cor_sample3 <- ash_cor3(cor_sample, nsamples)$ash_cor_PD;
fields::image.plot(as.matrix(ash_cor_sample3))

diags <- eigen(cor_sample)$values
ash_cov_sample3 <- diag(sqrt(diags))%*%ash_cor_sample3%*%diag(sqrt(diags))
ash_cov_sample2 <- diag(sqrt(diags))%*%ash_cor_sample2%*%diag(sqrt(diags))
ash_cov_sample <- diag(sqrt(diags))%*%ash_cor_sample%*%diag(sqrt(diags))
cov_sample <- diag(sqrt(diags))%*%cor_sample%*%diag(sqrt(diags))

eigen_ash_sample3 <- eigen(ash_cov_sample3, only.values = TRUE)
eigen_ash_sample2 <- eigen(ash_cov_sample2, only.values = TRUE)
eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)
eigen_sample <- eigen(cov_sample, only.values = TRUE)

library(ggplot2)

dim <- dim(cor_sample)[1]
dim <- 25
eigendata <- data.frame(
  eigenorder = 1:dim,
  sample_cov = sort(log(as.numeric(eigen_sample$values)+1),  decreasing=TRUE)[1:dim],
  ash_cov = sort(log(as.numeric(eigen_ash_sample$values)+1),  decreasing=TRUE)[1:dim],
  ash_cov2 = sort(log(as.numeric(eigen_ash_sample2$values)+1),  decreasing=TRUE)[1:dim],
  ash_cov3 = sort(log(as.numeric(eigen_ash_sample3$values)+1),  decreasing=TRUE)[1:dim])

colnames(eigendata) <- c("eigenorder",
                         "sample_cov",
                         "ash_cov",
                         "ash_cov2",
                         "ash_cov3")

library(ggplot2)
ggplot(eigendata, aes(eigenorder)) +
  geom_line(aes(y = sample_cov, colour = "sample cov")) +
  geom_line(aes(y = ash_cov, colour = "ash cov")) +
  geom_line(aes(y = ash_cov2, colour = "ash cov2"))+
  geom_line(aes(y = ash_cov3, colour = "ash cov 3"))+
  xlab("Order of eigenvalues (sorted)")+
  ylab("log(Eigenvalues + 1)")+
  scale_colour_manual(values=c("blue", "red", "green", "orange"))+
  ggtitle(paste0("Eigenvalues distribution n/p=", round(dim(cor_sample)[1]/nsamples, 4),"\n for different shrinkage methods"))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


#################  image plots for the correlation matrices  ##########################

library(fields)
set.seed(1)
par(mfrow=c(2,2))
cols = gray.colors(100)
image.plot(cov2cor(pop_cov$Sigma), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample2), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample3), col=cols, nlevel=100)

###############   PCA  plot  #######################################

eigen_ash <- eigen(ash_cov_sample)
plot(pr$x[,1], pr$x[,2])

######################   eigenvalues vs prcomp  ############################

pr <- prcomp(dat)
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");

plot(pr$x[,1], pr$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

mean_dat <- (rep(1,dim(dat)[1]))%*% t(colMeans(dat))
dat_adjusted <- dat - mean_dat
pr1 <- eigen(dat_adjusted %*% t(dat_adjusted))
eigs <- eigen(cov_sample)$values

library(expm)
projected_dat <- sqrtm(ash_cov_sample)%*% sqrtm(solve(cov_sample)) %*% dat
pr_projected <- prcomp(projected_dat)
plot(pr_projected$x[,1], pr_projected$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

projected_dat <- sqrtm(ash_cov_sample2)%*% sqrtm(solve(cov_sample)) %*% dat
pr_projected <- prcomp(projected_dat)
plot(pr_projected$x[,1], pr_projected$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

projected_dat <- sqrtm(ash_cov_sample3)%*% sqrtm(solve(cov_sample)) %*% dat
pr_projected <- prcomp(projected_dat)
plot(pr_projected$x[,1], pr_projected$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))
