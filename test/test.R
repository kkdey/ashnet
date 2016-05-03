dim <- 1000
nsamples <- 50
pop_cov <- ashnet::generate_cov(dim);
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")

generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov$Sigma);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
plot(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l")
cor_sample <- cov2cor(cov_sample)
# image(nearPD(as.matrix(cor_sample), conv.tol = 1e-06)$mat)

shafer_mat <- ashnet::shafer_strimmer_shrinker(generate_sample);
shafer_eigen <- eigen(shafer_mat$mat, only.values = TRUE)

plot(sort(pop_cov$eigen, decreasing = TRUE), type="l", col="black")
lines(sort(as.numeric(shafer_eigen$values), decreasing = TRUE), type="l", col="blue")
lines(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l", col="green")

ash_cor_sample <- ashnet::ash_cor(cor_sample, nsamples);
#ash_cor_sample[abs(ash_cor_sample)  < 0.0001] <- 0
#image(nearPD(as.matrix(ash_cor_sample), conv.tol = 1e-06)$mat)

# samp_length <- c(10, 50, 200, 300, 500, 5);
# cols <- c("red", "green", "yellow", "blue", "cyan", "brown")

# plot(sort(as.numeric(pop_cov$eigen), decreasing = TRUE), type="l", col="black")

# for(num in 1:length(samp_length)){
#  ash_cor_sample <- ash_cor(cor_sample, samp_length[num]);
#  ash_cov_sample <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample%*%diag(sqrt(diag(cov_sample)))
#  eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)
#  lines(sort(as.numeric(eigen_ash_sample$values), decreasing = TRUE), type="l", col=cols[num])
#}


ash_cov_sample <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample%*%diag(sqrt(diag(cov_sample)))
eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)
plot(sort(as.numeric(eigen_ash_sample$values), decreasing = TRUE), type="l", col="red")
lines(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l", col="green")
lines(sort(as.numeric(shafer_eigen$values), decreasing = TRUE), type="l", col="blue")
lines(sort(as.numeric(pop_cov$eigen), decreasing = TRUE), type="l", col="black")

ash_cov_master_sample <- ashnet::ash_cov_master(generate_sample);
shrinkage <- ash_cov_master_sample$shrink_intensity
ash_cov_master <- ash_cov_master_sample$ash_cov_ledoit_wolf;
eigen_ash_cov_master <- eigen(ash_cov_master);


library(ggplot2)

eigendata <- data.frame(
    eigen_order = 1:dim,
    ash_cov_LW = sort(log(as.numeric(eigen_ash_cov_master$values)+1), decreasing=TRUE),
    ash_cov = sort(log(as.numeric(eigen_ash_sample$values)+1), decreasing = TRUE),
    sample_cov = sort(log(as.numeric(eigen_sample$values)+1), decreasing=TRUE),
    shafer_strimmer_cov= sort(log(as.numeric(shafer_eigen$values)+1), decreasing = TRUE),
    pop_cov = sort(log(as.numeric(pop_cov$eigen)+1), decreasing=TRUE))

colnames(eigendata) <- c("eigenorder",
                         "ash_cov_LW",
                         "ash_cov",
                         "sample_cov",
                         "shafer_strimmer",
                         "pop_cov")


# test_data <- data.frame(
#   var0 = 100 + c(0, cumsum(runif(49, -20, 20))),
#   var1 = 150 + c(0, cumsum(runif(49, -10, 10))),
#   date = seq.Date(as.Date("2002-01-01"), by="1 month", length.out=100))

library(ggplot2)
ggplot(eigendata, aes(eigenorder)) +
  geom_line(aes(y = ash_cov_LW, colour = "ash cov LW")) +
  geom_line(aes(y = ash_cov, colour = "ash cov"))+
  geom_line(aes(y = sample_cov, colour = "sample cov"))+
  geom_line(aes(y = shafer_strimmer, colour = "shafer cov"))+
  geom_line(aes(y = pop_cov, colour = "pop cov"))+
  xlab("Order of eigenvalues (sorted)")+
  ylab("log(Eigenvalues + 1)")+
  scale_colour_manual(values=c("blue", "red", "grey", "green", "orange"))+
  ggtitle(paste0("Eigenvalues distribution n/p=", round(nsamples/dim, 4),"\n for different shrinkage methods"))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))




# plot(sort(log(as.numeric(eigen_ash_cov_master$values)+1),
#           decreasing = TRUE),
#           type="l",
#           col="red",
#           ylim=c(0,10))
# lines(sort(log(as.numeric(pop_cov$eigen)+1), decreasing = TRUE),
#       type="l", col="black")
# lines(sort(log(as.numeric(eigen_sample$values)+1), decreasing = TRUE),
#       type="l", col="green")
# lines(sort(log(as.numeric(eigen_ash_sample$values)+1), decreasing = TRUE),
#       type="l", col="blue")
