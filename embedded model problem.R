n.times = 100


res.all = matrix(0,n.times,7)
for(zz in c(1:n.times)){
  n1 = 100
  p = 2
  Beta_ini = c(-1,2)
  W = cbind(matrix(1,n1,1),matrix(rnorm(n1*(p-1),0.4,sqrt(0.4)),n1,(p-1)))
  W = as.matrix(W)
  
  X.sample = function(sigma,al,beta,W){
    a.1 = 1.1396; ka = 0.618
    la = (a.1/2*(2+al)/(a.1-al))^ka
    g.1.log = digamma(1) - digamma(la)
    g.1 = exp(g.1.log)
    g.2 = sqrt(trigamma(la) + trigamma(1))

    n = dim(W)[1]; mu = W%*%beta; y = runif(n);
    y1.tra =y^(-1/la)

    x.sam = mu - log((y1.tra-1)/g.1)*sigma/g.2
  }
  
  # X.sample = function(sigma,al,beta,W){
  #   n = dim(W)[1]; mu = W%*%beta; y = runif(n);
  #   y1.tra =y^(-1/al)
  #   x.sam = rep(0,n)
  #   x.sam[is.infinite(y1.tra)] = 1/al*log(y[is.infinite(y1.tra)])*sigma + mu[is.infinite(y1.tra)]
  #   x.sam[!is.infinite(y1.tra)] = -log(y1.tra[!is.infinite(y1.tra)]-1)*sigma + mu[!is.infinite(y1.tra)]
  #   x.sam
  # }
  
  ## data process
  beta = Beta_ini
  al = 1.13
  # al = 100
  sigma = 1
  
  y = X.sample(sigma,al,beta,W)
  
  
  loglikelihood.fun = function(sigma,al,beta,W,y){
    a.1 = 1.1396; ka = 0.618
    la = (a.1/2*(2+al)/(a.1-al))^ka
    g.1.log = digamma(1) - digamma(la)
    g.1 = exp(g.1.log)
    g.2 = sqrt(trigamma(la) + trigamma(1))
    ai = (y-W%*%beta)*g.2/sigma - g.1.log
    
    -sum(ai + (la+1) * log(1+exp(-ai)) - log(la*g.2/sigma))
  }
  
  W.GL = W
  X.GL = y
  
  reu1 = est_mean_reg(W.GL,X.GL)
  reu2 = est_beta_for_GL1(W.GL,X.GL)
  
  res.all[zz,] = c(reu1$logL, reu1$iter, reu1$times, reu2$logL , reu2$iter, reu2$times, reu2$alpha)
  
}

res.all


# res.all.old1 = res.all
# res.all.old = res.all
res1 = res.all[,1]-res.all[,4]
sum(res.all[,7]>=100)
par(mfrow=c(2,2))
hist(res1,main = '(a1)',xlab = 'difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
boxplot(res1,main = '(a2)',ylab = 'difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)


colMeans(res.all)
colMeans(res.all.old)



# # res.all.old1 = res.all
# # res.all.old = res.all
res1 = res.all.old[,1]-res.all.old[,4]
sum(res.all.old[,7]>=200)
par(mfrow=c(2,2))
hist(res1,main = '(a1)',xlab = 'difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
boxplot(res1,main = '(a2)',ylab = 'difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)

res2 = res.all.old1[,1]-res.all.old1[,4]
sum(res.all.old1[,7]>=200)
hist(res2,main = '(a3)',xlab = 'difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
boxplot(res2,main = '(a4)',ylab = 'difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)

# colMeans(res.all.old1)[1]-colMeans(res.all.old1)[4]

# write.csv(res.all.old,"/Users/caleblee/Desktop/skill summit time/2022年 正在投稿的文章/GL1_reg/programs/updated/embedded mean.csv")
# write.csv(res.all.old1,"/Users/caleblee/Desktop/skill summit time/2022年 正在投稿的文章/GL1_reg/programs/updated/embedded loc.csv")


write.csv(res.all.old,"/Users/caleblee/Desktop/skill summit time/2022年 正在投稿的文章/GL1_reg/programs/updated/embedded mean1.csv")
write.csv(res.all.old1,"/Users/caleblee/Desktop/skill summit time/2022年 正在投稿的文章/GL1_reg/programs/updated/embedded loc1.csv")

