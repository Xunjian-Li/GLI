est_beta_fun = function(n.repeat){
  
  a.1 = 1.1396; ka = 0.618
  a_1 = -2; a_2 = 1.139547
  
  X.sample = function(beta,sigma,al,W){
    la = (a.1/2*(2+al)/(a.1-al))^ka
    g.1.log = digamma(1) - digamma(la)
    g.1 = exp(g.1.log)
    g.2 = sqrt(trigamma(la) + trigamma(1))
    n = dim(W)[1]; mu = W%*%beta; y = runif(n);
    y1.tra =y^(-1/la)
    x.sam = mu - log((y1.tra-1)/g.1)*sigma/g.2
  }
  
  est_mean_reg = function(W,y){
    
    est_al_bisection = function(al, residual){
      kk.bisec = 0; al.all = 0; esp0 = 1; weight1 = 0.45
      f_al = fun.al(al, residual)
      while(abs(f_al) > 1e-10 & kk.bisec <1e4 & esp0>1e-10){
        if(f_al >= 0){
          al.new = weight1*al + (1-weight1)*a_2
          a_1 = al
        }else{
          al.new = weight1*al + (1-weight1)*a_1
          a_2 = al
        }
        f_al = fun.al(al.new, residual)
        esp0 = min(abs(al.new - a_2),abs(al.new - a_1))
        al = al.new
        kk.bisec = kk.bisec+1
        al.all[kk.bisec] = al.new
      }
      return(al)
    }
    
    fun.al = function(al, residual){
      flag = c(1:length(y))
      la = (a.1/2*(2+al)/(a.1-al))^ka
      g.1.log = digamma(1) - digamma(la)
      g.1 = exp(g.1.log)
      g.2 = sqrt(trigamma(la) + trigamma(1))
      ai.neg = - residual*(g.2/sigma) + g.1.log
      
      flag = flag[ai.neg>400]
      
      exp_ai <- exp(ai.neg)
      t11.1 = la*exp_ai-1
      t11.2 = exp_ai+1
      s1 = t11.1/t11.2
      t11.2.log = log(t11.2)
      
      s1[flag] = la
      t11.2.log[flag] = ai.neg
      
      t11.3 = psigamma(la, deriv = 2)/(2*sigma*g.2)
      s2 = residual*t11.3 + trigamma(la)
      s3 = - t11.2.log + 1/la + psigamma(la, deriv = 2)/(2*g.2^2)
      s4 = ka*(a.1+2)*(a.1/2)^ka*(al + 2)^(-0.382)/(a.1 - al)^1.618
      return(sum((s1*s2+s3)*s4))
    }
    
    W.GL = W
    X.GL = y
    
    start_time = Sys.time()
    esp = 1; al=0; beta =rep(0,p); sigma = 1; k =0; al =0.1; loglikelihood = 0; mark = 1
    B_mat = ginv(t(W.GL)%*%W.GL)
    n = dim(W.GL)[1]
    la = (a.1/2*(2+al)/(a.1-al))^ka
    g.1.log = digamma(1) - digamma(la)
    g.1 = exp(g.1.log)
    g.2 = sqrt(trigamma(la) + trigamma(1))
    residual <- y - W %*% beta
    ai = residual*(g.2/sigma) - g.1.log
    while(esp>1e-10 & k<10000){
      k = k+1
      exp_ai <- exp(ai)
      denominator <- 1 + exp_ai
      factor1 <- (exp_ai - la) / (denominator * g.2)  # 提取出第一部分的计算
      factor2 <- (la + 1) / (4 * sigma)            # 提取出第二部分的计算
      W_GL_beta <- W.GL %*% beta
      factor_combined <- factor1 + factor2 * W_GL_beta
      A <- crossprod(factor_combined, W.GL) / sigma
      beta.new = (4*sigma^2)/(la+1) *B_mat%*% t(A)
      
      ## update ai
      beta = beta.new
      residual <- y - W %*% beta
      ai = residual*g.2/sigma - g.1.log
      exp_ai <- exp(ai)
      denominator <- 1 + exp_ai
      term1 <- ((exp_ai - la) / denominator) / g.2
      term2 <- ((la + 1) * residual) / (4 * sigma)
      B <- -g.2^2 * mean((term1 - term2) * residual)
      C = -(la+1)*g.2^2/4*norm(residual, type = "2")^2 / length(residual)
      sigma.new = (-B+sqrt(B^2-4*C))/2
      
      ## update ai
      sigma = sigma.new
      ai = residual * g.2/sigma - g.1.log
      
      ## updata al
      al.new = est_al_bisection(al, residual)
      al = al.new
      la = (a.1/2*(2+al)/(a.1-al))^ka
      g.1.log = digamma(1) - digamma(la)
      g.1 = exp(g.1.log)
      g.2 = sqrt(trigamma(la) + trigamma(1))
      ai = residual*g.2/sigma - g.1.log
      loglikelihood[k] = - sum(ai + (la+1) * log(1+exp(-ai)) - log(la*g.2/sigma))
      
      if(k==1){ esp = 1 }else{ esp = abs(loglikelihood[k]-loglikelihood[k-1])/abs(loglikelihood[k]) }
    }
    end_time = Sys.time()
    aic0 = 2*(p+2)-2*loglikelihood[k]
    list(n, beta = t(beta), sigma = sigma, alpha = al, iter = k, times = as.numeric(difftime(end_time, start_time, units = "secs")),
         logL = loglikelihood[k], AIC = aic0)
  }
  
  p = dim(W)[2]
  n_row = dim(W)[1]
  
  beta0 = Beta_ini[c(1:p)];
  sigma0 = Beta_ini[p+1];
  al0 = Beta_ini[p+2];
  
  y = X.sample(beta0,sigma0,al0,W)
  est_beta1 = est_mean_reg(W,y)
  
  results0 = est_beta1
}

library(MASS)
library(parallel)

sample_number = 50
sample_size_mini = 200
start_time = Sys.time()
mark11 = 1
for(zzz.count in c(1:3)){
  if(zzz.count==1){p = 20; Beta_ini = c(rep(c(-1,2), p/2), 1, -0.5)} # bbe, sigma, alpha
  if(zzz.count==2){p = 20; Beta_ini = c(rep(c(-1,2),p/2), 1, 0)}
  if(zzz.count==3){p = 20; Beta_ini = c(rep(c(-1,2),p/2), 1, 0.5)}
  
  for(kkk.count in c(1:sample_number)){
    set.seed(100)
    # 设置矩阵维度和衰减系数
    alpha <- 0.1  # 调整衰减速率
    # 创建一个 100x100 的矩阵，元素通过指数衰减定义
    cov_matrix <- outer(1:(p-1), 1:(p-1), function(i, j) exp(-alpha * abs(i - j)))
    
    mean_vector <- rep(0,(p-1))
    W <- mvrnorm(n = kkk.count*sample_size_mini, mu = mean_vector, Sigma = cov_matrix)
    W = cbind(matrix(1,kkk.count*sample_size_mini,1),W)
    W = as.matrix(W)
    
    times = 10000
    
    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    clusterEvalQ(cl, library(MASS))
    clusterExport(cl, "W");
    clusterExport(cl, "Beta_ini");
    res0 <- parLapply(cl, 1:times, est_beta_fun)
    stopCluster(cl)
    res1<- matrix(unlist(res0),times,(p+7),byrow=T)
    
    ## estimation
    R_all.mse = res1[,2:(p+3)]
    bis = apply(R_all.mse,2,mean) - Beta_ini
    mse = (apply(R_all.mse,2,mean)-Beta_ini)^2 + apply(R_all.mse,2,var)
    avg_iter = mean(mean(res1[,(p+4)]))
    avg_time = mean(mean(res1[,(p+5)]))
    avg_LogL = mean(mean(res1[,(p+6)]))
    
    
    est_data = rep(0, (p+2)*2+6)
    est_data[1] = kkk.count*sample_size_mini
    est_data[c(2:(p+3))] = bis
    est_data[c((p+4):(2*p+5))] = mse
    est_data[2*p+6] = sum(abs(bis))
    est_data[2*p+7] = sum(mse)
    est_data[2*p+8] = avg_iter
    est_data[2*p+9] = avg_time
    est_data[2*p+10] = avg_LogL
    
    if(mark11==1){A1 = est_data}else{A1 = rbind(A1,est_data)}
    mark11 = mark11+1
  }
}
end_time = Sys.time()
print(end_time-start_time)

if(zzz.count==1){p = 20; Beta_ini = c(rep(c(-1,2),p/2),-0.5,1)}
if(zzz.count==2){p = 20; Beta_ini = c(rep(c(-1,2),p/2),0,1)}
if(zzz.count==3){p = 20; Beta_ini = c(rep(c(-1,2),p/2),0.5,1)}

p = 20
par(mfrow=c(3,2))
for(i in c(1:3)){
  ranks111 = c(1:sample_number) + sample_number*(i-1)
  # par(mfrow=c(1,1))
  # 使用 matplot 绘制多条线
  matplot(c(1:sample_number)*200,A1[ranks111,c(2:(p+3))], ylim = c(-0.008,0.008), type="l",lwd = 2, lty=1:(p+2), col=1:(p+2), xlab="Sample size", ylab="Bias", main="Bias of parameters",cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
  legend("topright", inset = c(0.05, 0.00), bty='n', col=1:(p+2), ncol=5,  lty=1:(p+2),lwd = 2,  legend=c(as.expression(lapply(1:p, function(i) bquote(beta[.(i)]))), 
                                                                                                           expression(alpha, delta)), cex=1)
  lines(c(0,10000),c(0,0))
  # 使用 matplot 绘制多条线
  
  matplot(c(1:sample_number)*200,log(A1[ranks111,c((p+4):(2*p+5))]), type="l",lwd = 2, lty=1:(p+2), col=1:(p+2), xlab="Sample size", ylab="log(MSE)", main="Logarithm of parameters' MSEs",cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
  legend("topright", inset = c(0.05, 0.00), bty='n', col=1:(p+2), ncol=5,  lty=1:(p+2),lwd = 2,  legend=c(as.expression(lapply(1:p, function(i) bquote(beta[.(i)]))), 
                                                                                                           expression(alpha, delta)), cex=1)
  # lines(c(0,10000),c(0,0))
}

par(mfcol=c(4,3))
for(i in c(1:3)){
  ranks111 = c(1:sample_number) + sample_number*(i-1)
  plot(c(1:sample_number)*200, log(A1[ranks111,46]),type='l', ylab="log(|Biases|)", xlab = "Sample size",main="Logarithm of absolute biases' sum",cex.lab=1.4,cex.axis=1.2,cex.main=1.5)
  print(A1[ranks111,46])
  plot(c(1:sample_number)*200, log(A1[ranks111,47]),type='l',main= "Logarithm of MSEs' sum", ylab="log(MSEs)", xlab = "Sample size",cex.lab=1.4,cex.axis=1.2,cex.main=1.5)
  plot(c(1:sample_number)*200, A1[ranks111,48],type='l',main="Average number of iterations", ylab = "Iterations", xlab = "Sample size",cex.lab=1.4,cex.axis=1.2,cex.main=1.5)
  plot(c(1:sample_number)*200, A1[ranks111,49],type='l',main="Average runtime", ylab = "Time (s)", xlab = "Sample size",cex.lab=1.4,cex.axis=1.2,cex.main=1.5)
  # plot(A1[ranks111,50],type='l',main="The average of logLiklihood", ylab = "LogLiklihood", xlab = "Sample size",cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
  
}
