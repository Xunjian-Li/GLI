################
confident_interval_fun = function(beta_mat){
  upper_beta = 0
  lower_beta = 0
  p = dim(beta_mat)[2]
  for(i in c(1:p)){
    upper_beta[i] = quantile(beta_mat[,i],0.975)
    lower_beta[i] = quantile(beta_mat[,i],0.025)
  }
  beta_co = rbind(lower_beta,upper_beta)
}

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
  
  est_par_pbb_BFGS = function(W,y){
    
    start_time = Sys.time()
    loglikelihood.fun = function(beta, sigma, al){
      a.1 = 1.1396; ka = 0.618
      la = (a.1/2*(2+al)/(a.1-al))^ka
      g.1.log = digamma(1) - digamma(la)
      g.1 = exp(g.1.log)
      g.2 = sqrt(trigamma(la) + trigamma(1))
      ai = (y-W%*%beta)*g.2/sigma - g.1.log
      -sum(ai + (la+1) * log(1+exp(-ai)) - log(la*g.2/sigma))
    }
    
    Partial = function(para){
      beta  = para[1:p]
      sigma = para[p+1]
      al    = para[p+2]
      n_row = length(y)
      a.1 = 1.1396; ka = 0.618
      la = (a.1/2*(2+al)/(a.1-al))^ka
      g.1.log = digamma(1) - digamma(la)
      g.1 = exp(g.1.log)
      g.2 = sqrt(trigamma(la) + trigamma(1))
      residual = y-W%*%beta
      ai = residual*(g.2/sigma) - g.1.log
      
      fun.al = function(al){
        la = (a.1/2*(2+al)/(a.1-al))^ka
        g.1.log = digamma(1) - digamma(la)
        g.1 = exp(g.1.log)
        g.2 = sqrt(trigamma(la) + trigamma(1))
        ai = residual*(g.2/sigma) - g.1.log
        s1 = (la-exp(ai))/(1+exp(ai))
        s2 = residual/(2*sigma*g.2)*psigamma(la, deriv = 2) + trigamma(la)
        s3 = -log(1+exp(-ai)) + 1/la + psigamma(la, deriv = 2)/(2*g.2^2)
        s4 = 1.37056*(al + 2)^(-0.382)/(1.1396 - al)^1.618
        sum((s1*s2+s3)*s4)
      }
      part1 = t(g.2/sigma*(1-(la+1)/(1+exp(ai))))%*%W
      part2 = -n_row/sigma + t(1-(la+1)/(1+exp(ai)))%*%residual*g.2/sigma^2
      part3 = fun.al(al)
      
      c(part1,part2,part3)
    }
    
    n = length(y); p = dim(W)[2]; error.min = 1e-10; esp.all = 1
    esp = 1; al=0; beta =rep(0,p); sigma = 1; k =0; al =0.1;loglikelihood = 0
    li.all = c(); 
    li.all.old = loglikelihood.fun(beta, sigma, al)
    I0 = diag(1,p+2)
    Gt = -I0
    
    i = 1
    para = c(beta, sigma, al)
    if(i==1){
      gradient = Partial(para)
      ## setting step size
      step.len = 0.01/max(abs(gradient))
      does1 = 1
      while(does1){
        ## first updating
        para.new = para - step.len*(-I0)%*%gradient
        li.all.new = loglikelihood.fun(beta, sigma, al)
        if(is.null(li.all.new) ||  li.all.new<li.all.old){
          step.len = 0.8*step.len
        }else{
          does1 = 0
        }
      }
      li.all.old = li.all.new
      li.all[i] = li.all.old
    }
    
    while(esp.all>error.min  && i<1e4){
      ## calculating gradient
      gradient.new = Partial(para.new)
      ## calculating Gt
      yt = as.matrix(gradient.new - gradient)
      st = matrix(c(para.new - para),p+2,1,byrow=F)
      Gt = (I0-st%*%t(yt)/as.numeric(t(st)%*%yt))%*%Gt%*%(I0-yt%*%t(st)/as.numeric(t(st)%*%yt)) + st%*%t(st)/as.numeric(t(st)%*%yt)
      para = para.new
      step.len = 1
      para.new.tra = para.new - step.len*Gt%*%gradient
      while(para.new.tra[p+1]<=0 || para.new.tra[p+2] > 1.1396 || para.new.tra[p+2] < -2){
        step.len = 0.5*step.len
        para.new.tra = para.new - step.len*Gt%*%gradient
      }
      li.all.new = loglikelihood.fun(para.new.tra[c(1:p)],para.new.tra[p+1],para.new.tra[p+2])
      
      while(is.null(li.all.new)  || li.all.new< li.all.old){
        step.len = 0.5*step.len
        para.new.tra = para.new - step.len*Gt%*%gradient
        li.all.new = loglikelihood.fun(para.new.tra[c(1:p)],para.new.tra[p+1],para.new.tra[p+2])
      }
      para.new = para.new.tra
      
      if(i>1){ esp.all = abs(li.all.new-li.all.old)/(abs(li.all.new)+1) }
      li.all.old = li.all.new
      gradient = gradient.new
      i = i+1
      li.all[i] = li.all.new
    }
    
    end_time = Sys.time()
    times = as.numeric(difftime(end_time, start_time, units = "secs"))
    p.dimen = p+2
    
    list(n, beta = para.new[c(1:p)], sigma = para.new[p+1], al = para.new[p+2], iter = i, times = times,
         logL = li.all.new, AIC = 2*p.dimen-2*li.all.new)
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
  beta_ori = Beta_ini
  beta0 = Beta_ini[c(1:p)];
  sigma0 = Beta_ini[p+1];
  al0 = Beta_ini[p+2];
  
  y = X.sample(beta0, sigma0, al0, W)
  est_beta1 = est_mean_reg(W,y)
  est_beta2 = est_par_pbb_BFGS(W,y)
  results0 = c(est_beta1,est_beta2)
  
  return(results0)
  
}

library(MASS)
library(parallel)

sample_size_mini = 10000

p.all = c(4,10,20,40,60,80,100,150,200,300,400,500,600)
sample_number = length(p.all)
set.seed(100)

start_time = Sys.time()
mark11 = 1
for(zzz.count in c(3:3)){

  for(kkk.count in c(1:sample_number)){
    
    # if(zzz.count==1){p = p.all[kkk.count]; Beta_ini = c(rep(c(-1,2),p/2), 1, -0.5)} # bbe, sigma, alpha
    # if(zzz.count==2){p = p.all[kkk.count]; Beta_ini = c(rep(c(-1,2),p/2), 1, 0)}
    if(zzz.count==3){p = p.all[kkk.count]; Beta_ini = c(rep(c(-1,2),p/2), 1, 0.5)}
    
    # 设置矩阵维度和衰减系数
    alpha <- 0.1  # 调整衰减速率
    cov_matrix <- outer(1:(p-1), 1:(p-1), function(i, j) exp(-alpha * abs(i - j)))

    mean_vector <- rep(0,(p-1))
    W <- mvrnorm(n = sample_size_mini, mu = mean_vector, Sigma = cov_matrix)
    W = cbind(matrix(1,sample_size_mini,1),W)
    W = as.matrix(W)

    times = 100

    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    clusterEvalQ(cl, library(MASS))
    clusterExport(cl, "W");
    clusterExport(cl, "Beta_ini");
    res0 <- parLapply(cl, 1:times, est_beta_fun)
    stopCluster(cl)
    res.all<- matrix(unlist(res0),times,(p+7)*2,byrow=T)
    
    ## estimation
    res1 = res.all[,1:(p+7)]
    R_all.mse = res1[,2:(p+3)]
    bis = apply(R_all.mse,2,mean) - Beta_ini
    mse = (apply(R_all.mse,2,mean)-Beta_ini)^2 + apply(R_all.mse,2,var)
    avg_iter = mean(res1[,(p+4)])
    avg_time = mean(res1[,(p+5)])
    avg_LogL = mean(res1[,(p+6)])
    est_data = rep(0, 4)
    est_data[1] = mean(mse)
    est_data[2] = avg_iter
    est_data[3] = avg_time
    est_data[4] = avg_LogL
    A1.1 = c(est_data, res1[,(p+4)], res1[,(p+5)], res1[,(p+6)])
    
    ## estimation
    res1 = res.all[,(p+8):((p+7)*2)]
    R_all.mse = res1[,2:(p+3)]
    bis = apply(R_all.mse,2,mean) - Beta_ini
    mse = (apply(R_all.mse,2,mean)-Beta_ini)^2 + apply(R_all.mse,2,var)
    
    avg_iter = mean(mean(res1[,(p+4)]))
    avg_time = mean(mean(res1[,(p+5)]))
    avg_LogL = mean(mean(res1[,(p+6)]))
    est_data = rep(0, 4)
    est_data[1] = mean(mse)
    est_data[2] = avg_iter
    est_data[3] = avg_time
    est_data[4] = avg_LogL
    A1.2 = c(est_data, res1[,(p+4)], res1[,(p+5)], res1[,(p+6)])

    if(mark11==1){B1 = c(p, A1.1,A1.2)}else{B1 = rbind(B1,c(p, A1.1,A1.2))}
    mark11 = mark11+1
  }
}
end_time = Sys.time()
print(end_time-start_time)

# computational_costs
B1[,c(c(1:5),c((6+3*times):(9+3*times)))]


