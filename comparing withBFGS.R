### convergence of BFGS
library(MASS)

X.sample = function(beta, sigma, al, W){
  a.1 = 1.1396; ka = 0.618
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
    flag = flag[ai.neg>400] # 如果ai.neg很大，避免求 exp(ai.neg)
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
    
    ## update beta
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
    
    ## update sigma
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
  times = as.numeric(difftime(end_time, start_time, units = "secs"))
  aic0 = 2*(p+2)-2*loglikelihood[k]
  list(beta = t(beta), sigma = sigma, alpha = al, iter = k, times = times, logL = loglikelihood[k], AIC = aic0)
}

est_par_pbb_BFGS = function(W,y){
  
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
    return(c(part1,part2,part3))
  }
  
  start_time = Sys.time()
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
  
  p.dimen = p+2
  list.all = list()
  list.all$beta = para.new[c(1:p)]
  list.all$sigma = para.new[p+1]
  list.all$al = para.new[p+2]
  list.all$times = as.numeric(difftime(end_time, start_time, units = "secs"))
  list.all$loglikelihood=li.all.new
  list.all$AIC = 2*p.dimen-2*li.all.new
  list.all$BIC = log(n)*p.dimen-2*li.all.new
  list.all$iter = i-1
  return(list.all)
}

set.seed(100)
n.times = 1000
n1 = 400
# p = 2; Beta_ini = c(-1,2)
p = 4; Beta_ini = c(-1,2,1,-2)
para.all = matrix(0,n.times,2*(p+4))
for(kk in c(1:n.times)){
  W = as.matrix(cbind(matrix(1,n1,1),matrix(rnorm(n1*(p-1),0.4,sqrt(0.4)),n1,(p-1))))
  ## data process
  beta = Beta_ini; sigma = 1; al = -0.5
  y = X.sample(beta,sigma,al,W)
  eat11.BFGS = est_par_pbb_BFGS(W,y)
  eat11.QLB = est_mean_reg(W,y)
  para.all[kk,] = c(eat11.BFGS$beta,eat11.BFGS$sigma,eat11.BFGS$al,eat11.BFGS$loglikelihood,eat11.BFGS$iter,
               eat11.QLB$beta,eat11.QLB$sigma,eat11.QLB$alpha,eat11.QLB$logL,eat11.QLB$iter)
}

colMeans(abs(para.all[,c(1:(p+2))] - matrix(c(beta,sigma,al), n.times, p+2,byrow = T)))
colMeans(abs(para.all[,c((p+5):(2*p+6))] - matrix(c(beta,sigma,al), n.times, p+2,byrow = T)))

# 设置 PDF 输出文件路径和长宽比
# pdf("*./difference.pdf", width = 10, height = 5)
res1 = para.all[,2*p+7]-para.all[,p+3]
par(mfrow=c(1,2))
hist(res1, main = '(a1)', xlab = 'Difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
boxplot(res1, main = '(a2)', ylab = 'Difference of loglikelihoods',cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
# dev.off()
