# y = x
library(MASS)

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



# 设置随机种子以便复现
# set.seed(123)
sample_size = 20000
p  = 10
# 设置矩阵维度和衰减系数
alpha <- 0.1  # 调整衰减速率
# 创建一个 100x100 的矩阵，元素通过指数衰减定义
cov_matrix <- outer(1:(p-1), 1:(p-1), function(i, j) exp(-alpha * abs(i - j)))

# 生成多元正态分布的数据（例如 100 个样本）
mean_vector <- rep(0,(p-1))
W <- mvrnorm(n = sample_size, mu = mean_vector, Sigma = cov_matrix)
W = as.matrix(cbind(matrix(1,sample_size,1),W))
## data process
beta = rep(c(-1,2),p/2); al = 0.5; sigma = 1
y = X.sample(beta, sigma, al, W)
reu1 = est_mean_reg(W,y)
reu1