# X: covariates
# y: estimated tau = estimated individual treatment effect
# nt: # of trees
# S: minumum size of terminal nodes - larger number causes smaller trees to be grown
# L: maximum number of terminal nodes
genrulesRF = function(X, y, nt,S ,L){
  N = dim(X)[1]
  sf = min(1, (11*sqrt(N)+1)/N)
  mn = 2+floor(rexp(1, 1/(L-2)))
  ns = S
  forest = randomForest(x = X, y=y, sampsize = sf*N ,replace=F, ntree =1, maxnodes=mn, nodesize = ns)
  for(i in 2:nt) {
    mn = 2+floor(rexp(1, 1/(L-2)))
    ns = S
    model1 = randomForest(x = X, y=y, sampsize = sf*N ,replace=F, ntree =1, maxnodes=mn, nodesize = ns)
    forest = combine(forest, model1)
  }
  treelist = RF2List(forest)
  rules = extractRules(treeList=treelist, X=X, ntree=nt, maxdepth=15)
  rules = c(rules)
  #rules = rules[take1(length(rules))]
  rulesmat = matrix(rules)
  colnames(rulesmat) = "condition"
  metric = getRuleMetric(rulesmat,X,y)
  pruned = pruneRule(metric, X, y, 0.025, typeDecay = 2)
  unique(pruned[,4])
}

genrulesRF.take1 = function(X, y, nt,S ,L){
  N = dim(X)[1]
  sf = min(1, (11*sqrt(N)+1)/N)
  mn = 2+floor(rexp(1, 1/(L-2)))
  ns = S
  forest = randomForest(x = X, y=y, sampsize = sf*N ,replace=F, ntree =1, maxnodes=mn, nodesize = ns)
  for(i in 2:nt) {
    mn = 2+floor(rexp(1, 1/(L-2)))
    ns = S
    model1 = randomForest(x = X, y=y, sampsize = sf*N ,replace=F, ntree =1, maxnodes=mn, nodesize = ns)
    forest = combine(forest, model1)
  }
  treelist = RF2List(forest)
  rules = extractRules(treeList=treelist, X=X, ntree=nt, maxdepth=15)
  rules = c(rules)
  rules = rules[take1(length(rules))]
  rulesmat = matrix(rules)
  colnames(rulesmat) = "condition"
  metric = getRuleMetric(rulesmat,X,y)
  pruned = pruneRule(metric, X, y, 0.025, typeDecay = 1)
  unique(pruned[,4])
}
take1 = function(len) {
  out = c()
  i = 0
  while (i < len){
    out = c(out, i+sample(1:2))
    i = i+2
  }
  out = out[1:len]
  out[seq(1, len, 2)]
}

createX = function(X, rules, t, corelim=1){
  Xr = matrix(0, nrow=dim(X)[1], ncol=length(rules))
  for (i in 1:length(rules)){
    Xr[eval(parse(text = rules[i])),i] = 1
  }
  
  Nr = dim(Xr)[2]
  ind = 1:Nr
  if(dim(X)[1]<200){
    t= 0.05
  }
  sup = apply(Xr, 2, mean)
  elim = which((sup<t)|(sup>(1-t)))
  
  if(length(elim)>0){
    ind = ind[-elim]
  }

  C = cor(Xr[,ind])

  diag(C) = 0
  Nr = dim(Xr[,ind])[2]
  elim=c()
  for(i in 1:(Nr-1)){
    #elim = c(elim, which(round(abs(C[i,(i+1):Nr]), digits=4)>=corelim) +i)
    elim = c(elim, which(round((C[i,(i+1):Nr]), digits=4)>=corelim) +i)
  }
  
  if(length(elim)>0){
    ind = ind[-elim]
  }else{
    ind = ind
  }
  #ind = ind[-elim]
  Xr = Xr[,ind]
  rules = rules[ind]
  list(data.matrix(Xr), rules)
}
createX.disjoint = function(X, rules, t, corelim=1){
  Xr = matrix(0, nrow=dim(X)[1], ncol=length(rules))
  for (i in 1:length(rules)){
    Xr[eval(parse(text = rules[i])),i] = 1
  }
  
  Nr = dim(Xr)[2]
  ind = 1:Nr
  if(dim(X)[1]<200){
    t= 0.05
  }
  sup = apply(Xr, 2, mean)
  elim = which((sup<t)|(sup>(1-t)))
  
  if(length(elim)>0){
    ind = ind[-elim]
  }
  
  C = cor(Xr[,ind])
  
  diag(C) = 0
  Nr = dim(Xr[,ind])[2]
  elim=c()
  for(i in 1:(Nr-1)){
    #elim = c(elim, which(round(abs(C[i,(i+1):Nr]), digits=4)>=corelim) +i)
    elim = c(elim, which(round((C[i,(i+1):Nr]), digits=4)>=corelim) +i)
  }
  
  if(length(elim)>0){
    ind = ind[-elim]
  }else{
    ind = ind
  }
  #ind = ind[-elim]
  
  ## remove coarser rules
  Xr = Xr[,ind]
  rules = rules[ind]
  Nr = dim(Xr)[2]
  Nunit = dim(Xr)[1]
  
  exc.vec = rep(0, Nr)
  for(i in 1:(Nr-1)){
    for(j in (i+1):Nr){
      temp.vec1 = Xr[,i]
      temp.vec2 = Xr[,j]
      vec.diff = temp.vec1 - temp.vec2
      num.diff.zero = sum(vec.diff == 0)
      num.diff.posi.one = sum(vec.diff == 1)
      num.diff.nega.one = sum(vec.diff == -1)
      ## if num.diff.zero = Nunit, rule i and rule j are identical
      ## if num.diff.posi.one = 0, rule j includes rule i => remove rule j
      ## if num.diff.nega.one = 0, rule i includes rule j => remove rule i
      
      if(num.diff.posi.one == 0){
        exc.vec[j] = 1
      }else if(num.diff.nega.one == 0){
        exc.vec[i] = 1
      }
    }
  }
  
  Xr = Xr[,exc.vec==0]
  rules = rules[exc.vec==0]
  list(data.matrix(Xr), rules)
}
createX.take1 = function(X, rules, t, corelim=1){
  Xr = matrix(0, nrow=dim(X)[1], ncol=length(rules))
  for (i in 1:length(rules)){
    Xr[eval(parse(text = rules[i])),i] = 1
  }
  
  Nr = dim(Xr)[2]
  ind = 1:Nr
  if(dim(X)[1]<200){
    t= 0.05
  }
  sup = apply(Xr, 2, mean)
  elim = which((sup<t)|(sup>(1-t)))
  
  if(length(elim)>0){
    ind = ind[-elim]
  }
  
  C = cor(Xr[,ind])
  
  diag(C) = 0
  Nr = dim(Xr[,ind])[2]
  elim=c()
  for(i in 1:(Nr-1)){
    elim = c(elim, which(round(abs(C[i,(i+1):Nr]), digits=4)>=corelim) +i)
    #elim = c(elim, which(round((C[i,(i+1):Nr]), digits=4)>=corelim) +i)
  }
  
  if(length(elim)>0){
    ind = ind[-elim]
  }else{
    ind = ind
  }
  #ind = ind[-elim]
  Xr = Xr[,ind]
  rules = rules[ind]
  list(data.matrix(Xr), rules)
}
rule_length = function(rule) {
  splitted = unlist(strsplit(rule, split = "&"))
  length(splitted)
}
calc_prior = function(rules, Xr,alpha, beta) {
  sup = apply(Xr, 2, function(x)mean(x>0))
  prior = c()
  len = unlist(lapply(rules, rule_length))
  sup = sqrt((sup)*(1-sup))
  for(i in 1:length(rules)){
    prior[i] = (len[i]^beta)/((2*sup[i])^(alpha))
    # if(len[i] <= 3){
    #   prior[i] = 1/((2*sup[i])^(alpha))
    # }else{
    #   prior[i] = (len[i]^beta)/((2*sup[i])^(alpha))
    # }
  }
  prior
}
penalty.func = function(rules, Xr,alpha) {
  sup = apply(Xr, 2, function(x)mean(x>0))
  len = unlist(lapply(rules, rule_length))
  sup = sqrt((sup)*(1-sup))
  power.penalty = (len^alpha) - 1
  prior = (1/(2*sup))^(power.penalty)
  prior
}

genrulesGBM = function(X, y, nt, S, L) {
  N = dim(X)[1]
  sf = min(1, (11*sqrt(N)+1)/N)
  mn = 2+floor(rexp(1, 1/(L-2)))
  ns = S
  dist = ifelse(is.numeric(y), "gaussian", "bernoulli")
  if (is.numeric(y)==F){
    y = as.numeric(y)-1
  }
  model1 = gbm.fit(x = X, y=y, bag.fraction = sf,n.trees =1, interaction.depth=(mn/2)
                   ,shrinkage = 0.01,distribution = dist, verbose = F, n.minobsinnode = ns)
  for(i in 2:nt) {
    mn = 2+floor(rexp(1, 1/(L-2)))
    model1$interaction.depth = (mn/2)
    model1 = gbm.more(model1, n.new.trees=1, verbose = F)
  }
  treelist = GBM2List(model1, X)
  rules = extractRules(treelist, X=X, ntree=nt, maxdepth=15)
  rules = c(rules)
  rules = rules[take1(length(rules))]
  rulesmat = matrix(rules)
  colnames(rulesmat) = "condition"
  metric = getRuleMetric(rulesmat,X,y)
  pruned = pruneRule(metric, X, y, 0.025, typeDecay = 1)
  unique(pruned[,4])
}

data.binary.vars = function(n, p, rho, prob){
  # n = sample size
  # p = # of variables
  # rho = correlation between variables
  # prob = probability of each variable 
  
  mu = rep(0, p)
  Sigma = matrix(rho, nrow = p, ncol = p) + diag(p)*(1-rho)
  rawvars = mvrnorm(n=n, mu=mu, Sigma=Sigma)
  pvars = pnorm(rawvars)
  binomvars = qbinom(pvars, 1, 0.5) 
  X = binomvars
  return(X)
}

remove.some.rules = function(X, rules, tau){
  # remove rules with the effects near the population mean
  m.tau = mean(tau)
  N = dim(X)[1]
  num.rules = dim(X)[2]
  
  exc = rep(NA, num.rules)
  for(i in 1:num.rules){
    sub.tau = tau[X[,i]==1]
    m.sub.tau = mean(sub.tau)
    sd.sub.tau = sd(sub.tau)
    lower.tau = m.sub.tau - sd.sub.tau/sqrt(length(sub.tau))
    upper.tau = m.sub.tau + sd.sub.tau/sqrt(length(sub.tau))
    exc[i] = (m.tau >= lower.tau)*(m.tau <= upper.tau)
  }
  
  new.X = X[,exc==0]
  new.rules = rules[exc==0]
  
  return(list(X = new.X, rules = new.rules))
}
glm.lasso <- function (x, y, q, type = c("conservative", "anticonservative"), 
                       ...) 
{
  if (!requireNamespace("glmnet", quietly = TRUE)) 
    stop("Package ", sQuote("glmnet"), " needed but not available")
  # if (is.data.frame(x)) {
  #   message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
  #   x <- model.matrix(~. - 1, x)
  # }
  if ("lambda" %in% names(list(...))) 
    stop("It is not permitted to specify the penalty parameter ", 
         sQuote("lambda"), " for lasso when used with stability selection.")
  type <- match.arg(type)
  if (type == "conservative") 
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q,
                                           intercept = FALSE, 
                                           standardize = FALSE,
                                           ...))
  if (type == "anticonservative") 
    fit <- glmnet::glmnet(x, y, dfmax = q - 1,
                          intercept = FALSE, 
                          standardize = FALSE,
                          ...)
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
}