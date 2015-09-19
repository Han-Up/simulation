## this program compares OLS, fixed effects (one-way, two-way),
## and First-Difference estimators w/ or w/o one-way fixed effects with simulated data.
## the progrma depends on data.table, lfe, doParallel, and foreach
## Han-Up Park, Sep 2015
## if packages missing in your workstation, run the following:
## install.package("data.table")
## install.package("lfe")
## install.package("foreach")
library(data.table)
library(lfe)
library(foreach)
library(doParallel)


## executing 1,000 simulations
#setup parallel backend to use 8 processors##!!set to the number of processing cores of your workstation
cl<-makeCluster(8)
registerDoParallel(cl)


case1<-foreach(k=1:1000 , .combine=rbind, .packages=c("lfe","data.table")) %dopar% {
	##T: # of years, N: # of firms
	##T=20, N=5000, Balanced panels
	N=5000
	T=20

	## year and firm indicators
	year <- factor(rep(seq(1,T,1),N))
	firm <- factor(rep(seq(1,N,1),T)); firm <- firm[order(firm)]



	## create covariates correlated with fixed effects
	## x is correlated with firm FE, x2 is correlated with year FE
	data = data.table(y_fe=rep(0,N*T),x=0,x2=0,fe_firm = 0,fe_year = 0,u=0,lag.y=0,lag.x=0,lag.x2=0,d.y_fe=0,d.x=0,d.x2=0)
     
	   ## assume effects for them, Fixed Effects for each firm and each year
	data[,firm.eff:=rnorm(1),by=firm]
	data[,year.eff:=rnorm(1),by=year]
	## normal x with mean of x's fixed effects
	data[, x:=rnorm(T,mean = -firm.eff, sd=1), by=firm]
	data[, x2:=rnorm(N,mean = -year.eff, sd=1), by=year]
	
	## the strictly exogenous error
	data[,u:= rnorm(length(x)),]
	
	## Case 1: only firm fixed effects exists
	data[,y_fe:= 0.5*x + firm.eff + u,]

	## first-differences
		## lagged observables
	data[, lag.y_fe:=c(NA, y_fe[-.N]), by=firm]
	data[, lag.x:=c(NA, x[-.N]), by=firm]
	data[, lag.x2:=c(NA, x2[-.N]), by=firm]
		## first differences i.e. change specification
	data[, d.y_fe:=y_fe - lag.y_fe, by=firm]
	data[, d.x:=x - lag.x, by=firm]
	data[, d.x2:=x - lag.x2, by=firm]
	## OLS, FE, FD
	ols= lm(y_fe ~ x + 1,data)
	fe = felm(y_fe ~ x  | firm | 0 | 0 ,data)
	fd = lm(d.y_fe ~ d.x , data)
	## although numerically FE equivalent to FD only in T=2, both are unbiased and consistent for largge N
	
	coef_x_ols <- coef(ols)[2]
	coef_x_fe <- coef(fe)[1]
	coef_x_fd <- coef(fd)[2]
	coef_adjr2_ols <- summary(ols)$adj.r.squared
	coef_adjr2_fe <- summary(fe)$P.adj.r.squared
	coef_adjr2_fd <- summary(fd)$adj.r.squared
    return(c(coef_x_ols,coef_x_fe,coef_x_fd,coef_adjr2_ols,coef_adjr2_fe,coef_adjr2_fd))
    }
	##results for case 1
##OLS coef, FE coef, FD coef, OLS R2, FE R2, FD R2
rbind(colMeans(case1),sqrt(diag(var(case1))))


## case2: both firm-fixed effects and year-fixed effects present

case2<-foreach(k=1:1000 , .combine=rbind, .packages=c("lfe","data.table")) %dopar% {
	##T: # of years, N: # of firms
	##T=20, N=5000, Balanced panels
	N=5000
	T=20

	## year and firm indicators
	year <- factor(rep(seq(1,T,1),N))
	firm <- factor(rep(seq(1,N,1),T)); firm <- firm[order(firm)]

	## create covariates correlated with fixed effects
	## x is correlated with firm FE, x2 is correlated with year FE
	data = data.table(y_fe=rep(0,N*T),x=0,x2=0,fe_firm = 0,fe_year = 0,u=0,lag.y=0,lag.x=0,lag.x2=0,d.y_fe=0,d.x=0,d.x2=0)
     
	   ## assume effects for them, Fixed Effects for each firm and each year
	data[,firm.eff:=rnorm(1),by=firm]
	data[,year.eff:=rnorm(1),by=year]
	## normal x with mean of x's fixed effects
	data[, x:=rnorm(T,mean = -firm.eff, sd=1), by=firm]
	data[, x2:=rnorm(N,mean = -year.eff, sd=1), by=year]
	
	## the strictly exogenous error
	data[,u:= rnorm(length(x)),]
	
	## Case 1: only firm fixed effects exists
	data[,y_fe:= 0.5*x + firm.eff + 0.5*x2 + year.eff + u,]

	## first-differences
		## lagged observables
	data[, lag.y_fe:=c(NA, y_fe[-.N]), by=firm]
	data[, lag.x:=c(NA, x[-.N]), by=firm]
	data[, lag.x2:=c(NA, x2[-.N]), by=firm]
		## first differences i.e. change specification
	data[, d.y_fe:=y_fe - lag.y_fe, by=firm]
	data[, d.x:=x - lag.x, by=firm]
	data[, d.x2:=x2 - lag.x2, by=firm]
	## OLS, FE, FD
	ols= lm(y_fe ~ x + x2 + 1,data)
	fe = felm(y_fe ~ x + x2 | firm + year| 0 | 0 ,data)
	fd = felm(d.y_fe ~ d.x +d.x2 | year|0|0, data)

	coef_x_ols1 <- coef(ols)[2]
	coef_x_ols2 <- coef(ols)[3]
	coef_x_fe1 <- coef(fe)[1]
	coef_x_fe2 <- coef(fe)[2]
	coef_x_fd1 <- coef(fd)[1]
	coef_x_fd2 <- coef(fd)[2]
	coef_adjr2_ols <- summary(ols)$adj.r.squared
	coef_adjr2_fe <- summary(fe)$P.adj.r.squared
	coef_adjr2_fd <- summary(fd)$P.adj.r.squared
    return(c(coef_x_ols1,coef_x_ols2,coef_x_fe1,coef_x_fe2,coef_x_fd1,coef_x_fd2,coef_adjr2_ols,coef_adjr2_fe,coef_adjr2_fd))
    }

	##results for case 2
##OLS X coef, OLS Z coef, FE X coef,FE Z coef, FD X coef, FD Z coef, OLS R2, FE R2, FD R2
rbind(colMeans(case2),sqrt(diag(var(case2))))
 
	## stop parallel processing backend
stopCluster(cl)

