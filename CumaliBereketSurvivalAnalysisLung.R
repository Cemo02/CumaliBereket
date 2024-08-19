library(survival)
library(survminer)
library(ggplot2)
library(readxl)

lung= read_excel('C:/Users/berek/Desktop/Dataset/lungcancer.xls')
lung
attach(lung)

#survival object Surv(time,event) lub Surv(time1,time2,event,type)
my.survival.object<-Surv(time,status)
my.survival.object
#kaplan-meier estimator
fit<-survfit(my.survival.object ~ 1)
plot(fit, xlab = "days", ylab="Survival")
summary(fit)
summary(fit)$surv     # returns the Kaplan-Meier estimate at each t_i
summary(fit)$time     # {t_i}
summary(fit)$n.risk   # {Y_i}
summary(fit)$n.event  # {d_i}
summary(fit)$std.err  # standard error of the K-M estimate at {t_i}
summary(fit)$lower    # lower pointwise estimates (alternatively, $upper)
str(fit)              # full summary of the fit object
str(summary(fit))     # full summary of the fit object
plot(fit, main="Kaplan-Meier estimate with 95% confidence bounds",
     xlab="time", ylab="survival function")

#KP for sex
fit_sex<-survfit(my.survival.object ~ sex)
ggsurvplot(fit_sex, data = lung)

#cumulative hazard function
my.fit<-summary(fit)
H.hat<--log(my.fit$surv)
H.hat<-c(H.hat, tail(H.hat,1))
plot(H.hat)

#Nelson-Aalen estimator
h.sort.of <- my.fit$n.event / my.fit$n.risk
H.tilde   <- cumsum(h.sort.of)
H.tilde   <- c(H.tilde, tail(H.tilde, 1))
plot(c(my.fit$time, 25), H.hat, xlab="time", ylab="cumulative hazard",
     main="comparing cumulative hazards", type="s")
points(c(my.fit$time, 250), H.tilde, lty=2, type="s")
legend("topleft", legend=c("H.hat","H.tilde"), lty=1:2)

#CoX Hazard
cox=coxph(Surv(time, status)~sex+age+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss
,data=lung)
summary(cox)

#coxph.fit <- coxph(my.survival.object ~sex+ age+thickness+ulcer, method="breslow")
#coxph.fit
#plot the baseline survival function
pl <- survfit(cox)
plot(pl,main="Cox PH",xlab = "Days", ylab="Survival")

#Log-Rank test
survdiff(Surv(time, status) ~ sex, rho=0) #log-rank or Mantel-Haenszel test
survdiff(Surv(time, status) ~ sex, rho=1) #Peto&Peto modification of the Gehan-Wilcoxon test


#AFT models
srFitWeib <- survreg(Surv(time, status) ~ sex+age+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, dist="weibull")
summary(srFitWeib)
srFitWeib$coeff    # covariate coefficients
srFitWeib$var      # variance-covariance matrix
srFitWeib$loglik   # log-likelihood
srFitWeib$scale    

srFitExp <- survreg(Surv(time, status) ~ sex+age+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, dist="exponential")
summary(srFitExp)
srFitExp$coeff    # covariate coefficients
srFitExp$var      # variance-covariance matrix
srFitExp$loglik   # log-likelihood
srFitExp$scale    

srFitLnorm <- survreg(Surv(time, status) ~ sex+age+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, dist="lognormal")
summary(srFitLnorm)
srFitLlog <- survreg(Surv(time, status) ~ sex+age+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, dist="loglogistic")
summary(srFitLlog)




















