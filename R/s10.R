suppressMessages({library(data.table);library(mice)})
A<-as.data.table(readRDS("A.rds"))
raw<-c("lact","creat","alb","hb","wbc","hr","spo2","gcs","rr","map")
m2<-c("age","male","dni","paco2","gcs","rr","map","vaso","pna","sepsis","chf","ckd","dbe")
m3<-c(m2,"lact","creat","alb","hb","wbc","hr","spo2","obes","dm","ster","loop")
dm<-A[,unique(c("fail7","hi","z",m3)),with=FALSE]
im<-mice(dm,m=50,maxit=10,printFlag=FALSE,seed=2026)
pl<-function(x,cv){fit<-with(im,glm(as.formula(paste("fail7~",x,"+",paste(cv,collapse="+"))),family=binomial)); s<-summary(pool(fit)); s<-s[s$term==x,]
 c(exp(s$estimate),exp(s$estimate-1.96*s$std.error),exp(s$estimate+1.96*s$std.error),s$p.value)}
M<-list(m3_hi=pl("hi",m3),m3_z=pl("z",m3),m2_hi=pl("hi",m2),m2_z=pl("z",m2))
cc<-na.omit(A[,c("fail7","hi","z",m3),with=FALSE]); M$cc_n<-nrow(cc); M$cc_ev<-sum(cc$fail7)
cf<-function(x){g<-glm(as.formula(paste("fail7~",x,"+",paste(m3,collapse="+"))),binomial,cc); s<-summary(g)$coef[x,]; c(exp(s[1]),exp(s[1]-1.96*s[2]),exp(s[1]+1.96*s[2]),s[4])}
M$cc_hi<-cf("hi"); M$cc_z<-cf("z")
ev_<-function(r) if(r<=1) 1 else r+sqrt(r*(r-1)); M$evalue<-c(ev_(sqrt(M$m3_hi[1])),ev_(sqrt(M$m3_hi[2])))
M$miss_m2<-sapply(c("gcs","rr","map"),function(v) 100*mean(is.na(A[[v]])))
saveRDS(M,"MI.rds"); print(lapply(M,function(x) signif(x,4)))
