suppressMessages({library(data.table);library(mice)}); set.seed(2026)
A<-readRDS("A.rds"); X<-readRDS("X.rds")
raw<-c("lact","creat","alb","hb","wbc","hr","spo2","gcs","rr","map")
dm<-A[,c("fail7","z","age","male","dni","paco2","vaso","pna","sepsis","chf","ckd","obes","dm","ster","loop",raw,"dbe"),with=FALSE]
im<-mice(dm,m=10,maxit=5,printFlag=FALSE); fit<-with(im,glm(as.formula(paste("fail7~z+age+male+dni+paco2+vaso+pna+sepsis+chf+ckd+obes+dm+ster+loop+",paste(raw,collapse="+"),"+dbe")),family=binomial)); s<-summary(pool(fit)); s<-s[s$term=="z",]
X$mice2<-c(exp(s$estimate),exp(s$estimate-1.96*s$std.error),exp(s$estimate+1.96*s$std.error),s$p.value)
X$nhosp<-uniqueN(A[db=="eICU-CRD"]$hosp)
X$comp<-A[,.(n=.N,fail=sum(fail7),intub=sum(intub7),death=sum(death7),dni=sum(dni),copd=sum(copd),abg=sum(src=="ABG")),by=db]
X$compQ<-A[,.(n=.N,intub=sum(intub7),death=sum(death7)),by=hi]
saveRDS(X,"X.rds"); print(signif(X$mice2,3)); print(X$nhosp); print(X$comp); print(X$compQ)
