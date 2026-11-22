suppressMessages({library(data.table);library(survival);library(splines);library(pROC);library(ggplot2);library(sandwich)})
A<-readRDS("A.rds"); A[,ag:=na-(cl+hco3s)]
M2<-"age+male+dni+paco2+gcs_i+rr_i+map_i+vaso+pna+sepsis+chf+ckd"; M3x<-"+lact_i+creat_i+alb_i+hb_i+wbc_i+hr_i+spo2_i+obes+dm+ster+loop"
fm<-function(db,M) paste(M,if(db=="eICU-CRD") "+sev_i" else if(db=="Pooled") "+dbe" else "")
orr<-function(g,k){s<-summary(g)$coef; c(exp(s[k,1]),exp(s[k,1]-1.96*s[k,2]),exp(s[k,1]+1.96*s[k,2]),s[k,4])}
lg<-function(f,d) suppressWarnings(glm(as.formula(f),binomial,d))
hr<-function(f,d,k,db){if(db=="Pooled") f<-paste(f,"+strata(db)"); c<-coxph(as.formula(f),d); c(exp(coef(c)[k]),exp(confint(c)[k,]),summary(c)$coef[k,5])}
X<-list()
for(db_ in c("Pooled","MIMIC-IV","eICU-CRD")){d<-if(db_=="Pooled") A else A[db==db_]; C2<-fm(db_,M2); C3<-paste(C2,M3x); o<-list()
 o$n<-c(nrow(d),sum(d$fail7),sum(d$hi),sum(d[hi==1]$fail7),sum(d[hi==0]$fail7),sum(d$hi==0))
 o$hiM1<-orr(lg("fail7~hi",d),"hi"); o$hiM2<-orr(lg(paste("fail7~hi+",C2),d),"hi"); o$hiM3<-orr(lg(paste("fail7~hi+",C3),d),"hi"); o$hiHR<-hr(paste("Surv(ft,fail7)~hi+",C2),d,"hi",db_)
 o$zM1<-orr(lg("fail7~z",d),"z"); o$zM2<-orr(lg(paste("fail7~z+",C2),d),"z"); o$zM3<-orr(lg(paste("fail7~z+",C3),d),"z"); o$zHR<-hr(paste("Surv(ft,fail7)~z+",C2),d,"z",db_)
 d$Qn<-as.numeric(d$Q); o$trend<-orr(lg(paste("fail7~Qn+",C2),d),"Qn")
 for(y in c("death7","death28","hdeath")){o[[paste0(y,"_hi")]]<-c(orr(lg(paste(y,"~hi+",C2),d),"hi"),100*mean(d[hi==0][[y]],na.rm=T),100*mean(d[hi==1][[y]],na.rm=T))}
 o$sd<-sd(d$dh); X[[db_]]<-o}
d<-A; C2<-fm("Pooled",M2)
X$s048<-{dd<-copy(d)[,dh2:=hco3-(24+0.48*(paco2-40))][,z:=-dh2/sd(dh2),by=db]; orr(lg(paste("fail7~z+",C2),dd),"z")}
X$copd<-{dd<-d[copd==1]; c(orr(lg(paste("fail7~z+",sub("chf","chf",C2)),dd),"z"),nrow(dd),sum(dd$fail7))}
X$noncopd<-{dd<-d[copd==0]; c(orr(lg(paste("fail7~z+",C2),dd),"z"),nrow(dd),sum(dd$fail7))}
X$excl_early<-{dd<-d[!(fail7==1 & ft<=6)]; c(orr(lg(paste("fail7~z+",C2),dd),"z"),nrow(dd))}
X$cox_z_strat<-X$Pooled$zHR
# RCS on raw dh pooled
b0<-glm(as.formula(paste("fail7~",C2)),binomial,d); bn<-glm(as.formula(paste("fail7~ns(dh,3)+",C2)),binomial,d); bl<-glm(as.formula(paste("fail7~dh+",C2)),binomial,d)
X$rcs_raw<-c(anova(b0,bn,test="LRT")[2,5],anova(bl,bn,test="LRT")[2,5])
ref<-median(d$dh); g<-data.table(dh=seq(quantile(d$dh,.02),quantile(d$dh,.98),length.out=200)); nd<-d[rep(1,200)]; nd[,dh:=g$dh]
for(v in c("age","paco2","gcs_i","rr_i","map_i")) set(nd,j=v,value=median(d[[v]])); for(v in c("male","dni","vaso","pna","sepsis","chf","ckd","dbe")) set(nd,j=v,value=0)
pr<-predict(bn,nd,se.fit=TRUE); nr<-copy(nd[1])[,dh:=ref]; p0<-predict(bn,nr)
X$rcsdat<-data.table(dh=g$dh,or=exp(pr$fit-p0),lo=exp(pr$fit-p0-1.96*pr$se.fit),hi=exp(pr$fit-p0+1.96*pr$se.fit)); X$ref<-ref
r1<-roc(d$fail7,d$z,quiet=TRUE); r2<-roc(d$fail7,d$phz,quiet=TRUE); X$roc<-list(a1=ci.auc(r1),a2=ci.auc(r2),p=roc.test(r1,r2)$p.value,r1=r1,r2=r2)
saveRDS(X,"X.rds")
s<-function(x) signif(x,3); for(k in c("Pooled","MIMIC-IV","eICU-CRD")){o<-X[[k]]; cat("\n",k,"n/ev/hi/evhi/evlo/lo",o$n,"\n"); for(n in c("hiM1","hiM2","hiM3","hiHR","zM1","zM2","zM3","zHR","trend","death7_hi","death28_hi","hdeath_hi")) cat(n,s(o[[n]]),"\n"); cat("sd",s(o$sd),"\n")}
cat("s048",s(X$s048),"\ncopd",s(X$copd),"\nnoncopd",s(X$noncopd),"\nexcl_early",s(X$excl_early),"\nrcs_raw",s(X$rcs_raw),"ref",X$ref,"\nroc",s(X$roc$a1),s(X$roc$a2),s(X$roc$p),"\n")
