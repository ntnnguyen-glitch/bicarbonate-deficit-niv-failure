suppressMessages({library(data.table);library(survival);library(mice)}); set.seed(2026)
source(textConnection(paste(readLines("s2.R")[1:10],collapse="\n")))
D<-readRDS("D.rds")$d; A0<-readRDS("A.rds")
A<-copy(A0)
M2<-"age+male+dni+paco2+gcs_i+rr_i+map_i+vaso+pna+sepsis+chf+ckd+dbe"; M3<-paste(M2,"+lact_i+creat_i+alb_i+hb_i+wbc_i+hr_i+spo2_i+obes+dm+ster+loop")
orr<-function(g,k){s<-summary(g)$coef; round(c(exp(s[k,1]),exp(s[k,1]-1.96*s[k,2]),exp(s[k,1]+1.96*s[k,2]),s[k,4]),4)}
lg<-function(f,d) suppressWarnings(glm(as.formula(f),binomial,d))
O<-list()
O$src<-A[,.N,by=.(db,src,hi)][order(db,src,hi)]
O$disch<-A[,.(n=.N,alive_disch_lt7=sum(hdeath==0 & disch_h<168,na.rm=T),alive_disch_lt28=sum(hdeath==0 & disch_h<672,na.rm=T),death_after_disch_le7=sum(!is.na(death_h)&death_h<=168&death_h>disch_h)),by=db]
A[,te:=pmin(intub_h,fifelse(death7==1,death_h,NA_real_),na.rm=TRUE)]
O$evt_after_disch<-A[fail7==1 & te>disch_h,.N,by=db]
O$evt_before_abg<-A[fail7==1 & te<abg_h,.N,by=db]
O$abg_timing<-A[,.(med=median(abg_h),q1=quantile(abg_h,.25),q3=quantile(abg_h,.75),post=sum(abg_h>0)),by=db]
# corrected follow-up: events at true time; non-events followed to 168 h (MIMIC, out-of-hospital death known) or censored at discharge (eICU)
A[,ft2:=fifelse(fail7==1,te,fifelse(db=="MIMIC-IV",168,pmin(168,pmax(disch_h,0.25))))]; A[ft2<=0,ft2:=0.25]
c1<-coxph(as.formula(paste("Surv(ft2,fail7)~hi+",sub("\\+dbe","",M2),"+strata(db)")),A); c2<-coxph(as.formula(paste("Surv(ft2,fail7)~z+",sub("\\+dbe","",M2),"+strata(db)")),A)
O$cox2<-rbind(hi=c(exp(coef(c1)["hi"]),exp(confint(c1)["hi",]),summary(c1)$coef["hi",5]),z=c(exp(coef(c2)["z"]),exp(confint(c2)["z",]),summary(c2)$coef["z",5]))
O$lr2<-c(pooled=survdiff(Surv(ft2,fail7)~Q,A)$pvalue,m=survdiff(Surv(ft2,fail7)~Q,A[db=="MIMIC-IV"])$pvalue,e=survdiff(Surv(ft2,fail7)~Q,A[db=="eICU-CRD"])$pvalue)
for(db_ in c("MIMIC-IV","eICU-CRD")){d<-A[db==db_]; f<-sub("\\+dbe",if(db_=="eICU-CRD") "+sev_i" else "",M2); c3<-coxph(as.formula(paste("Surv(ft2,fail7)~hi+",f)),d); O[[paste0("cox2_",db_)]]<-c(exp(coef(c3)["hi"]),exp(confint(c3)["hi",]),summary(c3)$coef["hi",5])}
# pre-NIV ABG only
P<-A[abg_h<=0]; O$pre<-list(n=nrow(P),ev=sum(P$fail7),hi=orr(lg(paste("fail7~hi+",M2),P),"hi"),z=orr(lg(paste("fail7~z+",M2),P),"z"))
# DNI documented before NIV
f2<-sub("dni","dni_pre",M2); O$dnipre<-list(n=sum(A$dni_pre),adj_hi=orr(lg(paste("fail7~hi+",f2),A),"hi"),excl_hi=orr(lg(paste("fail7~hi+",f2),A[dni_pre==0]),"hi"),excl_z=orr(lg(paste("fail7~z+",f2),A[dni_pre==0]),"z"),nexcl=nrow(A[dni_pre==0]))
# no post-baseline covariates (drop sedv) in Model 3
O$m3_nosed<-orr(lg(paste("fail7~hi+",sub("\\+sedv","",M3)),A),"hi")
# MICE binary
raw<-c("lact","creat","alb","hb","wbc","hr","spo2","gcs","rr","map")
dm<-A[,c("fail7","hi","age","male","dni","paco2","vaso","pna","sepsis","chf","ckd","obes","dm","ster","loop",raw,"dbe"),with=FALSE]
im<-mice(dm,m=10,maxit=5,printFlag=FALSE); fit<-with(im,glm(as.formula(paste("fail7~hi+age+male+dni+paco2+vaso+pna+sepsis+chf+ckd+obes+dm+ster+loop+",paste(raw,collapse="+"),"+dbe")),family=binomial)); s<-summary(pool(fit)); s<-s[s$term=="hi",]
O$mice_hi<-round(c(exp(s$estimate),exp(s$estimate-1.96*s$std.error),exp(s$estimate+1.96*s$std.error),s$p.value),4)
O$miss<-A[,.(lact=mean(is.na(lact)),alb=mean(is.na(alb))),by=hi]
O$lact_measured_vs_fail<-A[,.(fail=mean(fail7)),by=.(lact_measured=!is.na(lact))]
# DeltaHCO3 and pH in the same model
O$withph<-orr(lg(paste("fail7~z+phz+",M2),A),"z")
O$Rver<-R.version.string
O$inhosp_fail<-A[,.(n=sum((!is.na(intub_h)&intub_h<=168)|(hdeath==1&!is.na(death_h)&death_h<=168&death_h<=disch_h+24))),by=db]
saveRDS(O,"O.rds"); saveRDS(A,"A2.rds"); print(O)
