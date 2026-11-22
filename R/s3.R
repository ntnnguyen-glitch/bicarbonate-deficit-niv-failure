suppressMessages({library(data.table);library(survival);library(cmprsk);library(splines);library(WeightIt);library(MatchIt);library(sandwich);library(mice);library(pROC)})
set.seed(2026); source(textConnection(paste(readLines("s2.R")[1:10],collapse="\n")))
D<-readRDS("D.rds")$d; D[,dni_any:=dni][,dni:=dni_pre]; A<-prep(D,FALSE,TRUE); A[,ag:=na-(cl+hco3s)]
A[,z:=-dh/sd(dh),by=db]; A[,Q:=factor(cut(-dh,quantile(-dh,0:4/4),include.lowest=TRUE,labels=FALSE)),by=db]; A[,hi:=as.integer(Q=="4")]
A[,phz:=-(ph-mean(ph))/sd(ph),by=db]; A[,pcz:=(paco2-mean(paco2))/sd(paco2),by=db]
A[,te:=pmin(intub_h,fifelse(death7==1,death_h,NA_real_),na.rm=TRUE)]; A[,ft:=fifelse(fail7==1,te,fifelse(db=="MIMIC-IV",168,pmin(168,pmax(disch_h,0.25))))]; A[ft<=0,ft:=0.25]
A[,ev:=fifelse(intub7==1,1L,fifelse(death7==1,2L,0L))]; A[,intub48:=as.integer(!is.na(intub_h)&intub_h<=48)]
A[,icufd:=fifelse(death28==1,0,pmax(0,pmin(28,28-icuout_h/24)))]; A[,age65:=as.integer(age>=65)]
M2<-"age+male+dni+paco2+gcs_i+rr_i+map_i+vaso+pna+sepsis+chf+ckd"; M3x<-"+lact_i+creat_i+alb_i+hb_i+wbc_i+hr_i+spo2_i+obes+dm+ster+loop"
fm<-function(db,M) paste(M,if(db=="eICU-CRD") "+sev_i" else if(db=="Pooled") "+dbe" else "")
orr<-function(g,k){s<-summary(g)$coef; c(exp(s[k,1]),exp(s[k,1]-1.96*s[k,2]),exp(s[k,1]+1.96*s[k,2]),s[k,4])}
fmt<-function(x) sprintf("%.2f (%.2f–%.2f)",x[1],x[2],x[3]); lg<-function(f,d) suppressWarnings(glm(as.formula(f),binomial,d))
hrf<-function(f,d,k,db){if(db=="Pooled") f<-paste(f,"+strata(db)"); c<-coxph(as.formula(f),d); c(exp(coef(c)[k]),exp(confint(c)[k,]),summary(c)$coef[k,5])}
run<-function(d,db){C2<-fm(db,M2); C3<-paste(C2,M3x); o<-list(n=nrow(d)); o$rate<-c(fail=mean(d$fail7),intub=mean(d$intub7),death=mean(d$death7)); o$dh<-quantile(d$dh,c(.25,.5,.75))
 t1<-d[,.(n=.N,ev=sum(fail7),pct=100*mean(fail7),rng=sprintf("%.1f to %.1f",min(dh),max(dh))),by=Q][order(Q)]
 g1<-lg("fail7~Q",d); g2<-lg(paste("fail7~Q+",C2),d); g3<-lg(paste("fail7~Q+",C3),d)
 for(k in 2:4){nm<-paste0("Q",k); h<-hrf(paste("Surv(ft,fail7)~Q+",C2),d,nm,db); t1[k,`:=`(c1=fmt(orr(g1,nm)),c2=fmt(orr(g2,nm)),c3=fmt(orr(g3,nm)),p2=orr(g2,nm)[4],hr=fmt(h),ph=h[4])]}
 o$t1<-t1; d$Qn<-as.numeric(d$Q); o$trend<-orr(lg(paste("fail7~Qn+",C2),d),"Qn"); o$km<-survdiff(Surv(ft,fail7)~Q,d)$pvalue
 o$t1b<-rbindlist(lapply(list(c("z","−ΔHCO3 (mỗi 1-SD, bù trừ kém hơn)"),c("phz","pH thấp hơn (mỗi 1-SD)"),c("pcz","PaCO2 cao hơn (mỗi 1-SD)")),function(a){
  f2<-if(a[1]=="pcz") sub("paco2","1",C2) else C2; f3<-if(a[1]=="pcz") sub("paco2","1",C3) else C3
  gc<-lg(paste("fail7~",a[1]),d); g2<-lg(paste("fail7~",a[1],"+",f2),d); g3<-lg(paste("fail7~",a[1],"+",f3),d); h<-hrf(paste("Surv(ft,fail7)~",a[1],"+",f2),d,a[1],db)
  data.table(x=a[2],orc=fmt(orr(gc,a[1])),or2=fmt(orr(g2,a[1])),or3=fmt(orr(g3,a[1])),hr=fmt(h),auc=as.numeric(auc(d$fail7,d[[a[1]]],quiet=TRUE)),p=orr(g2,a[1])[4],p3=orr(g3,a[1])[4])}))
 o$delong<-roc.test(roc(d$fail7,d$z,quiet=TRUE),roc(d$fail7,d$phz,quiet=TRUE))$p.value
 o$beyond_ph<-orr(lg(paste("fail7~z+phz+",C2),d),"z")
 b0<-glm(as.formula(paste("fail7~",C2)),binomial,d); bl<-glm(as.formula(paste("fail7~z+",C2)),binomial,d); bn<-glm(as.formula(paste("fail7~ns(z,3)+",C2)),binomial,d)
 o$rcs<-c(anova(b0,bn,test="LRT")[2,5],anova(bl,bn,test="LRT")[2,5]); o$inc<-c(anova(b0,bl,test="LRT")[2,5],auc(d$fail7,fitted(b0),quiet=TRUE),auc(d$fail7,fitted(bl),quiet=TRUE))
 X<-model.matrix(as.formula(paste("~z+",C2)),d)[,-1]; X<-X[,apply(X,2,var)>0,drop=FALSE]; s<-summary(crr(d$ft,d$ev,X,failcode=1,cencode=0))$coef
 o$fg<-c(s["z",2],exp(s["z",1]-1.96*s["z",3]),exp(s["z",1]+1.96*s["z",3]),s["z",5]); o$gray<-cuminc(d$ft,d$ev,d$Q)$Tests[1,2]
 so<-c(intub48="Đặt NKQ 48 giờ",intub7="Đặt NKQ 7 ngày",death7="Tử vong 7 ngày",death28="Tử vong 28 ngày",hdeath="Tử vong nội viện")
 o$t2<-rbindlist(lapply(names(so),function(y){g<-lg(paste(y,"~Q+",C2),d); gz<-lg(paste(y,"~z+",C2),d)
  data.table(y=so[y],r=sprintf("%.1f%% vs %.1f%%",100*mean(d[Q=="1"][[y]],na.rm=TRUE),100*mean(d[Q=="4"][[y]],na.rm=TRUE)),q4=fmt(orr(g,"Q4")),pq=orr(g,"Q4")[4],sd=fmt(orr(gz,"z")),psd=orr(gz,"z")[4])}))
 gl<-lm(as.formula(paste("icufd~z+",C2)),d); s<-summary(gl)$coef["z",]; o$icufd<-c(s[1],s[1]-1.96*s[2],s[1]+1.96*s[2],s[4])
 ps<-as.formula(paste("hi~",C3)); for(est in c("ATE","ATO")){w<-weightit(ps,data=d,method="glm",estimand=est); ww<-pmin(w$weights,quantile(w$weights,.99))
  g<-suppressWarnings(glm(fail7~hi,quasibinomial,d,weights=ww)); v<-sqrt(diag(vcovHC(g,"HC0")))["hi"]; b<-coef(g)["hi"]; o[[est]]<-c(exp(b),exp(b-1.96*v),exp(b+1.96*v),2*pnorm(-abs(b/v)))
  o[[paste0("smd_",est)]]<-max(abs(cobalt::bal.tab(ps,data=d,weights=ww,s.d.denom="pooled",binary="std")$Balance$Diff.Adj[-1]))}
 mt<-matchit(ps,data=d,method="nearest",caliper=0.2); md<-match.data(mt); g<-suppressWarnings(glm(fail7~hi,quasibinomial,md,weights=weights)); v<-sqrt(diag(vcovCL(g,cluster=~subclass)))["hi"]; b<-coef(g)["hi"]
 o$psm<-c(exp(b),exp(b-1.96*v),exp(b+1.96*v),2*pnorm(-abs(b/v)),nrow(md)); o$hi<-orr(lg(paste("fail7~hi+",C3),d),"hi")
 raw<-c("lact","creat","alb","hb","wbc","hr","spo2","gcs","rr","map",if(db!="MIMIC-IV") "sev"); dm<-d[,c("fail7","z","age","male","dni","paco2","vaso","pna","sepsis","chf","ckd","obes","dm","ster","loop",raw,if(db=="Pooled") "dbe"),with=FALSE]
 im<-mice(dm,m=10,maxit=5,printFlag=FALSE); fit<-with(im,glm(as.formula(paste("fail7~z+age+male+dni+paco2+vaso+pna+sepsis+chf+ckd+obes+dm+ster+loop+",paste(raw,collapse="+"),if(db=="Pooled") "+dbe")),family=binomial)); s<-summary(pool(fit)); s<-s[s$term=="z",]
 o$mice<-c(exp(s$estimate),exp(s$estimate-1.96*s$std.error),exp(s$estimate+1.96*s$std.error),s$p.value)
 ev_<-function(r) if(r<=1) 1 else r+sqrt(r*(r-1)); o$evalue<-c(ev_(sqrt(o$hi[1])),ev_(sqrt(o$hi[2])))
 sens<-list("Loại bệnh nhân DNI"=list(d[dni_any==0],C2),"Chỉ khí máu trước NIV"=list(d[abg_h<=0],C2),"Hiệu chỉnh thêm anion gap"=list(d[!is.na(ag)],paste(C2,"+ag")),"Hiệu chỉnh thêm lactate"=list(d,paste(C2,"+lact_i")),
  "Loại bệnh nhân dùng lợi tiểu quai 24 h trước"=list(d[loop==0],C2),"Chỉ HCO3 khí máu động mạch"=list(d[src=="ABG"],C2),"Chỉ HCO3 huyết thanh"=list(d[!is.na(hco3s)][,z:=-(hco3s-(24+0.35*(paco2-40)))/sd(hco3s-(24+0.35*(paco2-40)))],C2),
  "Công thức kỳ vọng 0,4×ΔPaCO2"=list(copy(d)[,z:=-(hco3-(24+0.4*(paco2-40)))/sd(hco3-(24+0.4*(paco2-40)))],C2),"Model 3 đầy đủ"=list(d,C3))
 o$sens<-rbindlist(lapply(names(sens),function(k){a<-sens[[k]]; if(nrow(a[[1]])<100) return(NULL); g<-lg(paste("fail7~z+",a[[2]]),a[[1]]); r<-orr(g,"z"); data.table(a=k,n=nrow(a[[1]]),or=fmt(r),p=r[4])}))
 if(db!="MIMIC-IV"){g<-lg(paste("fail7~z+",C2),d); v<-sqrt(diag(vcovCL(g,cluster=~hosp)))["z"]; b<-coef(g)["z"]; o$sens<-rbind(o$sens,data.table(a="Sai số chuẩn robust cụm theo bệnh viện",n=nrow(d),or=fmt(exp(c(b,b-1.96*v,b+1.96*v))),p=2*pnorm(-abs(b/v))))}
 sg<-list("Giới"=c("male","Nữ","Nam"),"Tuổi"=c("age65","<65","≥65"),"COPD"=c("copd","Không","Có"),"Viêm phổi"=c("pna","Không","Có"),"Nhiễm khuẩn huyết"=c("sepsis","Không","Có"),"Suy tim"=c("chf","Không","Có"))
 o$sub<-rbindlist(lapply(names(sg),function(k){a<-sg[[k]]; C<-gsub(paste0("\\b",a[1],"\\+|\\+",a[1],"\\b"),"",C2); if(a[1]=="age65") C<-sub("age\\+","",C)
  pint<-anova(glm(as.formula(paste("fail7~z+",a[1],"+",C)),binomial,d),glm(as.formula(paste("fail7~z*",a[1],"+",C)),binomial,d),test="LRT")[2,5]
  rbindlist(lapply(0:1,function(l){dd<-d[get(a[1])==l]; r<-orr(lg(paste("fail7~z+",C),dd),"z"); data.table(g=k,l=a[l+2],n=nrow(dd),ev=sprintf("%d (%.1f%%)",sum(dd$fail7),100*mean(dd$fail7)),or=fmt(r),p=r[4],pint=pint)}))}))
 o}
R<-list(); for(db_ in c("MIMIC-IV","eICU-CRD","Pooled")){d<-if(db_=="Pooled") copy(A) else A[db==db_]; R[[db_]]<-run(d,db_); cat(db_,"done\n")}
R$flow<-readRDS("D.rds")$flow; R$nflow<-list(m=c(abg=nrow(D[db=="MIMIC-IV"]),hyp=nrow(D[db=="MIMIC-IV"&paco2>=45&!is.na(hco3)&!is.na(ph)]),acid=nrow(A[db=="MIMIC-IV"])),e=c(abg=nrow(D[db=="eICU-CRD"]),hyp=nrow(D[db=="eICU-CRD"&paco2>=45&!is.na(hco3)&!is.na(ph)]),acid=nrow(A[db=="eICU-CRD"])))
saveRDS(R,"R.rds"); saveRDS(A,"A.rds")
s<-function(x) signif(x,3); for(db in c("MIMIC-IV","eICU-CRD","Pooled")){o<-R[[db]]; cat("\n==",db,o$n,"\n"); print(o$t1[,.(Q,n,pct=round(pct,1),c2,c3,hr,p2=s(p2),ph=s(ph))]); print(o$t1b[,.(x,or2,or3,hr,auc=round(auc,3),p=s(p),p3=s(p3))])
 cat("trend",s(o$trend),"KM",s(o$km),"delong",s(o$delong),"beyondpH",s(o$beyond_ph),"rcs",s(o$rcs),"inc",s(o$inc),"\nFG",s(o$fg),"gray",s(o$gray),"icufd",s(o$icufd),"\n"); print(o$t2[,.(y,r,q4,pq=s(pq),sd,psd=s(psd))])
 cat("ATE",s(o$ATE),s(o$smd_ATE),"ATO",s(o$ATO),"PSM",s(o$psm),"MICE",s(o$mice),"hi",s(o$hi),"E",s(o$evalue),"\n"); print(o$sens[,.(a,n,or,p=s(p))]); print(o$sub[,.(g,l,n,or,p=s(p),pint=s(pint))])}
print(R$nflow)
