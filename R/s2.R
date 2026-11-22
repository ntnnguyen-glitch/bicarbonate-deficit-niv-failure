suppressMessages(library(data.table))
prep<-function(D,copd_only=TRUE,acid=TRUE){d<-copy(D)[paco2>=45 & !is.na(hco3) & !is.na(ph)]; if(acid) d<-d[ph<7.35]; if(copd_only) d<-d[copd==1]
 d[,exp_hco3:=24+0.35*(paco2-40)]; d[,dh:=hco3-exp_hco3]
 d[,cat3:=factor(fifelse(dh< -4,"Inadequate",fifelse(dh>4,"Excess","Appropriate")),levels=c("Appropriate","Inadequate","Excess"))]
 d[,intub7:=as.integer(!is.na(intub_h)&intub_h<=168)]; d[,death7:=as.integer(!is.na(death_h)&death_h<=168)]; d[,fail7:=as.integer(intub7==1|death7==1)]
 d[,death28:=as.integer(!is.na(death_h)&death_h<=672)]
 for(v in c("lact","creat","alb","hb","wbc","rr","hr","map","gcs","spo2","sev","wt","glu","k")){x<-d[[v]];x[is.na(x)]<-median(x,na.rm=TRUE);set(d,j=paste0(v,"_i"),value=x)}
 d[,dbe:=as.integer(db=="eICU-CRD")]; d}
D<-readRDS("D.rds")$d; d<-prep(D)
print(d[,.(n=.N,fail=mean(fail7),intub=mean(intub7),death=mean(death7),dh=median(dh),inad=mean(cat3=="Inadequate"),exc=mean(cat3=="Excess")),by=db])
C3<-"age+male+dni+obes+chf+ckd+pna+sepsis+dm+paco2+lact_i+creat_i+alb_i+hb_i+wbc_i+rr_i+hr_i+map_i+gcs_i+spo2_i+vaso+ster+loop+sedv"
r<-function(y,x,dd,f){g<-glm(as.formula(paste(y,"~",x,"+",f)),binomial,dd); s<-summary(g)$coef; k<-grep(paste0("^",x),rownames(s)); for(j in k) cat(sprintf("%-10s %-8s %-18s %5d OR %.2f (%.2f-%.2f) p=%.3g\n",unique(dd$db)[1],y,rownames(s)[j],nrow(dd),exp(s[j,1]),exp(s[j,1]-1.96*s[j,2]),exp(s[j,1]+1.96*s[j,2]),s[j,4]))}
d[,dhz:=-dh/sd(dh)]  # per SD lower compensation
for(db_ in c("MIMIC-IV","eICU-CRD")){dd<-d[db==db_]; f<-if(db_=="eICU-CRD") paste(C3,"+sev_i") else C3
 for(y in c("fail7","intub7","death7","death28")){r(y,"dhz",dd,"1"); r(y,"dhz",dd,f); r(y,"cat3",dd,f)}}
dd<-d; r("fail7","dhz",dd,paste(C3,"+dbe")); r("fail7","cat3",dd,paste(C3,"+dbe"))
# vs pH alone (head-to-head)
for(db_ in c("MIMIC-IV","eICU-CRD")){dd<-d[db==db_]; dd[,phz:=-(ph-mean(ph))/sd(ph)]; r("fail7","phz",dd,C3)}
