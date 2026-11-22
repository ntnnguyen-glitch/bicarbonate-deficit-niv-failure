suppressMessages({library(data.table);library(survival);library(ggplot2);library(pROC)})
A<-readRDS("A.rds"); R<-readRDS("R.rds"); X<-readRDS("X.rds"); MI<-readRDS("MI.rds")
# ---------------- Table 1
f1<-function(x,d=1) sprintf(paste0("%.",d,"f (%.",d,"f–%.",d,"f)"),median(x,na.rm=T),quantile(x,.25,na.rm=T),quantile(x,.75,na.rm=T))
fb<-function(x) sprintf("%d (%.1f)",sum(x,na.rm=T),100*mean(x,na.rm=T))
smd<-function(x,g,bin=FALSE){a<-x[g==1];b<-x[g==0]; if(bin){p1<-mean(a,na.rm=T);p0<-mean(b,na.rm=T); (p1-p0)/sqrt((p1*(1-p1)+p0*(1-p0))/2)} else (mean(a,na.rm=T)-mean(b,na.rm=T))/sqrt((var(a,na.rm=T)+var(b,na.rm=T))/2)}
A[,eicu:=as.integer(db=="eICU-CRD")]
V<-list(c("Data source, eICU-CRD","eicu","b"),c("Age, years","age","c0"),c("Male sex","male","b"),c("Do-not-intubate order","dni","b"),c("COPD","copd","b"),c("Obesity or hypoventilation","obes","b"),c("Congestive heart failure","chf","b"),c("Chronic kidney disease","ckd","b"),c("Diabetes mellitus","dm","b"),c("Pneumonia","pna","b"),c("Sepsis","sepsis","b"),
 c("pH","ph","c2"),c("PaCO2, mmHg","paco2","c0"),c("Bicarbonate, mmol/L","hco3","c1"),c("Expected bicarbonate, mmol/L","exp_hco3","c1"),c("ΔHCO3, mmol/L","dh","c1"),c("Lactate, mmol/L","lact","c1"),c("Creatinine, mg/dL","creat","c1"),c("Albumin, g/dL","alb","c1"),c("Hemoglobin, g/dL","hb","c1"),c("White cell count, ×10⁹/L","wbc","c1"),
 c("Heart rate, beats/min","hr","c0"),c("Respiratory rate, breaths/min","rr","c0"),c("Mean arterial pressure, mmHg","map","c0"),c("SpO2, %","spo2","c0"),c("Glasgow Coma Scale","gcs","c0"),c("Vasopressor at NIV start","vaso","b"),c("Systemic corticosteroid","ster","b"),c("Loop diuretic in prior 24 h","loop","b"),c("Sedative or opioid in first 6 h","sedv","b"),c("APACHE IV score (eICU-CRD only)","sev","c0"),
 c("NIV failure within 7 days","fail7","b"),c("Intubation within 7 days","intub7","b"),c("Death within 7 days","death7","b"),c("Death within 28 days","death28","b"),c("In-hospital death","hdeath","b"))
T1<-rbindlist(lapply(V,function(v){x<-A[[v[2]]]; g<-A$hi; b<-v[3]=="b"; dg<-if(b) 0 else as.integer(substr(v[3],2,2))
 fx<-function(i) if(b) fb(x[i]) else f1(x[i],dg); miss<-mean(is.na(x)); data.table(Variable=v[1],Overall=fx(rep(TRUE,length(x))),`Q1–Q3`=fx(g==0),`Q4`=fx(g==1),SMD=sprintf("%.2f",abs(smd(x,g,b))),Missing=if(miss>0) sprintf("%.1f",100*miss) else "0")}))
saveRDS(T1,"T1.rds"); print(T1,nrows=60)
TS3<-rbindlist(lapply(V[-1],function(v){x<-A[[v[2]]]; g<-A$eicu; b<-v[3]=="b"; dg<-if(b) 0 else as.integer(substr(v[3],2,2))
 fx<-function(i) if(b) fb(x[i]) else f1(x[i],dg); ms<-function(i) if(mean(is.na(x[i]))>0) sprintf("%.1f",100*mean(is.na(x[i]))) else "0"
 data.table(Variable=v[1],`MIMIC-IV`=fx(g==0),`eICU-CRD`=fx(g==1),SMD=if(v[2]=="sev") "–" else sprintf("%.2f",abs(smd(x,g,b))),`Missing MIMIC-IV, %`=ms(g==0),`Missing eICU-CRD, %`=ms(g==1))}))
TS3[Variable=="APACHE IV score (eICU-CRD only)",`:=`(`MIMIC-IV`="–",`Missing MIMIC-IV, %`="–")]; saveRDS(TS3,"TS3.rds"); print(TS3,nrows=60)
# ---------------- Figures
TH<-theme_classic(base_size=9,base_family="Helvetica")+theme(strip.background=element_blank(),strip.text=element_text(face="bold",size=9),legend.position="bottom",legend.title=element_text(size=8),legend.text=element_text(size=8),axis.text=element_text(colour="black"))
sv<-function(p,f,w,h){ggsave(paste0("fig/",f,".tiff"),p,width=w,height=h,dpi=600,compression="lzw",bg="white"); ggsave(paste0("fig/",f,".png"),p,width=w,height=h,dpi=300,bg="white")}
# Fig1 flow
fl<-R$flow; nf<-R$nflow
bx<-data.table(x=c(1,1,1,1,3,3,3,3,2),y=c(4,3,2,1,4,3,2,1,0),
 t=c(sprintf("MIMIC-IV v3.1\nFirst ICU stay with NIV, adults,\nno invasive ventilation before NIV\n(n = %s)",format(fl$m_niv,big.mark=",")),sprintf("Arterial PaCO2 within −6 to +2 h\nof NIV start (n = %s)",format(nf$m[["abg"]],big.mark=",")),sprintf("PaCO2 ≥45 mmHg with pH\nand bicarbonate available (n = %d)",nf$m[["hyp"]]),sprintf("pH <7.35\nMIMIC-IV cohort (n = %d)",nf$m[["acid"]]),
  sprintf("eICU-CRD v2.0\nFirst ICU stay with NIV, adults,\nno invasive ventilation before NIV\n(n = %s)",format(fl$e_niv,big.mark=",")),sprintf("Arterial PaCO2 within −6 to +2 h\nof NIV start (n = %s)",format(nf$e[["abg"]],big.mark=",")),sprintf("PaCO2 ≥45 mmHg with pH\nand bicarbonate available (n = %s)",format(nf$e[["hyp"]],big.mark=",")),sprintf("pH <7.35\neICU-CRD cohort (n = %d)",nf$e[["acid"]]),
  sprintf("Pooled analytic cohort (n = %s)\nΔHCO3 quartiles within each database\nQ1–Q3 (n = %s) vs Q4, largest deficit (n = %d)",format(nf$m[["acid"]]+nf$e[["acid"]],big.mark=","),format(X$Pooled$n[6],big.mark=","),X$Pooled$n[3])))
W<-1.55; H<-0.34; bx[,h:=ifelse(y==4,0.42,H)]
ex<-data.table(x=c(1,1,1,3,3,3),y=c(3.5,2.5,1.5,3.5,2.5,1.5),t=c(sprintf("No arterial PaCO2 in window: %s",format(fl$m_niv-nf$m[["abg"]],big.mark=",")),sprintf("PaCO2 <45 mmHg or missing pH/HCO3: %d",nf$m[["abg"]]-nf$m[["hyp"]]),sprintf("pH ≥7.35: %d",nf$m[["hyp"]]-nf$m[["acid"]]),
 sprintf("No arterial PaCO2 in window: %s",format(fl$e_niv-nf$e[["abg"]],big.mark=",")),sprintf("PaCO2 <45 mmHg or missing pH/HCO3: %s",format(nf$e[["abg"]]-nf$e[["hyp"]],big.mark=",")),sprintf("pH ≥7.35: %d",nf$e[["hyp"]]-nf$e[["acid"]])))
p<-ggplot()+geom_rect(data=bx,aes(xmin=x-W/2,xmax=x+W/2,ymin=y-h,ymax=y+h),fill="white",colour="black",linewidth=.35)+geom_text(data=bx,aes(x,y,label=t),size=2.35,lineheight=.95,family="Helvetica")+
 geom_segment(data=bx[y%in%c(4,3,2)],aes(x=x,xend=x,y=y-h,yend=y-1+H+0.05),arrow=arrow(length=unit(1.3,"mm"),type="closed"),linewidth=.35)+
 geom_segment(data=data.table(x=c(1,3)),aes(x=x,xend=ifelse(x<2,2-0.06,2+0.06),y=1-H,yend=0+H+0.05),arrow=arrow(length=unit(1.3,"mm"),type="closed"),linewidth=.35)+
 geom_label(data=ex,aes(x=x+0.03,y=y,label=t),hjust=0,size=2.1,linewidth=0.25,label.r=unit(0,"pt"),family="Helvetica",fill="white")+
 coord_cartesian(xlim=c(0.15,3.95),ylim=c(-0.45,4.5))+theme_void()
sv(p,"Figure_1",7,5.6)
# Fig2 KM
km<-function(d,lab){s<-survfit(Surv(ft/24,fail7)~Q,d); data.table(t=c(0,s$time),y=c(0,1-s$surv),g=c(NA,rep(paste0("Q",1:4),s$strata)),db=lab)[,g:=fifelse(is.na(g),"Q1",g)]}
mk<-function(d,lab){s<-survfit(Surv(ft/24,fail7)~Q,d); r<-rbindlist(lapply(1:4,function(qq){x<-d[as.character(Q)==as.character(qq)]; s2<-survfit(Surv(ft/24,fail7)~1,x); data.table(t=c(0,s2$time),y=c(0,1-s2$surv),g=paste0("Q",qq))})); r[,db:=lab]; r}
fd<-rbind(mk(A,"Pooled cohort"),mk(A[db=="MIMIC-IV"],"MIMIC-IV"),mk(A[db=="eICU-CRD"],"eICU-CRD")); fd[,db:=factor(db,levels=c("Pooled cohort","MIMIC-IV","eICU-CRD"))]
lp<-data.table(db=factor(c("Pooled cohort","MIMIC-IV","eICU-CRD"),levels=levels(fd$db)),lab=c("Log-rank P < 0.001","Log-rank P < 0.001","Log-rank P < 0.001"))
p<-ggplot(fd,aes(t,y,colour=g))+geom_step(linewidth=.6)+facet_wrap(~db)+geom_text(data=lp,aes(x=0.2,y=0.76,label=lab),inherit.aes=FALSE,hjust=0,size=2.8)+
 scale_y_continuous(labels=function(x) paste0(100*x,"%"),limits=c(0,0.8))+scale_x_continuous(breaks=0:7)+scale_colour_manual(values=c("#1b9e77","#7570b3","#e6ab02","#d95f02"),name="ΔHCO3 quartile (Q4 = largest deficit)")+labs(x="Days after NIV initiation",y="Cumulative incidence of NIV failure")+TH
sv(p,"Figure_2",7.2,3.3)
# Fig3 RCS
rd<-X$rcsdat; p<-ggplot(rd,aes(dh,or))+geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey80")+geom_line(linewidth=.7)+geom_hline(yintercept=1,linetype="dashed",colour="grey40",linewidth=.4)+
 scale_x_reverse()+scale_y_log10()+annotate("text",x=max(rd$dh),y=max(rd$hi)*0.85,hjust=0,size=2.8,label=sprintf("P overall < 0.001\nReference: ΔHCO3 = %.1f mmol/L (median)",X$ref))+
 labs(x="ΔHCO3, observed minus expected bicarbonate (mmol/L)\n← smaller bicarbonate deficit        larger bicarbonate deficit →",y="Adjusted odds ratio for NIV failure (log scale)")+TH
p<-p+scale_x_reverse()
sv(p,"Figure_3",4.6,3.8)
# Fig4 forest
fo<-function(lab,x,grp) data.table(lab=lab,or=x[1],lo=x[2],hi=x[3],p=x[4],grp=grp)
S<-R$Pooled$sens; gs<-function(k) {r<-S[a==k]; v<-as.numeric(strsplit(gsub("[()–]"," ",r$or)," +")[[1]]); c(v[1:3],r$p)}
F<-rbind(fo("Pooled, Model 2",X$Pooled$zM2,"Main estimate"),fo("MIMIC-IV, Model 2",X$`MIMIC-IV`$zM2,"Database-specific"),fo("eICU-CRD, Model 2",X$`eICU-CRD`$zM2,"Database-specific"),
 fo("Pooled, Model 3 (multiple imputation)",MI$m3_z,"Alternative adjustment"),fo("Hospital-clustered robust SE",gs("Sai số chuẩn robust cụm theo bệnh viện"),"Alternative adjustment"),fo("Pooled, Model 3 (median imputation)",X$Pooled$zM3,"Alternative adjustment"),fo("Additional adjustment for lactate",gs("Hiệu chỉnh thêm lactate"),"Alternative adjustment"),
 fo("Excluding do-not-intubate orders",gs("Loại bệnh nhân DNI"),"Restriction"),fo("Excluding prior loop diuretic",gs("Loại bệnh nhân dùng lợi tiểu quai 24 h trước"),"Restriction"),fo("Excluding failure within 6 h",X$excl_early,"Restriction"),
 fo("COPD only",X$copd,"Restriction"),fo("Arterial blood gas HCO3 only",gs("Chỉ HCO3 khí máu động mạch"),"Exposure definition"),fo("Serum bicarbonate only",gs("Chỉ HCO3 huyết thanh"),"Exposure definition"),
 fo("Expected slope 0.40 mmol/L per mmHg",gs("Công thức kỳ vọng 0,4×ΔPaCO2"),"Exposure definition"),fo("Expected slope 0.48 mmol/L per mmHg",X$s048,"Exposure definition"))
saveRDS(F,"F4.rds"); F[,lab:=factor(lab,levels=rev(lab))]; F[,txt:=sprintf("%.2f (%.2f–%.2f)",or,lo,hi)]; F[,grp:=factor(grp,levels=unique(grp))]
p<-ggplot(F,aes(or,lab))+geom_vline(xintercept=1,linetype="dashed",colour="grey40",linewidth=.4)+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.25,linewidth=.45)+geom_point(shape=15,size=1.9)+
 geom_text(aes(x=2.45,label=txt),hjust=0,size=2.6)+facet_grid(grp~.,scales="free_y",space="free_y",switch="y")+scale_x_log10(breaks=c(0.8,1,1.25,1.5,2),limits=c(0.8,3.4))+
 labs(x="Adjusted odds ratio per 1-SD lower ΔHCO3 (log scale)",y=NULL)+TH+theme(strip.placement="outside",strip.text.y.left=element_text(angle=0,hjust=1,face="bold",size=7.5),panel.spacing=unit(2,"pt"))
sv(p,"Figure_4",7,4.6)
# FigS1 ROC
r1<-X$roc$r1; r2<-X$roc$r2; rr<-rbind(data.table(fp=1-r1$specificities,tp=r1$sensitivities,m=sprintf("ΔHCO3, AUC %.3f",auc(r1))),data.table(fp=1-r2$specificities,tp=r2$sensitivities,m=sprintf("pH, AUC %.3f",auc(r2))))
p<-ggplot(rr[order(m,fp,tp)],aes(fp,tp,colour=m))+geom_path(linewidth=.6)+geom_abline(linetype="dashed",colour="grey50",linewidth=.4)+annotate("text",x=.55,y=.12,label=sprintf("DeLong P < 0.001"),size=2.8)+scale_colour_manual(values=c("#d95f02","#1b9e77"),name=NULL)+labs(x="1 − specificity",y="Sensitivity")+coord_equal()+TH
sv(p,"Figure_S1",4,4.2)
cat("done\n")
