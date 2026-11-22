suppressMessages({library(DBI);library(data.table)})
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),Sys.getenv("SOURCE_DB","hco3_source.duckdb"),read_only=TRUE); Q<-function(s) as.data.table(dbGetQuery(con,s))
H<-function(a,b) as.numeric(difftime(a,b,units="hours"))
# ================= MIMIC =================
m<-Q("select * from mc"); fl<-list(m_niv=nrow(m))
lb<-Q("select stay_id,specimen_id,itemid,charttime,valuenum v,t0 from mlab where valuenum is not null"); lb[,h:=H(charttime,t0)]
cat("50803 n:",lb[itemid==50803,.N]," 50882 n:",lb[itemid==50882,.N],"\n")
# baseline ABG: PaCO2 sample closest to t0 within -6..+2h, HCO3 from same specimen (50803) else serum 50882 within 2h
pc<-lb[itemid==50818&h>=-6&h<=2&v%between%c(10,200)][,d:=abs(h)][order(d)][,.SD[1],by=stay_id]
ph<-lb[itemid==50820,.(stay_id,specimen_id,ph=v)]; hc<-lb[itemid==50803,.(stay_id,specimen_id,hco3a=v)]
pc<-merge(merge(pc,ph,by=c("stay_id","specimen_id"),all.x=TRUE),hc,by=c("stay_id","specimen_id"),all.x=TRUE)
sh<-merge(pc[,.(stay_id,hp=h)],lb[itemid==50882,.(stay_id,h,hs=v)],by="stay_id")[abs(h-hp)<=4][,d:=abs(h-hp)][order(d)][,.SD[1],by=stay_id]
pc[,hco3s:=sh$hs[match(stay_id,sh$stay_id)]]; pc[,hco3:=fifelse(!is.na(hco3a),hco3a,hco3s)]; pc[,src:=fifelse(!is.na(hco3a),"ABG","serum")]
m<-merge(m,pc[,.(stay_id,paco2=v,ph,hco3,src,hco3s,abg_h=h)],by="stay_id")
near<-function(it,lo,hi,a=-24,z=2){x<-lb[itemid==it & v%between%c(lo,hi) & h>=a & h<=z][,d:=abs(h)][order(d)][,.SD[1],by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(lact=near(50813,0.1,30),creat=near(50912,0.1,20),bun=near(51006,1,250),alb=near(50862,1,6),hb=near(51222,3,20),wbc=near(51301,0.1,200),na=near(50983,100,180),cl=near(50902,60,140),k=near(50971,1.5,9),glu=near(50931,20,1500),pao2=near(50821,20,600,-6,2))]
ce<-Q("select stay_id,itemid,charttime,valuenum v,value,t0 from mce"); ce[,h:=H(charttime,t0)]
nc<-function(its,lo,hi,a=-6,z=2){x<-ce[itemid %in% its & v%between%c(lo,hi) & h>=a & h<=z][,d:=abs(h)][order(d)][,.SD[1],by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(hr=nc(220045,20,250),rr=nc(c(220210,224690),3,70),map=nc(c(220052,220181),30,200),spo2=nc(220277,50,100),fio2=nc(223835,21,100))]
g<-ce[itemid %in% c(220739,223900,223901)&h>=-6&h<=2][,d:=abs(h)][order(d)][,.SD[1],by=.(stay_id,itemid)][,.(gcs=if(.N==3) sum(v) else NA_real_),by=stay_id]; m[,gcs:=g$gcs[match(stay_id,g$stay_id)]]
wt<-ce[itemid %in% c(226512,224639)&v%between%c(30,300),.(w=median(v)),by=stay_id]; m[,wt:=wt$w[match(stay_id,wt$stay_id)]]
m[,dni:=as.integer(stay_id %in% ce[itemid==223758 & grepl("DNI|DNR / DNI|Comfort",value,ignore.case=TRUE)]$stay_id)]
m[,dni_pre:=as.integer(stay_id %in% ce[itemid==223758 & h<=0 & grepl("DNI|DNR / DNI|Comfort",value,ignore.case=TRUE)]$stay_id)]
dx<-Q("select hadm_id,icd_code from mdx"); fx<-function(p) as.integer(m$hadm_id %in% dx[grepl(p,icd_code)]$hadm_id)
m[,`:=`(copd=fx("^(491|492|496|J43|J44)"),obes=fx("^(2780|32723|E66|G4733)"),chf=fx("^(428|I50)"),ckd=fx("^(585|N18)"),pna=fx("^(48[0-6]|J1[2-8]|J69|507)"),sepsis=fx("^(99591|99592|78552|A4[01]|R652)"),dm=fx("^(250|E1[01])"))]
inp<-Q("select stay_id,itemid,starttime,endtime,t0 from min_"); inp[,`:=`(hs=H(starttime,t0),he=H(endtime,t0))]
m[,vaso:=as.integer(stay_id %in% inp[itemid %in% c(221906,221289,221662,222315,221749)&hs<=2&he>=-6]$stay_id)]
rx<-Q("select stay_id,med,charttime,t0 from mrx"); rx[,h:=H(charttime,t0)]
m[,`:=`(ster=as.integer(stay_id %in% rx[grepl("methylpred|prednis|hydrocort|dexameth",med)&h<=0]$stay_id),loop=as.integer(stay_id %in% rx[grepl("furosemide|bumetanide|torsemide",med)&h<=0&h>=-24]$stay_id),
        sedv=as.integer(stay_id %in% rx[grepl("midazolam|lorazepam|morphine|fentanyl|dexmedetom|propofol|haloperidol|quetiapine",med)&h<=6&h>=0]$stay_id))]
imv<-Q("select * from mimv"); ett<-Q("select * from mett")
it<-rbind(imv[,.(stay_id,t=vs)],ett[,.(stay_id,t)]); it<-merge(it,m[,.(stay_id,t0)],by="stay_id")[t>t0,.(t=min(t)),by=stay_id]
m[,intub_h:=H(it$t[match(stay_id,it$stay_id)],t0)]
dth<-as.POSIXct(m$deathtime,tz="UTC"); dd<-as.POSIXct(ifelse(is.na(m$dod),NA,paste(as.character(m$dod),"12:00:00")),tz="UTC"); dth[is.na(dth)]<-dd[is.na(dth)]
m[,death_h:=H(dth,t0)]; m[,hdeath:=hospital_expire_flag]; m[,icuout_h:=H(outtime,t0)]; m[,disch_h:=H(dischtime,t0)]; m[,male:=as.integer(gender=="M")]; m[,db:="MIMIC-IV"]; m[,sev:=NA_real_]; m[,hosp:=1L]
# ================= eICU =================
e<-Q("select * from ec"); fl$e_niv<-nrow(e)
el<-Q("select s,o,n,v,t0 from elab where v is not null"); el[,h:=(o-t0)/60]
pce<-el[n=="paCO2"&h>=-6&h<=2&v%between%c(10,200)][,d:=abs(h)][order(d)][,.SD[1],by=s]
pr<-function(nm,lo,hi,w=1){x<-merge(pce[,.(s,hp=h)],el[n %in% nm & v%between%c(lo,hi),.(s,h,v)],by="s")[abs(h-hp)<=w][,d:=abs(h-hp)][order(d)][,.SD[1],by=s]; x$v[match(pce$s,x$s)]}
pce[,`:=`(ph=pr("pH",6.5,7.8,0.25),hco3a=pr("HCO3",3,60,0.25),hco3s=pr("bicarbonate",3,60,4))]; pce[,hco3:=fifelse(!is.na(hco3a),hco3a,hco3s)]; pce[,src:=fifelse(!is.na(hco3a),"ABG","serum")]
e<-merge(e,pce[,.(s,paco2=v,ph,hco3,src,hco3s,abg_h=h)],by="s")
ne<-function(nm,lo,hi,a=-24,z=2){x<-el[n %in% nm & v%between%c(lo,hi) & h>=a & h<=z][,d:=abs(h)][order(d)][,.SD[1],by=s]; x$v[match(e$s,x$s)]}
e[,`:=`(lact=ne("lactate",0.1,30),creat=ne("creatinine",0.1,20),bun=ne("BUN",1,250),alb=ne("albumin",1,6),hb=ne("Hgb",3,20),wbc=ne("WBC x 1000",0.1,200),na=ne("sodium",100,180),cl=ne("chloride",60,140),k=ne("potassium",1.5,9),glu=ne("glucose",20,1500),pao2=ne("paO2",20,600,-6,2))]
nc_<-Q("select s,o,lab,nm,try_cast(v as double) v,t0 from enc"); nc_[,h:=(o-t0)/60]
cat("enc names:",paste(head(nc_[,.N,by=nm][order(-N)]$nm,25),collapse="|"),"\n")
nn<-function(pat,lo,hi){x<-nc_[grepl(pat,nm,ignore.case=TRUE)&v%between%c(lo,hi)&h>=-6&h<=2][,d:=abs(h)][order(d)][,.SD[1],by=s]; x$v[match(e$s,x$s)]}
e[,`:=`(hr=nn("^Heart Rate$",20,250),rr=nn("^Respiratory Rate$",3,70),map=nn("MAP|Mean",30,200),spo2=nn("O2 Saturation|SpO2",50,100),gcs=nn("^GCS Total$",3,15),fio2=NA_real_)]
ap<-Q("select s,max(apache) a from eap group by 1"); e[,sev:=ap$a[match(s,ap$s)]]; e[sev<0,sev:=NA]
edx<-Q("select s,lower(ds) ds from edx"); fe<-function(p) as.integer(e$s %in% edx[grepl(p,ds)]$s)
e[,`:=`(copd=fe("copd|chronic obstructive|emphysema|chronic bronchitis"),obes=fe("obes|hypoventilation|sleep apnea"),chf=fe("congestive heart failure|chf|heart failure"),ckd=fe("chronic kidney|ckd|chronic renal|esrd|end stage renal"),pna=fe("pneumonia"),sepsis=fe("sepsis|septic"),dm=fe("diabetes"))]
cp<-Q("select c.s,c.o,lower(c.cplitemvalue) v,e.t0 from ecp c join ec e on e.s=c.s"); e[,dni:=as.integer(s %in% cp[grepl("no intubation|comfort|do not intubate",v)]$s)]; e[,dni_pre:=as.integer(s %in% cp[grepl("no intubation|comfort|do not intubate",v) & o<=t0]$s)]
inf<-Q("select s,o,d,t0 from einf"); inf[,h:=(o-t0)/60]; e[,vaso:=as.integer(s %in% inf[grepl("norepi|levophed|epineph|dopamine|vasopressin|phenyleph",d)&h<=2&h>=-6]$s)]
md<-Q("select s,o,d,t0 from emed"); md[,h:=(o-t0)/60]
e[,`:=`(ster=as.integer(s %in% md[grepl("methylpred|solu-medrol|prednis|hydrocort|dexameth",d)&h<=0]$s),loop=as.integer(s %in% md[grepl("furosemide|lasix|bumetanide|bumex|torsemide",d)&h<=0]$s),
        sedv=as.integer(s %in% c(md[grepl("midazolam|lorazepam|morphine|fentanyl|dexmedetom|precedex|propofol|haloperidol",d)&h>=0&h<=6]$s,inf[grepl("midazolam|lorazepam|morphine|fentanyl|dexmedetom|precedex|propofol",d)&h>=0&h<=6]$s)))]
ep<-Q("select * from eep"); ep<-merge(ep,e[,.(s,t0)],by="s")[vs>t0,.(t=min(vs)),by=s]; e[,intub_h:=(ep$t[match(s,ep$s)]-t0)/60]
e[,hdeath:=as.integer(hs=="Expired")]; e[,disch_h:=(ho-t0)/60]; e[,death_h:=fifelse(hdeath==1,disch_h,NA_real_)]; e[,icuout_h:=(uo-t0)/60]
e[,`:=`(male=as.integer(gender=="Male"),age=as.numeric(ifelse(age=="> 89","90",age)),wt=w,db="eICU-CRD",stay_id=s)]
setnames(m,"t0","t0x",skip_absent=TRUE)
cols<-c("stay_id","db","hosp","abg_h","dni_pre","age","male","wt","paco2","ph","hco3","src","hco3s","lact","creat","bun","alb","hb","wbc","na","cl","k","glu","pao2","hr","rr","map","spo2","gcs","fio2","sev","dni","copd","obes","chf","ckd","pna","sepsis","dm","vaso","ster","loop","sedv","intub_h","death_h","hdeath","icuout_h","disch_h")
D<-rbind(m[,..cols],e[,hosp:=hosp][,..cols]); fl$m_abg<-nrow(m); fl$e_abg<-nrow(e)
saveRDS(list(d=D,flow=fl),"D.rds"); print(fl); print(D[,.(n=.N,copd=mean(copd),hco3=mean(!is.na(hco3)),abg=mean(src=="ABG",na.rm=T),intub=mean(intub_h<=168,na.rm=T)),by=db])
