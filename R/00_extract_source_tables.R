## Build source tables for the HCO3-NIV study directly from the PhysioNet CSV files.
## Set MIMIC_DIR (MIMIC-IV v3.1 root containing hosp/ and icu/) and EICU_DIR (eICU-CRD v2.0 root).
suppressMessages({library(DBI);library(duckdb)})
MIMIC_DIR<-Sys.getenv("MIMIC_DIR"); EICU_DIR<-Sys.getenv("EICU_DIR"); OUT<-Sys.getenv("SOURCE_DB","hco3_source.duckdb")
stopifnot(dir.exists(MIMIC_DIR),dir.exists(EICU_DIR))
con<-dbConnect(duckdb::duckdb(),OUT); X<-function(s) invisible(dbExecute(con,s)); N<-function(t) dbGetQuery(con,paste("select count(*) n from",t))$n
X("SET threads=8"); X("SET memory_limit='16GB'")
m<-function(p) sprintf("read_csv('%s/%s.csv.gz',header=true,auto_detect=true,sample_size=-1,ignore_errors=true)",MIMIC_DIR,p)
e<-function(p) sprintf("read_csv('%s/%s.csv.gz',header=true,all_varchar=true,ignore_errors=true)",EICU_DIR,p)
log<-function(...) cat(format(Sys.time(),"%H:%M:%S"),...,"\n")
## ---------------- MIMIC-IV ----------------
X(paste("create or replace table s0 as select i.*, p.gender, p.anchor_age+(year(i.intime)-p.anchor_year) age, p.dod, a.admittime,a.dischtime,a.deathtime,a.race,a.insurance,a.admission_type,a.hospital_expire_flag from",m("icu/icustays"),"i join",m("hosp/patients"),"p using(subject_id) join",m("hosp/admissions"),"a using(hadm_id)"))
X(paste("create or replace table dev as select stay_id,itemid,charttime,value from",m("icu/chartevents"),"where itemid in (226732,223758)")); log("dev",N("dev"))
X(paste("create or replace table imv as select stay_id,starttime vs,endtime ve from",m("icu/procedureevents"),"where itemid=225792"))
X("create or replace table mniv as select stay_id,min(charttime) t0 from dev where itemid=226732 and (value like 'Bipap%' or value like 'CPAP%') group by 1")
X("create or replace table mc as with a as (select s.*,n.t0,row_number() over(partition by s.subject_id order by n.t0) rn from s0 s join mniv n using(stay_id)
 where n.t0 between s.intime - interval 2 hour and s.outtime and s.age>=18
 and s.stay_id in (select stay_id from (select stay_id,row_number() over(partition by subject_id order by intime) k from s0) where k=1))
 select a.* exclude(rn) from a where rn=1
 and not exists(select 1 from imv i where i.stay_id=a.stay_id and i.vs<=a.t0)
 and not exists(select 1 from dev v where v.stay_id=a.stay_id and v.itemid=226732 and (v.value like 'Endotracheal%' or v.value like 'Tracheostomy%') and v.charttime<=a.t0)"); log("mc",N("mc"))
X(paste("create or replace table mlab as select l.hadm_id,l.specimen_id,l.itemid,l.charttime,l.valuenum,l.value,m.stay_id,m.t0 from",m("hosp/labevents"),"l join mc m using(hadm_id)
 where l.itemid in (50818,50820,50821,50803,50813,50912,51006,51301,51222,51265,50882,50862,50983,50902,50971,50931,50885) and l.charttime between m.t0 - interval 24 hour and m.t0 + interval 48 hour")); log("mlab",N("mlab"))
X(paste("create or replace table mce as select c.stay_id,c.itemid,c.charttime,c.valuenum,c.value,m.t0 from",m("icu/chartevents"),"c join mc m using(stay_id)
 where (c.itemid in (226512,224639,220045,220210,224690,220052,220181,220277,220739,223900,223901,223835,226732) and c.charttime between m.t0 - interval 24 hour and m.t0 + interval 48 hour) or c.itemid=223758")); log("mce",N("mce"))
X(paste("create or replace table min_ as select i.stay_id,i.itemid,i.starttime,i.endtime,m.t0 from",m("icu/inputevents"),"i join mc m using(stay_id) where i.endtime>=m.t0 - interval 24 hour and i.starttime<=m.t0 + interval 72 hour"))
X(paste("create or replace table mrx as select e.hadm_id,lower(e.medication) med,e.charttime,m.stay_id,m.t0 from",m("hosp/emar"),"e join mc m using(hadm_id)
 where e.event_txt in ('Administered','Started','Confirmed','Restarted') and e.charttime between m.t0 - interval 24 hour and m.t0 + interval 48 hour
 and regexp_matches(lower(e.medication),'furosemide|bumetanide|torsemide|methylpred|prednison|hydrocort|dexameth|midazolam|lorazepam|morphine|fentanyl|dexmedetom|propofol|haloperidol|quetiapine')"))
X("create or replace table mimv as select i.stay_id,i.vs,i.ve from imv i join mc using(stay_id)")
X("create or replace table mett as select v.stay_id,min(v.charttime) t from dev v join mc m using(stay_id) where v.itemid=226732 and (v.value like 'Endotracheal%' or v.value like 'Tracheostomy%') and v.charttime>m.t0 group by 1")
X(paste("create or replace table mdx as select d.hadm_id,trim(d.icd_code) icd_code,d.icd_version from",m("hosp/diagnoses_icd"),"d join mc using(hadm_id)")); log("MIMIC done")
## ---------------- eICU-CRD ----------------
NIV<-"regexp_matches(lower(v),'bipap|cpap|nippv|niv|non-invasive|noninvasive|bi-pap|avaps')"
X(paste("create or replace table rcd as select cast(patientunitstayid as bigint) s,cast(respchartoffset as int) o,respchartvalue v from",e("respiratoryCharting"),"where respchartvaluelabel in ('O2 Device','Oxygen Delivery Method')"))
X(paste("create or replace table ncd as select cast(patientunitstayid as bigint) s,cast(nursingchartoffset as int) o,nursingchartvalue v from",e("nurseCharting"),"where nursingchartcelltypevallabel='O2 Admin Device'"))
X(paste("create or replace table eniv as select s,min(o) t0 from (select * from rcd union all select * from ncd) where",NIV,"group by 1")); log("eniv",N("eniv"))
X(paste("create or replace table ep as select distinct cast(patientunitstayid as bigint) s,try_cast(ventstartoffset as int) vs from",e("respiratoryCare"),"where try_cast(ventstartoffset as int) is not null and try_cast(ventstartoffset as int)<>0"))
X(paste("create or replace table ec as with p as (select cast(patientunitstayid as bigint) s,uniquepid,try_cast(hospitalid as int) hosp,try_cast(replace(age,'> 89','90') as double) age,gender,
 try_cast(admissionweight as double) w,try_cast(admissionheight as double) h,try_cast(hospitaldischargeoffset as int) ho,hospitaldischargestatus hs,try_cast(unitdischargeoffset as int) uo,
 try_cast(hospitaladmitoffset as int) hao,try_cast(unitvisitnumber as int) uvn from",e("patient"),"),
 a as (select p.*,n.t0,row_number() over(partition by uniquepid order by p.hao desc,p.uvn) rn from p join eniv n using(s) where n.t0 between -120 and p.uo and p.age>=18)
 select a.* exclude(rn) from a where rn=1 and not exists(select 1 from ep x where x.s=a.s and x.vs<=a.t0)")); log("ec",N("ec"))
X(paste("create or replace table elab as select cast(l.patientunitstayid as bigint) s,cast(l.labresultoffset as int) o,l.labname n,try_cast(l.labresult as double) v,c.t0 from",e("lab"),"l join ec c on c.s=cast(l.patientunitstayid as bigint)
 where l.labname in ('pH','paCO2','paO2','lactate','creatinine','BUN','WBC x 1000','Hgb','bicarbonate','HCO3','albumin','sodium','chloride','potassium','glucose') and cast(l.labresultoffset as int) between c.t0-1440 and c.t0+2880"))
X(paste("create or replace table enc as select cast(n.patientunitstayid as bigint) s,cast(n.nursingchartoffset as int) o,n.nursingchartcelltypevallabel lab,n.nursingchartcelltypevalname nm,n.nursingchartvalue v,c.t0 from",e("nurseCharting"),"n join ec c on c.s=cast(n.patientunitstayid as bigint)
 where cast(n.nursingchartoffset as int) between c.t0-1440 and c.t0+2880")); log("enc",N("enc"))
X(paste("create or replace table einf as select cast(i.patientunitstayid as bigint) s,cast(i.infusionoffset as int) o,lower(i.drugname) d,c.t0 from",e("infusionDrug"),"i join ec c on c.s=cast(i.patientunitstayid as bigint) where cast(i.infusionoffset as int) between c.t0-1440 and c.t0+2880"))
X(paste("create or replace table emed as select cast(m.patientunitstayid as bigint) s,cast(m.drugstartoffset as int) o,lower(m.drugname) d,c.t0 from",e("medication"),"m join ec c on c.s=cast(m.patientunitstayid as bigint)
 where cast(m.drugstartoffset as int) between c.t0-1440 and c.t0+2880 and regexp_matches(lower(m.drugname),'furosemide|lasix|bumetanide|bumex|torsemide|methylpred|solu-medrol|prednison|hydrocort|dexameth|midazolam|lorazepam|morphine|fentanyl|dexmedetom|precedex|propofol|haloperidol')"))
X("create or replace table eep as select x.s,x.vs from ep x join ec using(s)")
X(paste("create or replace table edx as select cast(d.patientunitstayid as bigint) s,lower(d.diagnosisstring) ds from",e("diagnosis"),"d join ec c on c.s=cast(d.patientunitstayid as bigint)"))
X(paste("create or replace table eap as select cast(a.patientunitstayid as bigint) s,max(try_cast(a.apachescore as double)) apache from",e("apachePatientResult"),"a join ec c on c.s=cast(a.patientunitstayid as bigint) where a.apacheversion='IVa' group by 1"))
X(paste("create or replace table ecp as select cast(g.patientunitstayid as bigint) s,cast(g.cplitemoffset as int) o,g.cplgroup,g.cplitemvalue from",e("carePlanGeneral"),"g join ec c on c.s=cast(g.patientunitstayid as bigint)"))
log("eICU done"); dbDisconnect(con,shutdown=TRUE)
