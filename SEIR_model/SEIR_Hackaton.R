library(EpiModel)
library(dplyr)
library(tidyr)
library(data.table)
library(httr)
library(RCurl)
library(foreign)
library(lubridate)

# retrieve data and data prep

list_files_remote_git<-function(git_master, git_directory){
  req <- GET(paste0(git_master,"?recursive=1"))
  stop_for_status(req)
  filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)
  ll<-grep(git_directory, filelist, value = TRUE, fixed = TRUE)
  return(ll[grepl(".csv",ll)])
}
git_raw="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/"
git_api="https://api.github.com/repos/CSSEGISandData/COVID-19/git/trees/master"
git_directory="csse_covid_19_data/csse_covid_19_daily_reports"
ls=list_files_remote_git(git_api,git_directory)

get_data_from_dir_git <- function(git_raw,git_api,git_directory) {
  #listfiles 
  ll<-list_files_remote_git(git_api,git_directory)
  url_list=lapply(ll,function(x){paste0(git_raw,x)})
  
  # read the files into a list of data.frames
  data.list <- lapply(url_list, function(x){read.csv(text=getURL(x),header = T)})
  good_db<-sapply(data.list,function(x)grepl("FIPS",paste(names(x),collapse = " ")))
  
  data.list<-data.list[good_db]
  # concatenate into one big data.frame
  data.cat <- do.call(rbind, data.list)
  
  # aggregate
  # data.agg <- aggregate(value ~ index, data.cat, mean)
}

data_county<-get_data_from_dir_git(git_raw,git_api,git_directory)

data_county_us<-data_county[which(data_county$Country_Region=="US"),]
data_county_us$time_f1<-mdy_hm(data_county_us$Last_Update,tz="EST")
data_county_us$time_f2<-ymd_hms(data_county_us$Last_Update,tz="EST")
dt_tmp1<-data_county_us[which(!is.na(data_county_us$time_f2)),]
dt_tmp1$time<-dt_tmp1$time_f2
dt_tmp2<-data_county_us[which(is.na(data_county_us$time_f2)),]
dt_tmp2$time<-dt_tmp2$time_f1
data_county_us<-rbind(dt_tmp1,dt_tmp2)

ICUbeds=fread("/Users/giancarlo/Documents/GitHub/Local-COVID19/data/KHN_ICU_bed_county_analysis_2.csv", sep=",", header=T)
names(ICUbeds)[1]="FIPS"

full_data<-merge(data_county_us,ICUbeds[,c("FIPS","all_icu","Total_pop","60plus")], on=("FIPS"),all.x = T)
full_data<-as.data.table(full_data)
full_data<-full_data[,c("time_f1","time_f2"):=NULL]
full_data<-full_data[order(full_data$time),]
#epidemic modelling in one county test




model_L2<-function(county, detection_rate, inf_prob, disease_duration,start_time){
  
  dt=full_data[FIPS==county,]
  days=length(dt$Confirmed)
  pop=dt[1,]$Total_pop
  pop60=dt[1,]$`60plus`
  
  
  params=param.dcm(
    inf.prob = inf_prob,
    #inter.eff=0.3,
    #inter.start=1
    act.rate=30,
    rec.rate=1.0/disease_duration,
    a.rate=0,
    ds.rate = 0,
    di.rate = 0,
    dr.rate =0
    )
  init=init.dcm(
    s.num=pop,
    i.num=1,
    r.num=0)
  #  s.num.g2=pop60,
  #  i.num.g2=0,
  #  r.num.g2=0
  
  control=control.dcm(
    type="SIR",
    nsteps=30
  )
    
  model=dcm(params,init, control)
  infected=model$epi$i.num
  infected_avg_detected=infected*detection_rate
  l_overlap=start_time:(start_time+days-1)
  L2<-sum((dt$Confirmed-infected_avg_detected[l_overlap,])^2)
  model<-NULL
  return(L2)
}

modelSIR<-function(county, detection_rate, inf_prob, disease_duration, start_time){
  dt=full_data[FIPS==county,]
  days=length(dt$Confirmed)
  pop=dt[1,]$Total_pop
  pop60=dt[1,]$`60plus`
  
  
  params=param.dcm(
    inf.prob = inf_prob,
    #inter.eff=0.3,
    #inter.start=1
    act.rate=30,
    rec.rate=1.0/disease_duration,
    a.rate=0,
    ds.rate = 0,
    di.rate = 0,
    dr.rate =0
  )
  init=init.dcm(
    s.num=pop,
    i.num=1,
    r.num=0)
  #  s.num.g2=pop60,
  #  i.num.g2=0,
  #  r.num.g2=0
  
  control=control.dcm(
    type="SIR",
    nsteps=30
  )
  
  model=dcm(params,init, control)
  infected=model$epi$i.num
  infected_avg_detected=infected*detection_rate
  return(list("model"=model,"min_days"=start_time))
}

fitmodel_count_imp_sampling<-function(county, detection_rate,inf_prob,disease_duration,start_time,eps1,eps2,eps3,eps4,n){
  
   for(i in 1:n){
     print(paste0("iteration: ",i))
     if(i==1){
       invL2matrix<-lapply(detection_rate,function(x){lapply(inf_prob,function(y){lapply(disease_duration,function(z){lapply(start_time,function(w){data.frame("dt"=x,"inf_p"=y,"dur"=z,"start"=w,"inv_p"=1.0/(model_L2(county,x,y,z,w)))})})})})
       invL2matrix<-do.call(rbind,unlist(unlist(unlist(invL2matrix,recursive = F),recursive = F),recursive = F))
       }
      else{
       invL2matrix<-lapply(1:250,function(x)data.frame("dt"=newvals[i,]$dt_rs,"inf_p"=newvals[i,]$infp_rs,"dur"=newvals[i,]$dur_rs,"start"=newvals[i,]$start_rs, "inv_p"=1.0/(model_L2(county,newvals[i,]$dt_rs,newvals[i,]$infp_rs,newvals[i,]$dur_rs,newvals[i,]$start_rs)^2)))
       invL2matrix<-do.call(rbind,invL2matrix)
       }
    
     invL2matrix<-invL2matrix%>%na.omit
     invL2matrix
     invL2matrix$prob<-invL2matrix$inv_p/sum(invL2matrix$inv_p)
     # for(x in detection_rate){
     #   for(y in inf_prob){
     #     for(z in disease_duration){
     #      print(c(x,y,z,1.0/model_L2(county,x,y,z)))
     #     }
     #   }
     # }
     
     resampled<-sample(1:length(invL2matrix$prob), size =250, replace=T ,  prob=invL2matrix$prob)
     newvals<-invL2matrix[resampled,]
     newvals$dt_rs<-sapply(newvals$dt,function(x)runif(1,x-eps1,x+eps1))
     newvals$infp_rs<-sapply(newvals$inf_p,function(x)runif(1,x-eps2,x+eps2))
     newvals$dur_rs<-sapply(newvals$dur,function(x)sample(c(x,x-eps3,x+eps3),1))
     newvals$start_rs<-sapply(newvals$start,function(x)sample(c(x,x-eps4,x+eps4),1))
     eps1=eps1*0.95
     eps2=eps2*0.95
     eps3=eps3
     eps4=eps4
     # detection_rate<-newvals$dt_rs
     # inf_prob<-newvals$infp_rs
     # disease_duration<-newvals$dur_rs
    }
  return(newvals)
}

#NY 1861  5039  8218 11396 14572 17746 20926 24104 27286 29860 32067 34342 36686 39098
county=6037#sample(full_data$FIPS,1) 

eps1=0.05
detection_rate=c(0.01,1:9 *eps1)
eps2=0.005
inf_prob=1:25*eps2
eps3=1
disease_duration=13+1:1*eps3
eps4=1
start_time=7+(1:10)*eps4

fitm<-fitmodel_count_imp_sampling(county,detection_rate,inf_prob,disease_duration,start_time,eps1,eps2,eps3,eps4,20)

detection_rate<-median(fitm$dt_rs)
inf_prob<-median(fitm$infp_rs)
disease_duration=median(fitm$dur_rs)
start_time=median(fitm$start_rs)
m=modelSIR(county,detection_rate,inf_prob,disease_duration,start_time)


plot(m$model$epi$i.num$run1*detection_rate,type="l")
points(m$min_days:(m$min_days+length(full_data[FIPS==county,]$Confirmed)-1),full_data[FIPS==county,]$Confirmed,pch=1)

