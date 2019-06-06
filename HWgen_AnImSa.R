## An extreme heatwave stochastic weather generator based on analogues
## and a simplified "importance sampling" algorithm
## Pascal Yiou (LSCE), June 2016, Feb 2018, June 2019
## This code is distributed freely and as is under a CeCill License.
## Se lance par:
## qsub -q mediump -l nodes=1:ppn=12 /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa.sh
## or
## qsub -q mediump -l nodes=1:ppn=12 /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa_aver.sh
## or
## qsub -q mediump -l nodes=1:ppn=12 /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa_demo.sh (provided on github)

SI=Sys.info()
if(SI[[1]] == "Linux"){
##  Rsource="/home/users/yiou/RStat/"
##  ANAdir="/home/estimr1/yiou/estimr1/NCEP/"
##  TNdir = "/home/estimr1/yiou/estimr1/WEGE/"
  Tdir = "/home/estimr1/yiou/estimr1/ECAD/" ## Needs to be adapted
  OUTdir="/home/estimr2/yiou/IMPSAMP/" ## Needs to be adapted
}

library(parallel)
##R CMD BATCH "--args ${staname} ${varname} ${Lsim} ${mostart} ${daystart} ${nsim} ${yy0} ${yy1} ${fileanalo} ${JOBID}" /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa.R /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa${jobid}.Rout

args=(commandArgs(TRUE))
print(args)
if(length(args)>0){
    staname =args[1]
    varname =args[2]
    Lsim =as.integer(args[3])
    mo.start =as.integer(args[4])
    day.start =as.integer(args[5])
    nsim =as.integer(args[6])
    yymin=as.integer(args[7])
    yymax=as.integer(args[8])
    fileanalo=args[9]
    jobid=args[10]
    alpha.TN <<- as.numeric(args[11])
    weight.T = args[12]
}else{ ## Default option
    varname="TG"
    staname="Orly"
    Lsim = 90
    mo.start =06
    day.start =01
    nsim =100
    yymin=1948
    yymax=2018
    fileanalo="http://www-lscedods.cea.fr/estimr2/DASE_NK/ana_slp_surface_base_rms_NA_latest_-80.0_50.0_22.5_70.0_1_30_20.txt"
    jobid="test"
    alpha.TN <<- 0.5
    weight.T = 1
}

## Sets days in calendar year
## Creates l.mo, l.da et moda (list of days in year)
l.mo=c(1:12)
l.da=c(31,28,31,30,31,30,31,31,30,31,30,31)
moda=c()
for(i in l.mo){
  for(j in 1:l.da[i]){
    moda=c(moda,i*100+j)
  }
}
## ------------------------------------------------------------------------
## Read input data
## Read analog file
##example: fileanalo="http://www-lscedods.cea.fr/estimr2/DASE_NK/ana_slp_surface_base_rms_NA_latest_-80.0_50.0_22.5_70.0_1_30_20.txt"
analo = read.table(fileanalo,header=TRUE)
date.a = analo$date

date.a.cal=match(as.integer(substr(analo$date,5,8)),moda)


## Read temperature  data for Berlin-De Bilt-Paris (Orly)-Toulouse-Madrid
## This is a precomputed file from ECA&D. This input file is provided in the
## distribution
##varname="TG"
setwd(paste(Tdir,"ECA_",varname,sep=""))
filin = paste(varname,"_1948_BPDTM.dat",sep="")
TN = read.table(file=filin,header=TRUE)
TN[,2:ncol(TN)] = TN[,2:ncol(TN)]/10
name.sta=list(Berlin=list(ID="STAID004005",NAM="Berlin"),
              Paris=list(ID="STAID000038",NAM="Paris"), ## N'existe plus!
              DeBilt=list(ID="STAID000162",NAM="De.Bilt"),
              Madrid=list(ID="STAID000230",NAM="Madrid"),
              Toulouse=list(ID="STAID000033",NAM="Toulouse"),
              Orly =   list(ID="STAID011249",NAM="Orly"),
              BDOTM = list(ID="TN.aver",NAM="BDOTM"))
## Daily average over stations
TN.aver=apply(TN[,2:ncol(TN)],1,mean,na.rm=TRUE)
TN = cbind(TN,TN.aver)
## ------------------------------------------------------------------------


## Simulation parameters
## Poids pour suivre le cycle saisonnier (ou etre proche du jour calendaire)
alpha.cal = 0.3

## Simulation pour une ville
##e.g. staname="Orly"
##IDsta="STAID011249"
IDsta = name.sta[[staname]]$ID
##sta=6 ## Paris (Orly)
sta = which(names(TN) == IDsta)

## ------------------------------------------------------------------------
## Computer parameters, especially for running in // on the LSCE batch cluster
## Please adapt to your own cluster
ncpus = as.numeric(Sys.getenv(c("NCPU")))

print(paste(detectCores(),"cores detected"))
nproc=3
if(SI[[1]] == "Darwin") nproc=1 ## Pour le mac portable
if(SI[[4]] %in% paste("obelix",2:6,sep="")) nproc=3 # Pour les machines interactives
if(SI[[4]] %in% paste("obelix",10:51,sep="")){
  nproc=  as.numeric(Sys.getenv(c("NCPU"))) # Pour les machines BATCH du LSCE
}
ncpus = max(nproc,ncpus,na.rm=TRUE)
print(paste("Calcul sur",ncpus,"CPUs"))


## ------------------------------------------------------------------------
## Simulation stochastique avec des poids sur la distance au jour calendaire
## a simuler, et le rang de la temperature.
## lanamax est le nombre de jours analogues autour du jour cible
## Version 1: poids sur les rangs des temperatures analogues
## Prefered method
"simu.extrHW1" = function(t.start=20030601,Lsim=90,alpha.cal = 0.5,
  alpha.TN = 0.1,sta=7,lanamax=30)
  {
    I0 = which(analo$date == t.start)
    t0.cal = match(as.integer(substr(t.start,5,8)),moda)
    t0=t.start
    T.sim=c(TN[TN$Date == t0,sta])
    t.sim=c(t0)
    ndum=c()
    for(i in 1:Lsim){
      I = which(analo$date == t0)
      I1 = I +1
      t1=analo$date[I1]
      ana.d1 = c(analo$date[I0+i],unlist(analo[I1,2:21]))
      ana.d1 = intersect(ana.d1,TN$Date)
## Poids relatif au jour calendaire pour respecter le cycle saisonnier  
      ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
      diff.ana.cal = abs(ana.d1.cal - (t0.cal+i))
      weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
      sdum=sum(weight.cal,na.rm=TRUE)
      weight.cal = weight.cal/sdum
## Poids sur la valeur de la temperature (pour avoir l'analogue le plus chaud)  
      d1=ana.d1
      Idum = match(ana.d1,TN$Date)
      TN.d1 = TN[[IDsta]][Idum] ## Temperatures analogues
      TN.d1.sort=sort(TN.d1,index.return=TRUE,decreasing=TRUE,
                      na.last=TRUE,method="radix")
      weight.TN = exp(-alpha.TN*c(1:length(TN.d1)))
## Nombre d'analogues = length(TN.d1)
      sdum=sum(weight.TN,na.rm=TRUE)
      weight.TN[TN.d1.sort$ix] = weight.TN/sdum
## Produit des poids pour normalisation des probabilites
      weight.all=weight.cal * weight.TN
      weight.all[is.na(weight.all)]=0
      weight.all = weight.all/sum(weight.all,na.rm=TRUE)
      ndum=c(ndum,length(which(weight.all >= 0.05)))
## Echantillonnage aleatoire      
      d1.max = sample(d1,size=1,prob=weight.all)
      T.sim = c(T.sim,TN[TN$Date == d1.max,sta])
      t.sim = c(t.sim, d1.max)
      t0=ifelse(idyn0,d1.max,t1)
    }
    return(cbind(t.sim,T.sim))
}
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
## Simulation stochastique avec des poids sur la distance au jour calendaire
## a simuler, et la "valeur" de la temperature.
## lanamax est le nombre de jours analogues autour du jour cible
## Version 2: poids a 0 pour alpha.TN analogues
"simu.extrHW2" = function(t.start=20030601,Lsim=90,alpha.cal = 0.5,
  alpha.TN = 1,sta=7,lanamax=30)
  {
    I0 = which(analo$date == t.start)
    t0.cal = match(as.integer(substr(t.start,5,8)),moda)
    t0=t.start
    T.sim=c(TN[TN$Date == t0,sta])
    t.sim=c(t0)
    ndum=c()
    for(i in 1:Lsim){
      I = which(analo$date == t0)
      I1 = I +1
      t1=analo$date[I1]
      ana.d1 = c(analo$date[I0+i],unlist(analo[I1,2:21]))
      ana.d1 = intersect(ana.d1,TN$Date)
## Poids relatif au jour calendaire pour respecter le cycle saisonnier  
      ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
      diff.ana.cal = abs(ana.d1.cal - (t0.cal+i))
      weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
      sdum=sum(weight.cal,na.rm=TRUE)
      weight.cal = weight.cal/sdum
## Poids sur la valeur de la temperature (pour avoir l'analogue le plus chaud)  
      d1=ana.d1
      Idum = match(ana.d1,TN$Date)
      TN.d1 = TN[[IDsta]][Idum] ## Temperatures analogues
## On enleve les alpha.TN analogues les moins chauds      
     TN.d1.sort=sort(TN.d1,index.return=TRUE,
                     na.last=TRUE,method="radix")
      weight.TN = rep(1,times=length(TN.d1))
      weight.TN[TN.d1.sort$ix[1:alpha.TN]]=0
## Nombre d'analogues = length(TN.d1)
      sdum=sum(weight.TN,na.rm=TRUE)
## Produit des poids pour normalisation des probabilites
      weight.all=weight.cal * weight.TN
      weight.all[is.na(weight.all)]=0
      weight.all = weight.all/sum(weight.all,na.rm=TRUE)
      ndum=c(ndum,length(which(weight.all >= 1/21)))
## Echantillonnage aleatoire      
      d1.max = sample(d1,size=1,prob=weight.all)
      T.sim = c(T.sim,TN[TN$Date == d1.max,sta])
      t.sim = c(t.sim, d1.max)
      t0=ifelse(idyn0,d1.max,t1)
    }
    return(cbind(t.sim,T.sim))
}

## ------------------------------------------------------------------------
## Simulation stochastique avec des poids sur la distance au jour calendaire
## a simuler, et la valeur de la temperature.
## lanamax est le nombre de jours analogues autour du jour cible
## Version 3: poids sur les ecarts au maximum des temperatures analogues
"simu.extrHW3" = function(t.start=20030601,Lsim=90,alpha.cal = 0.5,
  alpha.TN = 0.1,sta=7,lanamax=30)
  {
    I0 = which(analo$date == t.start)
    t0.cal = match(as.integer(substr(t.start,5,8)),moda)
    t0=t.start
    T.sim=c(TN[TN$Date == t0,sta])
    t.sim=c(t0)
    ndum=c()
    for(i in 1:Lsim){
      I = which(analo$date == t0)
      I1 = I +1
      t1=analo$date[I1]
      ana.d1 = c(analo$date[I0+i],unlist(analo[I1,2:21]))
      ana.d1 = intersect(ana.d1,TN$Date)
## Poids relatif au jour calendaire pour respecter le cycle saisonnier  
      ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
      diff.ana.cal = abs(ana.d1.cal - (t0.cal+i))
      weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
      sdum=sum(weight.cal,na.rm=TRUE)
      weight.cal = weight.cal/sdum
## Poids sur la valeur de la temperature (pour avoir l'analogue le plus chaud)  
      d1=ana.d1
      Idum = match(ana.d1,TN$Date)
      TN.d1 = TN[[IDsta]][Idum] ## Temperatures analogues
      TN.d1.max=max(TN.d1,na.rm=TRUE)
      weight.TN=exp(-alpha.TN*(TN.d1.max-TN.d1))
## Nombre d'analogues = length(TN.d1)
      sdum=sum(weight.TN,na.rm=TRUE)
      weight.TN=weight.TN/sdum
## Produit des poids pour normalisation des probabilites
      weight.all=weight.cal * weight.TN
      weight.all[is.na(weight.all)]=0
      weight.all = weight.all/sum(weight.all,na.rm=TRUE)
      ndum=c(ndum,length(which(weight.all >= 1/21)))
## Echantillonnage aleatoire      
      d1.max = sample(d1,size=1,prob=weight.all)
      T.sim = c(T.sim,TN[TN$Date == d1.max,sta])
      t.sim = c(t.sim, d1.max)
      t0=ifelse(idyn0,d1.max,t1)
    }
    return(cbind(t.sim,T.sim))
  }

## ------------------------------------------------------------------------
## Definition de la fonction de simulation
fun.name = paste("simu.extrHW",weight.T,sep="")
print(paste("Applying",fun.name))
SIMU.FUNC = match.fun(fun.name)

## ------------------------------------------------------------------------
## Wrapper pour les simulations stochastiques
"wrap.extrHW" = function(k)
  {
      XX = SIMU.FUNC(t.start=t.start,Lsim=Lsim,sta=sta,
                       alpha.TN=alpha.TN,alpha.cal=alpha.cal)
    return(XX)
  }
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Wrapper pour le calcul des moyennes de temperature
"wrap.mean" = function(i)
  {
    mm = mean(Xdum[[i]][,2])
    return(mm)
  }
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Simulations pour la reanalyse
"simu.yy" = function(idyn=TRUE,yymin=1948,yymax=2018,
                     mo0=mo.start,day0=day.start)
{
    idyn0 <<- idyn
    l.X.mean.dyn=list()
    l.T.mean.dyn=list()
    l.X.dyn = list()
    for(yy in yymin:yymax){
        print(paste("Processing",yy))
        t.start <<- (yy*100+mo0)*100+day0
        Xdum <<- mclapply(seq(1,nsim,by=1),wrap.extrHW,mc.cores=ncpus)
        l.X.dyn[[as.character(yy)]]=Xdum
        X.mean = mclapply(seq(1,nsim,by=1),wrap.mean,mc.cores=ncpus)
        l.X.mean.dyn[[as.character(yy)]]=unlist(X.mean)
        iTN=which(TN$Date==t.start)
        l.T.mean.dyn[[as.character(yy)]]=mean(TN[c(iTN:(iTN+Lsim)),sta],
                                              na.rm=TRUE)
    }
    return(list(l.X.mean=l.X.mean.dyn,
                l.T.mean=l.T.mean.dyn,
                l.X=l.X.dyn,ymin=yymin,ymax=yymax))
}
## ------------------------------------------------------------------------

simu.dyn=simu.yy(idyn=TRUE,yymin=yymin,yymax=yymax)
simu.sta=simu.yy(idyn=FALSE,yymin=yymin,yymax=yymax)

setwd(OUTdir)
fname=paste(varname,"-",staname,"-m",mo.start,"d",day.start,"_HW-animpsa_cal",alpha.cal,"_",varname,alpha.TN,"meth",weight.T,"-",jobid,".Rdat",sep="")

save(file=fname,simu.dyn,alpha.cal,alpha.TN,simu.sta,args,weight.T)

q("no")
## END of SCRIPT
## ------------------------------------------------------------------------


## Transformee de Laplace de la distribution uniforme entre 0 et b
"LU" = function(s,b)
{
    if(s>0){
        f = (1 - exp(-s*b))/(s*b)
    }else{
        f=1
    }
    return(f)
}

erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

"LG"=function(s)
{
    f=exp(s^2/4)*erfc(s*0.5)
    return(f)
}
