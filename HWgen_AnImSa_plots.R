## Visualisations des vagues de chaleur extremes generees par importance
## sampling par analogues
## Par Pascal Yiou (LSCE), Jan. 2019
SI=Sys.info()
if(SI[[1]] == "Darwin"){
  Rsource="/Users/yiou/programmes/RStat/"
  ANAdir="/Users/yiou/data/EUCLEIA/"
  OUTdir="/Users/yiou/data/EUCLEIA/"
}
if(SI[[1]] == "Linux"){
  Rsource="/home/users/yiou/RStat/"
  ANAdir="/home/estimr1/yiou/estimr1/NCEP/"
  NCEPdir="/home/scratch01/nkadyg/A2C2/NCEP/"
  TNdir = "/home/estimr1/yiou/estimr1/WEGE/"
  Tdir = "/home/estimr1/yiou/estimr1/ECAD/"
  OUTdir="/home/estimr2/yiou/IMPSAMP/"
}
## args=(commandArgs(TRUE))
## print(args)
## if(length(args)>0){
##     filin = args[1]
## }else{
##     filin = paste(OUTdir,"TG-Orly-m6d1_extrHW-wege_cal0.3_TG0.5-518447.Rdat",
##                   sep="")
## }

## Lecture des reanalyses (SLP, Z500, Z500 detrended)
## SLP: sim_slp_NA_1948-2018.nc
## Z500: sim_hgt_NA_1948-2018.nc
## Z500 detrended: sim_hgt_dtr_NA_1948-2018.nc
library(ncdf4)
source(paste(Rsource,"readextranc.R",sep=""))
source(paste(Rsource,"imagecont.R",sep=""))

args=(commandArgs(TRUE))
print(args)
if(length(args)>0){
    filin =args[1]
    varana =args[2]
    sufana =args[3]
    filesim = args[4]
}else{
#filin = paste(NCEPdir,"sim_slp_NA_1948-2018.nc",sep="")
#filin = paste(NCEPdir,"sim_hgt_NA_1948-2018.nc",sep="")
filin = paste(NCEPdir,"sim_hgt_dtr_NA_1948-2018.nc",sep="")
# varana="slp"
    varana = "hgt"
##    sufana = " "
    sufana = "detr"
    filesim= "/home/estimr2/yiou/IMPSAMP/TG-BDOTM-m6d1_HW-animpsa_cal0.3_TG0.5meth1-518262.Rdat"
}


## filin = paste(NCEPdir,"sim_slp_NA_1948-2018.nc",sep="")
## filin = paste(NCEPdir,"sim_hgt_NA_1948-2018.nc",sep="")
## filin = paste(NCEPdir, "sim_hgt_dtr_NA_1948-2018.nc",sep="")
## varana="slp" ## "hgt"
## varana="hgt"
## sufana="detr"
## varana="hgt"
nc = nc_open(filin)
datNCEP=lirevarnc(nc,varana)
dat.NCEP.dum=sousseasmean(datNCEP$dat,datNCEP$conv.time,l.year=c(1970:1999))
datNCEP$anom=dat.NCEP.dum$anom
datNCEP$seascyc=dat.NCEP.dum$seascyc
nc_close(nc)

## Lecture des observations
varname="TG"
setwd(paste(Tdir,"ECA_",varname,sep=""))
filin = paste(varname,"_1948_BPDTM.dat",sep="")
TN = read.table(file=filin,header=TRUE)
TN[,2:ncol(TN)] = TN[,2:ncol(TN)]/10
## name.sta=list(Berlin=list(ID="STAID004005",NAM="Berlin"),
##               Paris=list(ID="STAID000038",NAM="Paris"), ## N'existe plus!
##               DeBilt=list(ID="STAID000162",NAM="De.Bilt"),
##               Madrid=list(ID="STAID000230",NAM="Madrid"),
##               Toulouse=list(ID="STAID000033",NAM="Toulouse"),
##               Orly =   list(ID="STAID011249",NAM="Orly"))
name.sta=list(Berlin=list(ID="STAID004005",NAM="Berlin"),
  Paris=list(ID="STAID000038",NAM="Paris"), ## N'existe plus!
  DeBilt=list(ID="STAID000162",NAM="De.Bilt"),
  Madrid=list(ID="STAID000230",NAM="Madrid"),
  Toulouse=list(ID="STAID000033",NAM="Toulouse"),
  Orly =   list(ID="STAID011249",NAM="Orly"),
  BDOTM = list(ID="TN.aver",NAM="BDOTM"))
## Calcul d'une moyenne des stations
TN.aver=apply(TN[,2:ncol(TN)],1,mean,na.rm=TRUE)
TN = cbind(TN,TN.aver)
## Calcul du cycle saisonnier
mmdd = TN$Date %% 10000

## Lecture des simulations
##save(file=fname,simu.dyn,alpha.cal,alpha.TN,simu.sta,args)
setwd(OUTdir)
##list.in = system("ls TG-BDOTM*m6*.Rdat",intern=TRUE)
filin=filesim
##for(filin in list.in){
    print(paste("Processing",filin))
    setwd(OUTdir)
    load(filin)

    yymin=as.integer(args[7])
    yymax=as.integer(args[8])
    staname =args[1]
    varname = args[2] ## attention, args est reecrit avec load(filin)!
    Lsim =as.integer(args[3])
    mo.start =as.integer(args[4])
    day.start =as.integer(args[5])
    nsim =as.integer(args[6])
    yymin=as.integer(args[7])
    yymax=as.integer(args[8])
    jobid=args[10]
    
    yy=yymin:yymax
## Calcul du cycle saisonnier
    TN.seascyc = tapply(TN[[name.sta[[staname]]$ID]],mmdd,mean,na.rm=TRUE)

## Premier jour de la saison consideree
    idum=which(names(TN.seascyc) == as.character(mo.start*100+1))

    TN.mean.sort=sort(unlist(simu.dyn$l.T.mean),
                      decreasing=TRUE,index.return=TRUE)
    nyear = length(simu.dyn$l.T.mean)
    ## Selection d'annees:
    ## La plus chaude, la plus froide, la mediane et 2018
    ## Si l'ete le plus chaud est 2018, on montre le 2eme plus chaud
    y1=ifelse(yy[TN.mean.sort$ix[1]]==2018,2,1)
    l.yref=c(yy[TN.mean.sort$ix[c(y1,nyear,floor(nyear/2))]],
             2018)
    
## Composites de SLP pour chaque l.ref
    SLP.lref=c()
    for(y in l.yref){
        ddum1=y*10000+601
        ddum2=y*10000+831
        dum=apply(datNCEP$anom[datNCEP$time >= ddum1 &
                               datNCEP$time <= ddum2,],2,mean)
        SLP.lref=cbind(SLP.lref,dum)                       
    }
    
## Composites de SLP pour l'ensemble des simulations dynamiques,
##pour chaque l.ref
    SLP.lref.dyn=c()
    for(y in l.yref){
        sy=as.character(y)
        SLP.dum=c()
        for(i in 1:nsim){
            d.ref=simu.dyn$l.X[[sy]][[i]][,1]
            dum = apply(datNCEP$anom[datNCEP$time %in% d.ref,],2,mean)
            SLP.dum=cbind(SLP.dum,dum)
        }
        SLP.lref.dyn=cbind(SLP.lref.dyn,apply(SLP.dum,1,mean))
    }
    
## Composites de SLP pour l'ensemble des simulations statiques,
##pour chaque l.ref
    SLP.lref.sta=c()
    for(y in l.yref){
        sy=as.character(y)
        SLP.dum=c()
        for(i in 1:nsim){
            d.ref=simu.sta$l.X[[sy]][[i]][,1]
            dum = apply(datNCEP$anom[datNCEP$time %in% d.ref,],2,mean)
            SLP.dum=cbind(SLP.dum,dum)
        }
        SLP.lref.sta=cbind(SLP.lref.sta,apply(SLP.dum,1,mean))
    }
    
# image.cont(datNCEP$lon,datNCEP$lat,apply(datNCEP$anom[datNCEP$time >= 19630601 & datNCEP$time <= 19630831,],2,mean))
ZZ=unlist(simu.sta$l.T.mean)
XX=unlist(lapply(simu.sta$l.X.mean,median))
YY=unlist(lapply(simu.dyn$l.X.mean,median))
cor.test(XX,ZZ)
cor.test(YY,ZZ)

## Time series of mean summer temperatures for all years (with boxplots) 
##filout=paste(varname,"-",staname,"-extr-stadyn_pdf.pdf",sep="")
    filout=paste(OUTdir,varname,"-",staname,"-m",mo.start,"d",day.start,
                 "_extrHW-wege_cal",
                 alpha.cal,"_meth",weight.T,"-",alpha.TN,"-pdf-",
                 jobid,".pdf",sep="")
    pdf(filout)
    ## layout(matrix(c(1:4,5,5,5,5),2,4,byrow=TRUE))
    ## for(j in 1:4){
    ##     image.cont.c(datNCEP$lon,datNCEP$lat,SLP.lref[,j],
    ##                  xlab="",ylab="",mar=c(3,3,1,1))
    ## }
    par(mar=c(4,4,1,1))
    rangeplot=round(range(c(unlist(simu.sta$l.T.mean),
                            unlist(simu.dyn$l.X.mean)),na.rm=TRUE))
    plot(c(yymin,yymax),rangeplot,type="n",xlab="Years",
         ylab=paste(varname,"(C)"))
    boxplot(simu.sta$l.X.mean,at=c(yymin:yymax)-0.2,
            add=TRUE,axes=FALSE,col="blue")
    boxplot(simu.dyn$l.X.mean,at=c(yymin:yymax)+0.2,
            add=TRUE,axes=FALSE,col="red")
    lines(c(yymin:yymax),unlist(simu.sta$l.T.mean))
    abline(v=l.yref,lty=2,col=c("red","blue","green","grey"))

    legend("bottomright",lwd=c(1,5,5),col=c("black","blue","red"),
           legend=c("Obs.","Static","Dynamic"),bty="n")
    legend("bottomleft",
           legend=paste(staname,"t start = ",mo.start,"/",day.start),
           bty="n")
    dev.off()


## Trace de trajectoires typiques (max, min, median et 2018)
    yref=2002
    filout=paste(OUTdir,varname,"-",staname,"-m",mo.start,"d",day.start,"_extrHW-wege_cal",alpha.cal,"_meth",weight.T,"-",alpha.TN,"-ts-",jobid,".pdf",sep="")
    pdf(filout,width=10)

    layout(matrix(1:4,2,2,byrow=TRUE))
    par(mar=c(4,4,1,1))
    i=1
    for(yref in l.yref){
## Ordonner les simulations par leur moyenne sur yref
        X.sort=sort(simu.dyn$l.X.mean[[as.character(yref)]],
                    index.return=TRUE)

        IJJA.yref=which(TN[,1] == (yref*100+mo.start)*100+1)
        IJJA = which(floor(TN$Date / 100) %% 100 %in% c(6,7,8))
        rangeT=range(TN[[name.sta[[staname]]$ID]][IJJA],na.rm=TRUE)
        plot(TN[[name.sta[[staname]]$ID]][IJJA.yref:(IJJA.yref+90)],
             type="l",lwd=3,
             xlab=paste("Days since ",yref,"/",mo.start,"/1",sep=""),
             ylab=paste(varname," (",staname,")",sep=""),ylim=rangeT)
        lines(TN.seascyc[idum:(idum+92)],col="blue",lwd=3)
        lines(simu.dyn$l.X[[as.character(yref)]][[X.sort$ix[100]]][,2],
              col="red",lwd=2)
        lines(simu.dyn$l.X[[as.character(yref)]][[X.sort$ix[95]]][,2],
              col="orange",lty=2)
        lines(simu.dyn$l.X[[as.character(yref)]][[X.sort$ix[50]]][,2],
              col="brown",lwd=2)
        legend("bottom",bty="n",ncol=2,
               legend=c("Seasonal cycle","Observations","Max simulation",
                        "q95 simulation","q50 simulation"),
               lwd=c(3,3,2,1,2),lty=c(1,1,1,2,1),
               col=c("blue","black","red","orange","brown"))
        legend("topleft",bty="n",letters[i])
        i=i+1
    }
    dev.off()

## Maps of SLP for the 4 identified years for: observations and optimized for
## high temperatures
    legcol=c("red","blue","green","gray")
    filout=paste(OUTdir,varname,"-",staname,"-m",mo.start,"d",day.start,
                 "_extrHW-wege_cal",
                 alpha.cal,"_meth",weight.T,"-",alpha.TN,"_ATM",
                 varana,sufana,
                 "-",jobid,".pdf",sep="")
    pdf(filout,width=12)
    layout(matrix(1:12,3,4,byrow=TRUE))
    k=1
    for(j in 1:4){
        image.cont.c(datNCEP$lon,datNCEP$lat,SLP.lref[,j],
                     xlab="",ylab="",mar=c(3,3,1,1))
        legend("bottomleft",legend=l.yref[j],bg="white",text.col=legcol[j])
        legend("top",legend="Reanalysis",bg="white")
        legend("topleft",letters[k],bg="white")
        k=k+1
    }
    
    for(j in 1:4){
        image.cont.c(datNCEP$lon,datNCEP$lat,SLP.lref.sta[,j],
                     xlab="",ylab="",mar=c(3,3,1,1))
        legend("bottomleft",legend=l.yref[j],bg="white",text.col=legcol[j])
        legend("top",legend="Sta. Simulations",bg="white")
         legend("topleft",letters[k],bg="white")
       k=k+1
    }
    
    for(j in 1:4){
        image.cont.c(datNCEP$lon,datNCEP$lat,SLP.lref.dyn[,j],
                     xlab="",ylab="",mar=c(3,3,1,1))
        legend("bottomleft",legend=l.yref[j],bg="white",text.col=legcol[j])
        legend("top",legend="Dyn. Simulations",bg="white")
         legend("topleft",letters[k],bg="white")
       k=k+1
    }
    
    dev.off()
  
##}



q("no")
