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

setwd(OUTdir)
list.in = system("ls *meth1*518262.Rdat",intern=TRUE)
list.in = system("ls *BDOTM*meth1*518262.Rdat",intern=TRUE)

filin="TG-BDOTM-m6d1_HW-animpsa_cal0.3_TG0.5meth1-518262.Rdat"
## Lecture des donnees de simulations
i=1
l.alpha=c()
Tobs=c()
Tsimdyn=list()
Tsimsta=list()
for(filin in list.in){
    load(filin)
    l.alpha[[i]]=alpha.TN
    Tobs=unlist(simu.dyn$l.T.mean)
    Tsimdyn[[i]]=unlist(simu.dyn$l.X.mean)
    Tsimsta[[i]]=unlist(simu.sta$l.X.mean)
    i=i+1
}

sd.obs=sd(Tobs)
m.obs=mean(Tobs)
i0=which(l.alpha==0)
sd.dyn0=sd(Tsimdyn[[i0]])
m.dyn0=mean(Tsimdyn[[i0]])

proba.RTdyn=c()
proba.RTsta=c()
for(i in 1:length(l.alpha)){
    dum=1-pnorm(mean(Tsimdyn[[i]]),mean=m.obs,sd=sd.obs)
    proba.RTdyn=c(proba.RTdyn,dum)
    dum=1-pnorm(mean(Tsimsta[[i]]),mean=m.obs,sd=sd.obs)
    proba.RTsta=c(proba.RTsta,dum)
}

l.proba=1-pnorm(c(20:24),mean=m.obs,sd=sd.obs)
l.probs=c(0.5,0.7,0.9,0.99,0.99999)
Tdum=qnorm(l.probs,mean=m.obs,sd=sd.obs)

## Trace des resultats en boxplots
filout=paste("TG_",args[1],"-m",args[4],"-meth1-518262.pdf",sep="")
pdf(filout)
par(mar=c(4,4,1,3))
plot(c(-0.2,1.1),c(14,26),type="n",xlab="alpha",ylab="TG (C)",axes=FALSE)
abline(h=Tobs[["2003"]],lty=3)
boxplot(Tobs,add=TRUE,boxwex=0.2,at=-0.1)
for(i in 1:length(l.alpha)){
    boxplot(Tsimdyn[[i]],at=l.alpha[i],axes=FALSE,add=TRUE,
            col="red",boxwex=0.2)
    boxplot(Tsimsta[[i]],at=l.alpha[i],axes=FALSE,add=TRUE,
            col="blue",boxwex=0.2)
}
axis(side=1,at=l.alpha)
axis(side=2)
axis(side=4,at=Tdum,labels=round(1/(1-l.probs)))
##axis(side=3,at=l.alpha,labels=round(1/proba.RTsta))
box()
dev.off()

"erf" <- function(x) 2 * pnorm(x * sqrt(2)) - 1
"erfc" <- function(x) 1-erf(x)
"erf.inv" <- function(x) qnorm((x + 1)/2)/sqrt(2)

"laplace.T" = function(x)
{
    l = exp(0.25*x^2) * erfc(-x/2)
    return(l)
}
