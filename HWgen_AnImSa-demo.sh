#!/bin/sh -l
## Simulations de vagues de chaleur par Analogue Importance Sampling
## Pascal Yiou, LSCE, June 2018
## Se lance par:
## qsub -q mediump -l nodes=1:ppn=12 /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa-demo.sh
## File path needs to be adapted to local machine
module load R/3.3.2
NCPU=`wc -l < $PBS_NODEFILE`
export NCPU

JOBID=`echo ${PBS_JOBID} | cut -d. -f1`

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nStarting script at: ${start_date}\n"
echo -e ${JOBID}


varname=TG
Lsim=90
daystart=01
nsim=100
yy0=1948

## weighTN=1 -> alphaTN=0.1 Poids sur les rangs de T
## weighTN=2 -> alphaTN=1   On enleve alphaTN analogues
## weighTN=3 -> alphaTN=0.1 Poids sur les valeurs de T
weighTN=1
alphaTN=0.5

fileanalo="http://birdy.lsce.ipsl.fr:8096/outputs/4355b192-8208-11e9-8b7d-00304844f2cc/output.txt"
JOBID=518262
jobid=${JOBID}
##-----------------------------------------------------------------------

## Vagues de chaleur en ete
mostart=06
yy1=2018

staname=Orly
R CMD BATCH "--args ${staname} ${varname} ${Lsim} ${mostart} ${daystart} ${nsim} ${yy0} ${yy1} ${fileanalo} ${JOBID} ${alphaTN} ${weighTN}" /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa.R /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa${jobid}.Rout


## Vagues de chaleur en hiver
mostart=12
yy1=2016

staname=Orly
R CMD BATCH "--args ${staname} ${varname} ${Lsim} ${mostart} ${daystart} ${nsim} ${yy0} ${yy1} ${fileanalo} ${JOBID} ${alphaTN} ${weighTN}" /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa.R /home/users/yiou/RStat/A2C2/IMPSAM/HWgen_AnImSa${jobid}.Rout

##-----------------------------------------------------------------------

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nEnding script at: ${start_date}\n"
