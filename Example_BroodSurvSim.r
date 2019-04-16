
library(simfxns)
library(purrr)
library(jagsUI)




brood.sim<-function(broods,visits,meanD, meanS){
broods<-broods
visits<-visits
meanS=meanS
p=meanD
hatch<-rbinom(n=broods,15,0.85)+1 						# number of eggs that hatch, preventing any 0's from occurring								#covariate value
scov<-runif(n=broods,-2,2) 								#covariate value
est.surv<-plogis(qlogis(meanS)+0.8*scov)
chicks.d1<-rbinom(n=broods,size=hatch,prob=est.surv)	# N chicks surviving to 1st detection

#detection probability

alpha0<-qlogis(p)
dcov<-array(runif(n=broods*visits,-1,1),dim=c(broods,visits))
alpha1<--1

obs.p<-plogis(alpha0+alpha1*dcov)

obsv.data<-matrix(NA,nrow=broods,ncol=visits)	
		#observation data
for (v in 1:visits){
	obsv.data[,v]<-rbinom(n=broods,size=chicks.d1,prob=obs.p[,v])
	}


#format simulated data to a matrix	
test.data<-cbind(hatch,scov,obsv.data, dcov,chicks.d1)

return(test.data)

}



make.jagdata<-function(x){

broods=nrow(x)
visits=ncol(x)-5
X=x[,3:4]
hatched = x[,1]
scov = x[,2]
pcov=x[,5:6]


return(list(broods=broods,
			visits=visits,
			X=X,
			hatched=hatched,
			scov=scov,
			pcov=pcov))

}


JAG.initial<-function(x){
library(purrr)
	C = map2_int(apply(x[,3:4], 1,max),x[,1], ~sample(.x:.y,size=1))
	return(
	function()
	list(s.int=runif(1,-2,2),
		d.int=runif(1,-2,2),
		h1=runif(1,-4,4),
		p1=runif(1,-4,4),
		C=C)
		)
	}

	
	
newfolder<-'C:/Users/tlyons4/Desktop/FinalSims'


dir.create(newfolder) #create the folder to hold models and results




sink("C:/Users/tlyons4/Desktop/FinalSims/JAGS_broodsurv_LOGIS")
	cat("
model {
		
		for (i in 1:broods) 
	   	{
		#the model
			C[i]~dbin(surv[i], hatched[i])				#C-chicks surviving to be detected
				logit(surv[i])<-s.int+h1*scov[i] 	#model for survival 
				

				for (j in 1:visits){
				
				X[i,j]~dbin(det[i,j],C[i])		#count on visit j
				logit(det[i,j])<-d.int+p1*pcov[i,j]			#model for detection
				
				
				}
				}

		
	#priors
		 s.int~dlogis(0,1)
		 h1~dlogis(0,1)
		 p1~dlogis(0,1)
		 d.int~dlogis(0,1)
		 
		 
		 # alternatives
		 # s.int~dnorm(0,0.5)
		 # h1~dnorm(0,0.5)
		 # p1~dnorm(0,0.5)
		 # d.int~dnorm(0,0.5)
	
		
logit(trueD)<-d.int
logit(trueS)<-s.int		
	
	}



",fill=TRUE)
sink()



PTM<-c('trueD','trueS','h1','p1')

 


#Above does not change between sims





#change this to some location where you want results to be saved

#----------standard data generation
#---------brood size 20 detecttion 0.2


sim.folder<-paste0(newfolder,'/brood20_D2')

dir.create(sim.folder)

out<-sim.folder





repeat.sims(SIMdata.fxn=brood.sim,
			SIMdata.params=list(broods=20,
								visits=2,
								meanD=0.2,
								meanS=0.6),
			JAGdata.fxn=make.jagdata,
			JAGS.initial=JAGS.initial,
			PTM=PTM,
			chains=3,
			thin=1,
			iter=50000,
			burn=20000,
			adapt=1000,
			outfolder.path=out,
			out_label='',
			JAGmodel.path=paste0(newfolder,'/JAGS_broodsurv_LOGIS'),
			n.sims=300)

