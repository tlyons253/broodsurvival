library(simfxns)

broodsim_cheaphack2<-function(DSR,broods,meanD,meanS, minage,maxage,visits){

library(purrr)
broods<-broods
visits<-visits
DSR<-DSR
p<-meanD
maxage<-maxage
minage<-minage
surv15d<-meanS


hatch<-rbinom(broods,20,0.7)

sim.true<-matrix(NA,nrow=broods,ncol= (visits+1))
sim.true[,1]<-hatch

#scov<-runif(n=broods,-2,2) 								#covariate value
#surv15d<-plogis(qlogis(meanS)+0.8*scov)



age.obs<-matrix(NA,nrow=broods,ncol= (visits))
for(i in 1:broods){
age.obs[i,]<-sort(sample(size=visits,
				seq(minage,maxage,1),
				replace=FALSE))
	}
	
expdays<-age.obs-15
p.surv<-matrix(NA,nrow=broods,ncol=visits)

for(i in 1:broods){
	for(v in 1:visits){
	
p.surv[i,v]<-surv15d*(DSR^expdays[i,v])
#add in indexing for surv15d if it gets brood-level covariates

sim.true[i,v+1]<-rbinom(1,prob=p.surv[i,v],size=sim.true[i,v])
	}
} 



#Model for detection
		alpha0<-qlogis(p)
		dcov<-array(runif(n=broods*visits,-1,1),dim=c(broods,visits))
		alpha1<--1
		p.det<-plogis(alpha0+alpha1*dcov)


sim.obsv<-matrix(NA,nrow=broods,ncol=visits+1)	
sim.obsv[,1]<-hatch
 
		#observation data
for (i in 1:broods){
	for (v in 1:visits){
	sim.obsv[i,v+1]<-rbinom(1,size=sim.true[i,v+1],prob=p.det[i,v])
	}
}
				

test.data<-list(hatch=hatch,
				sim.obsv=sim.obsv,
				dcov=dcov,
				maxage=maxage,
				visits=visits,
				broods=broods,
				expdays=expdays,
				sim.true=sim.true)

return(test.data)

}






make.jagdata<-function(x){

return(list(
broods=x$broods,

X=x$sim.obsv,
hatched = x$hatch,
expdays= x$expdays,
dcov=x$dcov,
visits=x$visits


))

}



C = array(map2_int(apply(test.data$sim.obs[,2:3], 1,max,na.rm=TRUE),
				x$hatch,
				~sample(.x:.y,size=1)),dim=c(x$broods,x$maxage))
	C[,1]<-NA



 JAG.initial<-function(x){
	library(purrr)
	
	C = array(map2_int(apply(x$sim.obsv[,2:3], 1,max,na.rm=TRUE),
				x$hatch,
				~sample(.x:.y,size=1)),dim=c(x$broods,x$visits + 1))
	C[,1]<-NA
	return(
	function()


	list(s.int=runif(1,-2,2),
		d.int=runif(1,-2,2),
		dsr_lin=runif(1,-2,2),
		p1=runif(1,-2,2),
		C=C
				)
				)
	
	}

	


#Above does not change between sims





#change this to some location where you want results to be saved
newfolder<-'C:/Users/Tim/Desktop/logexp'


dir.create(newfolder) #create the folder to hold models and results

#redo to state-space


sink("C:/Users/Tim/Desktop/logexp/broodsurvival_logexp_cheaphack.txt")
	cat("
model {
		
		
		
		for (i in 1:broods) {
		
		
			C[i,1]<-hatched[i]
			
			
			for (j in 1: visits) {
		
		#the model
		
			C[i,j+1]~dbin(surv[i,j],C[i,j])#C-chicks surviving to be detected
		
				#model for survival 
				
		
				logit(S15d[i,j])<-s.int	#model for survival to 15d (age at first count)
				logit(DSR[i,j])<-dsr_lin	#model for daily survival after 15d
				
				surv[i,j]<-S15d[i,j]*pow(DSR[i,j],expdays[i,j])
	
			#observation model	
				X[i,j]~dbin(det[i,j],C[i,j])		#count on visit j
				logit(det[i,j])<-d.int+p1*dcov[i,j]			#model for detection
				
				
				}
				}

		
	#priors
		 s.int~dlogis(0,1)
		dsr_lin~dlogis(0,1)
		 p1~dlogis(0,1)
		 d.int~dlogis(0,1)
		 		 
	
logit(trueD)<-d.int
logit(trueS)<-s.int
logit(dsr)<-dsr_lin
	
	}



",fill=TRUE)
sink()

PTM<-c('trueD','dsr','trueS')

sim.folder<-paste0(newfolder,'/results')

dir.create(sim.folder)

out<-sim.folder



library(jagsUI)


test.data<-broodsim_cheaphack2(DSR=0.995,broods=100,meanD=0.8, meanS=0.65,
					minage=15,maxage=22,visits=5)

mydata<-make.jagdata(test.data)

inits=JAG.initial(test.data)



ni<-50000
nt=1
nb=20000



jag<-jags(data=mydata,
				inits=inits,
				 parameters.to.save=PTM,
				 model.file="C:/Users/Tim/Desktop/logexp/broodsurvival_logexp_cheaphack.txt",
				 n.chains=3,n.thin=nt,n.iter=ni,n.burnin=nb,n.adapt=1000,
				 DIC=F,parallel=TRUE)

jag

plot(jag)


