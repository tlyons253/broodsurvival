library(simfxns)


broodsim_ss<-function(SD1,broods,meanD,maxage){

library(purrr)
broods<-30
visits<-2
surv_day1<-0.97
p<-0.6
maxage<-22

hatch<-rbinom(broods,20,0.7)

sim.true<-matrix(NA,nrow=broods,ncol=maxage)
p.surv<-matrix(NA,nrow=broods,ncol=maxage)
sim.true[,1]<-hatch

age<-seq(1,maxage,1)

for(i in 1:broods){
	for(j in 1: (maxage-1)){
		p.surv[i,j]=plogis((qlogis(surv_day1)-0.1)+0.1*age[j])
			
		sim.true[i,j+1]<-rbinom(1,size=sim.true[i,j],p.surv[i,j])
		}
}
								

age.obs<-matrix(NA,nrow=broods,ncol=visits)


age.obs[,1]<-sample(size=broods,
				seq(15,20,1),
				replace=TRUE)
	
age.obs[,2]<-map2_int(rep(1,broods),
					ifelse(22-age.obs[,1]>3,
						3,
						22-age.obs[,1]),
					~sample(.x:.y,size=1))+age.obs[,1]


dcov<-matrix(NA,nrow=broods,ncol=maxage)

darray<-array(runif(n=broods*visits,-1,1),dim=c(broods,visits))

for(i in 1:broods){
	for (v in 1:2){
	
dcov[i,age.obs[i,v]]<-darray[i,v]
}}
		
		alpha0<-qlogis(p)
		alpha1<--1
		p.det<-plogis(alpha0+alpha1*darray)


sim.obsv<-matrix(NA,nrow=broods,ncol=maxage)	

 
		#observation data
for (i in 1:broods){
	for (v in 1:visits){
	sim.obsv[i,age.obs[i,v]]<-rbinom(n=1,
									size=sim.true[i,age.obs[i,v]],
									prob=p.det[,v])
	}
}
sim.obsv[,1]<-hatch				

test.data<-list(hatch=hatch,
				sim.obsv=sim.obsv,
				dcov=dcov,
				maxage=maxage,
				visits=visits,
				broods=broods,
				age.obs=age.obs)

return(test.data)

}




make.jagdata<-function(x){

broods=x$broods

X=x$sim.obsv
hatched = x$hatch

dcov=x$dcov
maxage=x$maxage
Age=x$age.obs
ageseq=seq(1,maxage,1)


return(list(broods=broods,
		
			X=X,
			hatched=hatched,
			dcov=dcov,
			maxage=maxage,
			Age=Age,
			ageseq=ageseq))

}







 JAG.initial<-function(x){
	library(purrr)
	C = array(map2_int(apply(x$sim.obs[,2:22], 1,max,na.rm=TRUE),
				x$hatch,
				~sample(.x:.y,size=1)),dim=c(x$broods,x$maxage))
	C[,1]<-NA
	return(
	function()


	list(s.int=runif(1,-2,2),
		d.int=runif(1,-2,2),
		B_age=runif(1,-2,2),
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


sink("C:/Users/Tim/Desktop/logexp/broodsurvival_logexp_SS.txt")
	cat("
model {
		
		for (i in 1:broods) {
		
		
			C[i,1]<-hatched[i]
			
			
			for (j in 1: (maxage-1)){
		
		#the model
		
			C[i,j+1]~dbin(surv[j],C[i,j])#C-chicks surviving to be detected
		
				#model for survival 
				}
				
				
			for (j in 1:2){

			#observation model	
				X[i,Age[i,j]]~dbin(det[i,Age[i,j]],C[i,Age[i,j]])		#count on visit j
				
				logit(det[i,Age[i,j]])<-d.int+p1*dcov[i,Age[i,j]]			#model for detection
				
				
				}
				}
		for (j in 1: (maxage-1)){
			logit(surv[j])<-s.int
				}
		
	#priors
		 s.int~dlogis(0,1)
		#B_age~dlogis(0,1)
		 p1~dlogis(0,1)
		 d.int~dlogis(0,1)
		 		 
	
logit(trueD)<-d.int
logit(dsr)<-s.int
s22<-prod(surv[1: (maxage-1)])
	
	}



",fill=TRUE)
sink()

PTM<-c('trueD','dsr','p1','s22')





library(jagsUI)

test.data<-broodsim_ss(0.97,30,0.6,22)

mydata<-make.jagdata(test.data)

inits=JAG.initial(test.data)



ni<-50000
nt=1
nb=20000



jag<-jags(data=mydata,
				inits=inits,
				 parameters.to.save=PTM,
				 model.file="C:/Users/Tim/Desktop/logexp/broodsurvival_logexp_SS.txt",
				 n.chains=3,n.thin=nt,n.iter=ni,n.burnin=nb,n.adapt=1000,
				 DIC=F,parallel=TRUE)



jag

plot(jag)