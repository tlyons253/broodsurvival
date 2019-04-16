
#function to simulate brood survival data
#survival is a linear function of age (days since hatchting)
#restricted to 2 visits, with a max of 3 days between visits
broodsim_ss<-function(SD1,broods,meanD,maxage){

library(purrr)
broods<-broods
visits<-2
surv_day1<-SD1
p<-meanD
maxage<-maxage

hatch<-rbinom(broods,20,0.7)

sim.true<-matrix(NA,nrow=broods,ncol=maxage)
p.surv<-matrix(NA,nrow=broods,ncol=maxage)
sim.true[,1]<-hatch

age<-seq(1,maxage,1)
#survival model
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



#detection model
		alpha0<-qlogis(p)
		dcov<-array(runif(n=broods*visits,-1,1),dim=c(broods,visits))
		alpha1<--1
		p.det<-plogis(alpha0+alpha1*dcov)


sim.obsv<-matrix(NA,nrow=broods,ncol=maxage)	

 
		#observation data
for (i in 1:broods){
	for (v in 1:visits){
	sim.obsv[i,age.obs[i,v]]<-rbinom(n=1,
									size=sim.true[i,age.obs[i,v]],
									prob=p.det[,v])
	}
}
				

test.data<-list(hatch=hatch,
				sim.obsv=sim.obsv,
				dcov=dcov,
				maxage=maxage,
				visits=visits,
				broods=broods)

return(test.data)

}


test<-broodsim_dsr(SD1=0.97,broods=30,meanD=0.6,maxage=22)

#==============================================================

#same function as above, but allows for variable number of visits
#no control over time between visits
# must ensure that visits < maxage-minage
broodsim_ss2<-function(SD1,broods,meanD,minage,maxage,visits){

library(purrr)
broods<-broods
visits<-visits
surv_day1<-SD1
p<-meanD
maxage<-maxage
minage<-minage
# meanlag<-meanlag
# maxlag<-maxlag

hatch<-rbinom(broods,20,0.7)

sim.true<-matrix(NA,nrow=broods,ncol=maxage)
p.surv<-matrix(NA,nrow=broods,ncol=maxage)
sim.true[,1]<-hatch

age<-seq(1,maxage,1)


#model for survival
for(i in 1:broods){
	for(j in 1: (maxage-1)){
		p.surv[i,j]=plogis((qlogis(surv_day1)-0.1)+0.1*age[j])
			
		sim.true[i,j+1]<-rbinom(1,size=sim.true[i,j],p.surv[i,j])
		}
}
								
age.obs<-matrix(NA,nrow=broods,ncol=visits)
for(i in 1:broods){
age.obs[i,]<-sort(sample(size=visits,
				seq(minage,maxage,1),
				replace=FALSE))
	}
	


#Model for detection
		alpha0<-qlogis(p)
		dcov<-array(runif(n=broods*visits,-1,1),dim=c(broods,visits))
		alpha1<--1
		p.det<-plogis(alpha0+alpha1*dcov)


sim.obsv<-matrix(NA,nrow=broods,ncol=maxage)	

 
		#observation data
for (i in 1:broods){
	for (v in 1:visits){
	sim.obsv[i,age.obs[i,v]]<-rbinom(n=1,
									size=sim.true[i,age.obs[i,v]],
									prob=p.det[,v])
	}
}
				

test.data<-list(hatch=hatch,
				sim.obsv=sim.obsv,
				dcov=dcov,
				maxage=maxage,
				visits=visits,
				broods=broods)

return(test.data)

}


broodsim_ss2(SD1=0.97,broods=30,meanD=0.6,minage=15,maxage=22,visits=2)



#==============================================================

#same function as above, but allows for variable number of visits
#no control over time between visits
# must ensure that visits < maxage-minage
#provides additional output for use with the cheaphack model
#survival sim is a biologically plausible model that has dsr varying by age
broodsim_cheaphack<-function(SD1,broods,meanD,minage,maxage,visits){

library(purrr)
broods<-broods
visits<-visits
surv_day1<-SD1
p<-meanD
maxage<-maxage
minage<-minage
# meanlag<-meanlag
# maxlag<-maxlag

hatch<-rbinom(broods,20,0.7)

sim.true<-matrix(NA,nrow=broods,ncol=maxage)
p.surv<-matrix(NA,nrow=broods,ncol=maxage)
sim.true[,1]<-hatch

age<-seq(1,maxage,1)


#model for survival
for(i in 1:broods){
	for(j in 1: (maxage-1)){
		p.surv[i,j]=plogis((qlogis(surv_day1)-0.1)+0.1*age[j])
			
		sim.true[i,j+1]<-rbinom(1,size=sim.true[i,j],p.surv[i,j])
		}
}
								
age.obs<-matrix(NA,nrow=broods,ncol=visits)
for(i in 1:broods){
age.obs[i,]<-sort(sample(size=visits,
				seq(minage,maxage,1),
				replace=FALSE))
	}
	
expdays<-age.obs-15

#Model for detection
		alpha0<-qlogis(p)
		dcov<-array(runif(n=broods*visits,-1,1),dim=c(broods,visits))
		alpha1<--1
		p.det<-plogis(alpha0+alpha1*dcov)


sim.obsv<-matrix(NA,nrow=broods,ncol=visits)	

 
		#observation data
for (i in 1:broods){
	for (v in 1:visits){
	sim.obsv[i,v]<-rbinom(n=1,size=sim.true[i,age.obs[i,v]],
									prob=p.det[,v])
	}
}
				

test.data<-list(hatch=hatch,
				sim.obsv=sim.obsv,
				dcov=dcov,
				maxage=maxage,
				visits=visits,
				broods=broods,
				expdays=expdays)

return(test.data)

}


broodsim_cheaphack(SD1=0.97,broods=30,meanD=0.6,minage=15,maxage=22,visits=2)

#==============================================================

#same function as above, but allows for variable number of visits
#no control over time between visits
# must ensure that visits < maxage-minage
#provides additional output for use with the cheaphack model
#this version follows the model treating survival as a 2part process
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


broodsim_cheaphack2(DSR=0.995,broods=30,meanD=0.6, meanS=0.65,
					minage=15,maxage=22,visits=2)









