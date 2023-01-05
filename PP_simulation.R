# Load necessary libraries for simulation
library(ggplot2)
library(reshape2)
library(wesanderson)
library(RColorBrewer)
library(data.table)

## Parameters

N<-1000000 # population size
s<-5 # reinfection efficiency
D<-64 # dilution factor

## An example of simulation for 40 generations

i<-1
pop<-matrix(0,nrow=N,ncol=2)
colnames(pop)<-c("wt","mt")
pop[-1,1]<-pop[1,2]<-1
pop.dat<-as.data.table(pop)
readout<-c(i,pop.summary(pop.dat))

while(i<40){
	i<-i+1
	pop.dat<-cyc(pop.dat)
	readout<-rbind(readout,c(i,pop.summary(pop.dat)))
	if(i %% 6 == 0){
		pop.dat<-pop.dat[sample(1:nrow(pop.dat),round(nrow(pop.dat)/D,digits=0))]
	}
	print(i)
}

## The output matrix - 1) number of generations, 2) host population size, 3) wild-type genotype frequency, 4) mutation genotype frequency, 5) average copy number of the phage-plasmid 

readout 


# Necesssary functions 
## 1) Function simulating segregational drift with Binomial sampling 
segregation<-function(y,cutoff){	
	target<-c(rep(0,2*y[1]),rep(1,2*y[2]))
	
	if(y[3]<50001){
		daughter1.mt<-colSums(replicate(y[3],sample(target,y[4])))
	}
	
	if(y[3]>50000){
		pool<-combn(target,y[4])		
		ncom<-ncol(pool)
		map<-rmultinom(1,y[3],rep(1/ncom,ncom))
		nmt<-colSums(pool)
		map.dat<-as.data.table(cbind(map,nmt))
		colnames(map.dat)<-c("pick","nmt")
		map.dat<-map.dat[map.dat$pick>0]
		daughter1.mt<-map.dat[,rep(nmt,pick)]
	}
	
	daughter1.wt<-y[4]-daughter1.mt
	daughter2.wt<-(2*y[1])-daughter1.wt
	daughter2.mt<-y[4]-daughter2.wt
	offspring<-rbind(t(rbind(daughter1.wt,daughter1.mt)),t(rbind(daughter2.wt,daughter2.mt)))
	rownames(offspring)<-NULL
	offspring<-as.data.frame(offspring)
	
	return(offspring)
}

## 2) Function simulating the life cycle of phage-plasmids
cyc<-function(pop.dat){
	kill<-pop.dat[pop.dat$wt<1,.N]
	burst<- kill * s
	
	pop.dat<-pop.dat[wt>0]
	popsize<-nrow(pop.dat)
	
	if(kill>0){
		host<-sample(1:popsize,burst,replace=TRUE)
		host.unique<-unique(host)
		pop.dat.sub<-pop.dat[host.unique]
		pop.dat.others<-pop.dat[-host.unique]
		news<-as.vector(table(as.integer(factor(host))))	
		pop.dat.sub[,"mt":=(pop.dat.sub[,"mt"]+news)]	
		pop.dat<-rbindlist(list(pop.dat.sub,pop.dat.others))	
	}
	
	pop.dat1<-pop.dat[mt<1,]
	offspring1<-rbindlist(list(pop.dat1,pop.dat1))
	
	pop.dat2<-pop.dat[mt>0,]
	stat<-pop.dat2[,.N,by=.(wt,mt)]
	stat$copies<-stat$wt+stat$mt
	stat<-as.matrix(stat)	
	offspring2<-rbindlist(apply(stat,1,segregation))
	colnames(offspring2)<-c("wt","mt")
		
	out<-rbindlist(list(offspring1,offspring2))
	return(out)
}

## 3)Function for generating standard format output
pop.summary<-function(pop.dat){
	popsize<-nrow(pop.dat)
	colsum<-colSums(pop.dat)
	freq<-colsum/sum(colsum)
	load<-sum(colsum)/popsize
	output<-c(popsize,freq,load)
	names(output)<-c("host","wt","mt","copy")
	return(output)
}




