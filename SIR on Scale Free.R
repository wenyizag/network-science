

rm(list=ls())
library("igraph")
#――――――――――――――――――――――――――――――――――――――――――――――
#initialize the matrix with 3 nodes
a=matrix(0,nrow=3,ncol=3)
a[2,3]=1
a[1,2]=1
a[1,3]=1
a
xulie=c(1,2, 1,3, 2,3)

for(i in 3 :100){
  a=rbind(cbind(a,rep(0,i)),rep(0,i+1))
  
  #preferential attachment
  c=apply(a,2,sum)
  d=sum(apply(a,2,sum))
  e=c/d
  x=sample(1:(length(e)-1),1,replace = T ,prob =e[1:(length(e)-1)])
  a[x,i+1]=1
  
  #Symmetric matrices
  for(k in 1:i) {
    for(m in (k+1):(i+1)){
      a[m,k]=a[k,m]
    }
  }
  xulie=c(xulie,i+1,which(a[,i+1]!=0))
  g <- graph(xulie, directed=F)
}
par(mfrow=c(2,1))
#plot(u, layout=layout.fruchterman.reingold)
plot(g, layout=layout.fruchterman.reingold)
#_______________________________________________________


closeAllConnections()
graphics.off()
rm(list=ls()) #clear all variables

library(igraph) # Load the igraph package

#### MODEL PARAMETERS ####
N = 1000; #number of nodes
av.dg = 8; #average degree
m = av.dg/2; # parameter of the BA model
q = 0.1 # rewiring probability in the WS model
p = av.dg/N # probability in the ER model

## MODELS ###
# BA network
G <- barabasi.game(N, m = av.dg/2, power=1.5, directed = FALSE)
#str = "BA"


#### SIR MODEL ####
# states: S:0 I:1 R:2
## Parameters of the SIR model
mu = 1 # probability of recovering
beta = 0.4 # probability of transmission

targetnodes = seq(1,10) # node that will be used as the seed nodes
Tmax = 20
Ninf = matrix(0,nrow = length(targetnodes), ncol = Tmax) # matrix that stores the number of infected nodes at each time step
# Ninf[i,t] yields the number of infected nodes at time t when the infection starts on i
Nrec = matrix(0,nrow = length(targetnodes), ncol = Tmax)
Nsus = matrix(0,nrow = length(targetnodes), ncol = Tmax)
for(i in targetnodes){
  vstates = matrix(0, nrow = N, ncol = 1)
  vstates[i] = 1
  vinfected = which(vstates %in% 1)
  t = 1
  # while there are infected nodes
  while(length(vinfected) > 0){ 
    vinfected = which(vstates %in% 1)
    # try to recover all infected nodes at step (t+1)
    for(j in vinfected){
      if(runif(1,0,1) <= mu){
        vstates[j] = 2
      }
    }
    # try to infect all the neighbors of infected nodes at step t
    for(j in vinfected){
      ng = neighbors(G, j)
      for(k in ng){
        if(runif(1,0,1) <= beta){
          if(vstates[k] == 0){# infect only susceptible nodes
            vstates[k] = 1
          }
        }
      }
    }
    Nsus[i, t] = length(which(vstates %in% 0))/N
    Nrec[i, t] = length(which(vstates %in% 2))/N
    Ninf[i, t] = length(which(vstates %in% 1))/N # store the fraction of infected nodes at time t
    #print(paste('t:', t, 'rhoi', length(which(vstates %in% 1))/N))
    t = t + 1
  }
}
rhoi = colMeans(Ninf) # average number if infected nodes from the result of each seed node
t = seq(1,length(rhoi)) # time steps

plot(t, rhoi, xlab = "Day", ylab = "Fraction of infected nodes BA alfa=1.5",
     col = 'gray', lwd=2,ylim = c(0,0.5), xlim = c(0,Tmax), pch = 21,  bg = "red", type="o")

avgRec = colMeans(Nrec)
t = seq(1,length(avgRec)) # time steps
plot(t, avgRec, xlab = "Day", ylab = "Fraction of recovered nodes BA alfa=1.5",
     col = 'gray', lwd=2,ylim = c(0,1), xlim = c(0,Tmax), pch = 21,  bg = "green", type="o")

avgSus = colMeans(Nsus)
t = seq(1,length(avgSus)) # time steps
plot(t, avgSus, xlab = "Day", ylab = "Fraction of susceptible nodes BA alfa=1.5",
     col = 'gray', lwd=2,ylim = c(0,1), xlim = c(0,Tmax), pch = 21,  bg = "blue", type="o")

