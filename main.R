library(tidyverse);library(MASS);require(zipfR);require(fOptions)
require(Matrix);library(parallel);library(foreach);library(doMC);library(rms)
library(abind);library(flipTime);library(Triangulation);library(BPST)
library(rARPACK);library(data.table);library(Rfast);library(readr);library(refund)
library(bigmemory); library(bigalgebra)

source("source_functions.R")

nnt = 49 # number of triangles
n = N = 2000 # number of women
nboot = 50 # number of bootstraps
n1 = n2 = ncr = 40 # number of pixels in each direction
npix=n1*n2
u1=seq(0,1,length.out=n1)
v1=seq(0,1,length.out=n2)
uu=rep(u1,each=n2)
vv=rep(v1,times=n1)
Z=as.matrix(cbind(uu,vv))


if(nt == 49){
  load("V1.rda");load("Tr1.rda")
  ind.inside=BPST::inVT(Brain.V1, Brain.Tr1, Z[,1], Z[,2])$ind.inside
  V.est<-Brain.V1; Tr.est<-Brain.Tr1;
}
if(nt == 80){
  load("V2.rda");load("Tr2.rda")
  ind.inside=BPST::inVT(Brain.V2, Brain.Tr2, Z[,1], Z[,2])$ind.inside
  V.est<-Brain.V2; Tr.est<-Brain.Tr2;
}
if(nt == 144){
  load("V3.rda");load("Tr3.rda")
  ind.inside=BPST::inVT(Brain.V3, Brain.Tr3, Z[,1], Z[,2])$ind.inside
  V.est<-Brain.V3; Tr.est<-Brain.Tr3;
}

d.est = 2; r = 1
Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
B <- as.matrix(Bfull.est$B) # DIM: npts * ntriangle*((d+1)(d+2)/2)
ind.inside <- Bfull.est$Ind.inside
Q2 <- Bfull.est$Q2
K <- Bfull.est$K

surv = read_rds('sim_data')

est = est.GCV(surv, max.step = 1e3, lambda = c(0, 1, 10, 100, 1000, 1e4))

# estimated bias for NDE and IDE
NDE.bias = g.bias = (exp(est[[1]]) - NDE.true)/NDE.true
IDE.bias = ab.bias = (exp(est[[2]] %*% est[[3]]) - IDE.true)/IDE.true

# bootstrapping for standard error for the effects
IDE.all = NDE.all = NULL
for(nb in 1:nboot){
  
  samp.boot = sample(unique(surv$id), size = nrow(surv), replace = T)
  boot.df = surv[samp.boot,]
  
  est = est.GCV(data = boot.df, max.step = 1e3, lambda = c(0, 1, 10, 100, 1000, 1e4))
  
  alpha.b = as.vector(est[[2]]); beta.b = as.vector(est[[3]])
  
  NDE.b = exp(est[[1]])
  IDE.b = as.numeric(exp(alpha.b %*% beta.b))
  
  IDE.all = c(IDE.all, IDE.b) 
  NDE.all = c(NDE.all, NDE.b)
  
}




