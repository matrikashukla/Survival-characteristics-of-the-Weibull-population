n<-57
#SHAPE PARAMETER
alp<-1.6
#SCALE PARAMETER
be<-3.4
#LEVEL OF SIGNIFICANCE
ls<-0.05
x<-rweibull(n,alp,be)
x
hist(x,prob=T)
# LOG LIKELIHOOD FUNCTION
logL<-function(theta){
-(n*log(theta[1])-n*theta[1]*log(theta[2])
+sum((theta[1]-1)*log(x))-sum((x/theta[2])^(theta[1])))
}

#TO OPTIMIZE LOG LIKELIHOOD FUNCTION
opt<-nlm(logL,c(.1,.1))
opt
alp_ml<-opt$estimate[1];alp_ml
be_ml<-opt$estimate[2];be_ml

#TESTING OF HYPOTHESIS
#Ho:THE RANDOM SAMPLE BELONGS TO THE WEIBULL DISTRIBUTION
#KOLMOGOROV SMIRNOV TEST
tst<-ks.test(x,"pweibull",alp_ml,be_ml)$p.value;tst
if(tst<ls){
print("REJECT Ho")
}else{
print("CANNOT REJECT Ho")
}

#SIMULATION
m<-40
reject<-c()
al_sim<-c()
be_sim<-c()
for(i in 1:n){
y<-rweibull(n,alp,be)
logL<-function(theta){
-(n*log(theta[1])-n*theta[1]*log(theta[2])
+sum((theta[1]-1)*log(y))-sum((y/theta[2])^(theta[1])))
}
opt_sim<-nlm(logL,c(.1,.1))
a_sim<-opt_sim$estimate[1]
b_sim<-opt_sim$estimate[2]
al_sim[i]<-a_sim
be_sim[i]<-b_sim
tst<-ks.test(y,"pweibull",a_sim,b_sim)$p.value
if(tst<=ls){
reject[i]<-0
}else{
reject[i]<-1
}
}
prob<-length(reject[reject==1])/m;prob

#SIMULATED PARAMETERS
alp_sim<-mean(al_sim);alp_sim
beta_sim<-mean(be_sim);beta_sim

#HISTOGRAM
hist(x,prob=T,ylim=c(0,0.3))
lines(density(x),col="red")
tr<-seq(min(x),max(x),length=n)
lines(tr,dweibull(tr,alp,be),col="green")
lines(tr,dweibull(tr,alp_ml,be_ml),col="blue")
lines(tr,dweibull(tr,alp_sim,beta_sim),col=15)
legend("topright",legend=c(expression("density"),eval(substitute(expression(paste(alpha,"=",a1,", ",beta,"=",b1)),list(a1=alp,b1=be))),
eval(substitute(expression(paste(alpha[ml],"=",a2,", ",beta[ml],"=",b2)),list(a2=alp_ml,b2=be_ml))),
eval(substitute(expression(paste(alpha[sim],"=",a3,", ",beta[sim],"=",b3)),list(a3=alp_sim,b3=beta_sim)))),
col=c("red","green","blue",15),lty=1)



#PDF
p_1<-dweibull(tr,alp,be);p_1
p_2<-dweibull(tr,alp_ml,be_ml);p_2
p_3<-dweibull(tr,alp_sim,beta_sim);p_3
plot(tr,p_1,type="l",col="green",main="PDF",xlab="time",ylab="pdf")
lines(tr,p_2,col="blue")
lines(tr,p_3,col=15)
legend("topright",legend=c(eval(substitute(expression(paste(alpha,"=",a1,", ",beta,"=",b1)),list(a1=alp,b1=be))),
eval(substitute(expression(paste(alpha[ml],"=",a2,", ",beta[ml],"=",b2)),list(a2=alp_ml,b2=be_ml))),
eval(substitute(expression(paste(alpha[sim],"=",a3,", ",beta[sim],"=",b3)),list(a3=alp_sim,b3=beta_sim)))),
col=c("green","blue",15),lty=1)

#CDF
c_1<-pweibull(tr,alp,be);c_1
c_2<-pweibull(tr,alp_ml,be_ml);c_2
c_3<-pweibull(tr,alp_sim,beta_sim);c_3

plot(tr,c_1,type="l",col="green",main="CDF",xlab="time",ylab="cdf")
lines(tr,c_2,col="blue")
lines(tr,c_3,col=15)
legend("bottomright",legend=c(eval(substitute(expression(paste(alpha,"=",a1,", ",beta,"=",b1)),list(a1=alp,b1=be))),
eval(substitute(expression(paste(alpha[ml],"=",a2,", ",beta[ml],"=",b2)),list(a2=alp_ml,b2=be_ml))),
eval(substitute(expression(paste(alpha[sim],"=",a3,", ",beta[sim],"=",b3)),list(a3=alp_sim,b3=beta_sim)))),
col=c("green","blue",15),lty=1)


#RELIABILITY
r1<-(1-c_1);r1
r2<-(1-c_2);r2
r3<-(1-c_3);r3

plot(tr,r1,type="l",col="green",main="RELIABILITY",xlab="Time",ylab="reliability")
lines(tr,r2,col="blue")
lines(tr,r3,col=15)
legend("topright",legend=c(eval(substitute(expression(paste(alpha,"=",a1,", ",beta,"=",b1)),list(a1=alp,b1=be))),
eval(substitute(expression(paste(alpha[ml],"=",a2,", ",beta[ml],"=",b2)),list(a2=alp_ml,b2=be_ml))),
eval(substitute(expression(paste(alpha[sim],"=",a3,", ",beta[sim],"=",b3)),list(a3=alp_sim,b3=beta_sim)))),
col=c("green","blue",15),lty=1)


#HAZARD
h1<-p_1/r1;h1
h2<-p_2/r2;h2
h3<-p_3/r3;h3


plot(tr,h1,type="l",col="green",main="HAZARD",xlab="Time",ylab="hazard")
lines(tr,h2,col="blue")
lines(tr,h3,col=15)
legend("topleft",legend=c(eval(substitute(expression(paste(alpha,"=",a1,", ",beta,"=",b1)),list(a1=alp,b1=be))),
eval(substitute(expression(paste(alpha[ml],"=",a2,", ",beta[ml],"=",b2)),list(a2=alp_ml,b2=be_ml))),
eval(substitute(expression(paste(alpha[sim],"=",a3,", ",beta[sim],"=",b3)),list(a3=alp_sim,b3=beta_sim)))),
col=c("green","blue",15),lty=1)

#  CUMMULATIVE HAZARD
H1<-(-log(r1));H1
H2<-(-log(r2));H2
H3<-(-log(r3));H3

plot(tr,H1,type="l",col="green",main="CUMULATIVE HAZARD",xlab="time",ylab="H(t)")
lines(tr,H2,col="blue")
lines(tr,H3,col=15)
legend("topleft",legend=c(eval(substitute(expression(paste(alpha,"=",a1,", ",beta,"=",b1)),list(a1=alp,b1=be))),
eval(substitute(expression(paste(alpha[ml],"=",a2,", ",beta[ml],"=",b2)),list(a2=alp_ml,b2=be_ml))),
eval(substitute(expression(paste(alpha[sim],"=",a3,", ",beta[sim],"=",b3)),list(a3=alp_sim,b3=beta_sim)))),
col=c("green","blue",15),lty=1)