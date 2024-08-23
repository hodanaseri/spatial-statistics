library(gstat)
data(coalash)
head(coalash)
plot(coalash$x,coalash$y)
library(sp)
library(geoR)
coordinates(coalash)=c('x','y')

model=variogram(coalash~1,data=coalash)
fitmodel=fit.variogram(model,vgm('Exp'))
df=data.frame(coalash$coalash,coordinates(coalash))
df2=as.data.frame(coalash)
x=c(1,2)
y=c(16,8)
df3=data.frame(x,y)
coordinates(df3)=c('x','y')
krig=krige(coalash~1,locations = ~x+y,data=df2,newdata=df3,model=fitmodel)
krig@data$var1.pred

x=c(4,14)
y=c(2,9)
df3=data.frame(x,y)
coordinates(df3)=c('x','y')
krig=krige(coalash~1,locations = ~x+y,data=df2,newdata=df3,model=fitmodel)
krig@data$var1.pred


mean(coalash$coalash)
m=coalash$coalash<9.778558
summary(m)
p=variogram(m~1,locations = coalash)
plot(p)
pfit=fit.variogram(p,vgm('Exp'))
plot(p,model=pfit)

krig2=krige(m~1,locations = ~x+y, data=df2,newdata=df3,model=pfit)
krig2

var=coalash$x
d=c(0,45,90,135)
exp=variogram(var~1,data=coalash,alpha=d)
plot(exp)

z=cbind(coalash$x,coalash$x^2,coalash$x^3,coalash$x^4)
reg=lm(coalash~z,data=coalash)
coef(reg)
summary(reg)
plot(coalash$x,coalash$coalash-fitted(reg))
v=(variogram(coalash~coalash-fitted(reg),data=coalash))
plot(v)

library(RandomFields)
library(geoR)
df=as.geodata(coalash)


plot(variogram(coalash~1,coalash,alpha=c(0,45,90,135)))
hscat(coalash~1,coalash,breaks=c(0,5,10,15,20))
p=variogram(coalash~1,data=coalash)
plot(p)
r=variogram(coalash~1,data=coalash,estimator='cressie')

plot(p,col='blue',pch=19,main='classic vs robust')
plot(r,col='red',pch=19,add=TRUE)
par(mfrow=c(1,2))
classic=fit.variogram(p,vgm('Exp'))
robust=fit.variogram(p,vgm('Exp'),fit.method = 1)
classicplot=plot(p,classic)
lines(classic,col='red',lwd=2)

plot(p,classic)
lines(p,classic,col='red')
plot(p,robust)
plot(p,pch=19,col='black',main='classic vs robust')
plot.new()

lines(robust,col='blue',lwd=2)

lin=fit.variogram(p,vgm('Lin'))
plot(p,lin)
attributes(lin)$SSErr
sph=fit.variogram(p,vgm('Sph'))
plot(p,sph)
attributes(sph)$SSErr
exp=fit.variogram(p,vgm('Exp'))
plot(p,exp)
attributes(exp)$SSErr
gau=fit.variogram(p,vgm('Gau'))
plot(p,gau)
attributes(gau)$SSErr
wav=fit.variogram(p,vgm('Wav'))
plot(p,wav)
attributes(wav)$SSErr

hist(coalash$x)
mean=mean(coalash$x)
sd=sd(coalash$x)
median=median(coalash$x)

classic=variogram(coalash~1,data=coalash)
robust=variogram(coalash~1,data=coalash,methods='Robust',gridded=FALSE)
plot(classic)
plot(robust)

lin2=fit.variogram.gls(log(x)~1,data(coalash),vgm(1,'Sph',900,1))
plot(p,lin2)

library(robustbase)
r.lmrob=lmrob(coalash~x+y, coalash)
summary(r.lmrob)

library(georob)
plot(sample.variogram(residuals(r.lmrob), locations=coalash[, c("x","y")],
                        + lag.dist.def=1, max.lag=10, estimator="matheron"), pch=1, col="black",
       + main="sample variogram of residuals coalash~x+y")
plot(sample.variogram(residuals(r.lmrob), locations=coalash[, c("x","y")],
                        + lag.dist.def=1, estimator="qn"), pch=2, col="blue", add=TRUE)
plot(sample.variogram(residuals(r.lmrob), locations=coalash[, c("x","y")],
                        + lag.dist.def=1, estimator="ch"), pch=3, col="cyan", add=TRUE)
plot(sample.variogram(residuals(r.lmrob), locations=coalash[, c("x","y")],
                        + lag.dist.def=1, estimator="mad"), pch=4, col="orange", add=TRUE)
legend("bottomright", pch=1:4, col=c("black", "blue", "cyan", "orange"),
         + legend=paste(c("method-of-moments", "Qn", "Cressie-Hawkins", "MAD"),
                        + "estimator"), bty="n")

r.georob.m0.exp.c2=georob(coalash~x+y, coalash, locations=~x+y,
                             + variogram.model="RMexp", param=c(variance=0.1, nugget=0.9, scale=1))
summary(r.georob.m0.exp.c2)