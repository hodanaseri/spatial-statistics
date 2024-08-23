library(gstat)
data(coalash)

threshold=10
coalash$indicator=ifelse(coalash$coalash>threshold,1,0)
model=variogram(indicator~1,coalash)
plot(model)
fit=fit.variogram(model,model=vgm(model='Exp'))
plot(model,fit)
new=expand.grid(x=seq(min(coalash$x),max(coalash$x)),
                      y=seq(min(coalash$y),max(coalash$y)))
coordinates(new)= ~x+y
result=krige(indicator~1,coalash,new,model=fit)
result
library(ggplot2)
df=as.data.frame(result)
ggplot(df,aes(x=x,y=y))+geom_tile(aes(fill=var1.pred))+
  scale_fill_gradient(low='blue',high='red')+
  labs(title='indicator kriging',x='x',y='y',fill='probability')

robustp=krige(indicator~1,coalash,new,model=fit,nmax=30,debug.level=-1,block=c(1,1))
robustp
robustdf=as.data.frame(robustp)
ggplot(robustdf,aes(x=x,y=y))+geom_tile(aes(fill=var1.pred))+
  scale_fill_gradient(low='blue',high='red')+
  labs(title='robust kriging',x='x',y='y',fill='probability')
