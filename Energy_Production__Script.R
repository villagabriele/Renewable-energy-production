library(forecast)
library(datasets)
library(stats)
library(ggplot2)
library(tseries)
library(lmtest)

#COMPILAMOS LOS DATOS
energy<-Renewable_Energy_Production[,2]
tail(energy)

#no hay datos no disponibles

#Tenemos datos mensuales desde el aÃ±o 1973 hasta febrero de 2023


#LO CONVERTIMOS EN SERIE TEMPORAL
en<-ts(energy,start=c(1973,1),frequency=12)
length(en) #tenemos 602 datos

autoplot(en)
#observamos una tendencia creciente, 
#la media va cambiando, por tanto no es estacionaria
#aun asi realizamos los test de estacionariedad:

adf.test(en)
#p-value= 0.8913 > 0.05 ---> Aceptamos H0, la serie no es estacionaria

kpss.test(en)
#p-value=  0.01 < 0.05 ---> Rechazamos H0--> la serie no es estacionaria

#Por tanto, sabemos que despuÃ©s tendremos que diferenciar


#Como la producciÃ³n de energÃ­a renovable depende del clima
#pensamos que porbablemenmte haya estacionalidad

#MODELIZACION DE LA PARTE ESTACIONAL
tsdisplay(en,lag.max=168,ci.type="ma")
#efectivamente se observa estacionalidad cada 12

#Ajustemos una rejilla para verlo mejor
par1<-par(mfcol=c(1,2), lab=c(14,5,7)) 

acf(en, main= "FAC Energia",168, xlim= c(0.55,14), ylim= c(0,1), ci.type="ma")
grid()

pacf(en, main= "FACP Energia",168, xlim= c(0.55,14), ylim= c(-0.4,1))
grid()

par(par1)
#la rejilla indica si se escapa o no
#ahora mismo estamos modelizando unicamente la componente estacional
#solo medimos los retardos que caen sobre la rejilla,los que estan en discontinua
#pues ahora mismo solo nos interesan los retardos multiplos de 12

#en la FAC se salen muchos, hasta el 4 estan fuera
#se observa una persistencia


#vamos a diferenciar la serie en la parte estacional:

den<- diff(en,lag=12)
autoplot(den) # vemos que la media esta centrada en el 0

tsdisplay(den,lag.max=168,ci.type="ma")
par1<-par(mfcol=c(1,2), lab=c(14,5,7)) 
acf(den, main= "FAC Energia_dif_12",168, xlim= c(0.55,14), ylim= c(-0.2,0.9), ci.type="ma")

grid()
pacf(den, main= "FACP Energia_dif_12",168, xlim= c(0.55,14), ylim= c(-0.2,0.6))
grid()
par(par1)
#ha desaparecido la persistencia
#y en fac no se sale ningun retardo

#proponemos como primer modelo:
e0<-Arima(en,order=c(0,0,0), 
          seasonal = list(order=c(0,1,0),period=12),
          include.constant = T)
print(e0)
coeftest(e0)


par1<-par(mfcol=c(1,2), lab=c(14,5,7)) 
acf(residuals(e0), main= "FAC Res E0",168, xlim= c(0.55,14))

grid()
pacf(residuals(e0), main= "FACP Res E0",168, xlim= c(0.55,14))
grid()
par(par1)

#Vemos que en FAC se sale el primer retardo y en PFAC tambien
#pensamos que hay corte en fac y decaimiento en facp
#pensamos en modelizarlo con un medias moviles, en principio de orden 1

e1<-Arima(en,order=c(0,0,0), 
          seasonal = list(order=c(0,1,1),period=12),
          include.constant = T)
print(e1)
#el coeficiente para la parte de medias moviles es -0.1268  
#si aÃ±ado y quito el doble de 0.0413 no entra el 1
coeftest(e1)


par1<-par(mfcol=c(1,2), lab=c(14,5,7)) 
acf(residuals(e1), main= "FAC Res E1",168, xlim= c(0.55,14))

grid()
pacf(residuals(e1), main= "FACP Res E1",168, xlim= c(0.55,14))
grid()
par(par1)



fac_E1<-acf(e1$residuals, 168, plot=F)
#lo escribimos en el script y vemos que hay que quitar el primero
fac_E1<-fac_E1$acf[-1]
#nos interesean los multiplos de 12, en valor abs, que escapen a las bandas de confianza
abs(fac_E1[1:14*12])>qnorm(0.975)*1/sqrt(602)
#asi sabemos bn cuale estan fuera y cuales dentro, x si no lo vemos bn en la grafica
which(abs(fac_E1[1:14*12])>qnorm(0.975)*1/sqrt(602))
#en fac res se escapan el 6,10, 11,14
#pero como no son de los primeros y aun queda modelizar la parte regular,
#podemos darlo como valido


pfac_E1<-pacf(e1$residuals, 168, plot=F)
#lo escribimos en el script y vemos que hay que quitar el primero
pfac_E1<-pfac_E1$acf[-1]
#nos interesean los multiplos de 12, en valor abs, que escapen a las bandas de confianza
abs(pfac_E1[1:14*12])>qnorm(0.975)*1/sqrt(602)
#asi sabemos bn cuale estan fuera y cuales dentro, x si no lo vemos bn en la grafica
which(abs(pfac_E1[1:14*12])>qnorm(0.975)*1/sqrt(602))
#se escapa el priemro de los retardos, pero no hemos encontrado ningun otro
#modelo en el que no se escape
#como lo que realmente nos importa es la FAC, aceptaremos este modelo como valido

checkresiduals(e1,lag.max = 168)


tsdisplay(residuals(e1),main = "Residuos E1", lag.max = 168, ci.type = "ma")


#MODELIZACIÃ“N DE LA PARTE REGULAR

#Como se sigue observando falta de estacionariedad, pues la varianza cambia
#diferenciamos ahora en la parte regular y debido al corte en FAC
#y decaimiento en PACF que observamos en las graficas anteriores
#propusimos como modelo para la parte regular (0,1,1), pero se salen los 
#primeros residuos en ACF y PACF:

e<-Arima(en,order=c(0,1,1),seasonal = list(order=c(0,1,1),period=12),include.constant=T)
print(e)
coeftest(e)

ggtsdisplay(residuals(e),main="Residuos E",lag.max=168)

#Por tanto decidimos aumentar la parte de medias moviles

e2<-Arima(en,order=c(0,1,2),seasonal = list(order=c(0,1,1),period=12),include.constant=T)

print(e2)
coeftest(e2)

ggtsdisplay(residuals(e2),main="Residuos E2",lag.max=168)

par1<-par(mfcol=c(1,2), lab=c(14,5,7)) 
acf(residuals(e2), main= "FAC Res E2",168, xlim= c(0.55,14),ylim=c(-0.15,0.15))

grid()
pacf(residuals(e2), main= "FACP Res E2",168, xlim= c(0.55,14))
grid()
par(par1)


#el primer retardo que sale en FAC y en PFAC es el 6, por lo que parece un 
#buen modelo


#Sin embargo, continuamos probando para ver si conseguimos uno mejor

e3= Arima(en,order=c(3,1,3),seasonal = list(order=c(0,1,1),period=12),include.constant=T)

print(e3)
coeftest(e3)

ggtsdisplay(residuals(e3),main="Residuos E3",lag.max=168)


par2<-par(mfcol=c(1,2), lab=c(14,5,7)) 
acf(residuals(e3), main= "FAC Res E3",168, xlim= c(0.55,14),ylim=c(-0.15,0.15))
pacf(residuals(e3), main= "FACP Res E3",168, xlim= c(0.2,14),ylim=c(-0.15,0.15))
par(par2)
#El primer retardo que sale en la FAC es el 17

#Parece que esta es una mejor modelizaciÃ³n para nuestra serie

#Comrprobemos que no hay problemas de invertibilidad
Phi<-e3$coef[-length(e3$coef)]
abs(polyroot(c(1,-Phi+qnorm(0.025)*coeftest(e3)[1:3,2]*c(1,1,1)))) 
abs(polyroot(c(1,-Phi+qnorm(0.025)*coeftest(e3)[1:3,2]*c(1,1,-1)))) 
abs(polyroot(c(1,-Phi+qnorm(0.025)*coeftest(e3)[1:3,2]*c(1,-1,1)))) 
abs(polyroot(c(1,-Phi+qnorm(0.025)*coeftest(e3)[1:3,2]*c(1,1,1)))) 

#Como vemos que hay raices menores que 1, rechazamos este modelo 
#y retomamos e2

#Estudiemos la invertibilidad para el modelo e2

Theta<- e2$coef[-length(e2$coef)]
abs(polyroot(c(1,Theta+qnorm(0.025)*coeftest(e2)[1:2,2]*c(1,1)))) 
abs(polyroot(c(1,Theta+qnorm(0.025)*coeftest(e2)[1:2,2]*c(1,-1)))) 
abs(polyroot(c(1,Theta+qnorm(0.025)*coeftest(e2)[1:2,2]*c(-1,1)))) 
abs(polyroot(c(1,Theta+qnorm(0.025)*coeftest(e2)[1:2,2]*c(-1,-1)))) 

#No hay porblemas de invertibilidad

#El modelo quedarÃ­a asi:
checkresiduals(e2,lag.max=168)


#VALIDACION:

#Ausencia de correlacion:
e2$var.coef
#se verifica, pues son todos cercanos a 0

#Normalidad:
ks.test(e2$residuals,'pnorm', mean(e2$residuals),sd(e2$residuals))
#aqui h0 es que los residuos siguen una distribucion normal, 
#tenemos pvalor>0.05--> aceptamos que esos residuos tengan una distribucion normal

#Veamos tambien los Graficos q-q:
qqnorm(e2$residuals)
qqline(e2$residuals)

#PREDICCION

#En este punto, hemos intentado prever cuál será la producción
#de energía renovable desde marzo hasta diciembre de 2023

#primero representamos el gráfico de los datos (en rojo) y del modelo (en azul)
plot(en, col="red", main = "Energia")
lines(fitted(e2),col="blue")
legend(1978,1100, legend=c("Datos","Modelo"), col= c('red','blue'), lty=c(1,1))

#predecimos la produccion de energia renovable en los proximo 5 meses:
predict(e2,10)
forecast(e2,h=10) #aquì podemos ver los intervalos de confianza tambièn
ajuste <- fitted(e2)
prde2<-forecast(ajuste,h=10) #prediccion hasta diciembre 2023 (final del año)

#grafico de las predicciones:
plot(prde2)