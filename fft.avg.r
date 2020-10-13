library(tuneR)
library(seewave)
library(ggplot2)

source('utils.r')

setWavPlayer('/usr/bin/afplay')

G3 <- 196
D4 <- 293.66
A4 <- 440
E5 <- 659.25

fmin <- 0.
fmax <- 8
amin <- -120


load.snd <- TRUE
#load.snd <- FALSE
if(load.snd)
{
# 	all <- readWave('Violins/All.wav')
# 	cao.g.4mn <- readWave('Violins/G Cao 4mn.wav')
# 	g.repeat <- readWave('Violins/G Cao 4mn.wav')
	salomon.g <- readWave('Salomon/Salomon_G.wav')
}


comp.spec <- TRUE
#comp.spec <- FALSE
if(comp.spec)
{
	spec.salomon.g <- fft.avg(salomon.g, fmax=8)
}

spec.tmp <- spec.salomon.g

f <- spec.tmp$freq
a <- spec.tmp$spec.array

idx <- which(f<=4)
f <- f[idx]
a <- a[idx,]

a.avg <- apply(a,1,median)
a.sd <- apply(a,1,sd)
ymin<-min(a.avg-a.sd)
ymax<-max(a.avg+a.sd)


plot(f, rep(0,length(f)), type='n', ylim=c(ymin, ymax), xlab='frequency kHz', ylab='dB')
lines(f, a.avg, col='blue')
lines(f,a.avg+a.sd, col='gray', lty=2)
lines(f,a.avg-a.sd, col='gray', lty=2)

ymin <- min(a)
ymax <- max(a)

quartz()
plot(f, rep(0, length(f)), type='n', ylim=c(ymin, ymax))
for(i in 1:ncol(a))
{
	lines(f, a[,i], col=i)
}



