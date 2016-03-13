#posterior distribution histograms for alpha and beta at T1 and T2 (to be inserted in manuscript)

t1=read.table('T1.chain',head=T)
t2=read.table('T2.chain',head=T)
thin=seq(1000,10000,10)
alpha.bin=seq(0,3,.25)
beta.bin=seq(0,2,.25*2/3)
alpha.yrange=c(0,1.75)
beta.yrange=c(0,2.25)

pdf('post_dist_hist.pdf',width=12,height=12)
par(mfrow=c(2,2),cex=1.4,mar=c(4,4,1,1),mgp=c(2,1,0))

hist(t1$alpha[thin],main="",prob=T,breaks=alpha.bin,xlab="",col="lightblue",ylim=alpha.yrange)
mtext(side=3,adj=.5,line=-1,expression(paste("Parameter ",alpha," at T1")),cex=1.6)

hist(t1$beta[thin],main="",prob=T,breaks=beta.bin,xlab="",col="lightyellow",ylim=beta.yrange)
mtext(side=3,adj=.5,line=-1,expression(paste("Parameter ",beta," at T1")),cex=1.6)

hist(t2$alpha[thin],main="",prob=T,breaks=alpha.bin,xlab="",col="lightblue",ylim=alpha.yrange)
mtext(side=3,adj=.5,line=-1,expression(paste("Parameter ",alpha," at T2")),cex=1.6)

hist(t2$beta[thin],main="",prob=T,breaks=beta.bin,xlab="",col="lightyellow",ylim=beta.yrange)
mtext(side=3,adj=.5,line=-1,expression(paste("Parameter ",beta," at T2")),cex=1.6)
dev.off()