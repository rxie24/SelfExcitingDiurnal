
#CREATE FIGURE 2

#Set directory for results and to save figure
wd <- "/home/xieryan/Dissertation1/Simulations2/"

#Subplot 1
five_testlls1 <- readRDS(paste0(wd,"five_testlls1_2.rds"))
chp_testlls1 <- five_testlls1[[1]]
hp_testlls1 <- five_testlls1[[2]]
pp_testlls1 <- five_testlls1[[3]]
cp_testlls1 <- five_testlls1[[4]]
lse_testlls1 <- five_testlls1[[5]]

#Subplot 2
five_testlls2 <- readRDS(paste0(wd,"five_testlls2_2.rds"))
chp_testlls2 <- five_testlls2[[1]]
hp_testlls2 <- five_testlls2[[2]]
pp_testlls2 <- five_testlls2[[3]]
cp_testlls2 <- five_testlls2[[4]]
lse_testlls2 <- five_testlls2[[5]]

#Subplot 3
five_testlls3 <- readRDS(paste0(wd,"five_testlls3_2.rds"))
chp_testlls3 <- five_testlls3[[1]]
hp_testlls3 <- five_testlls3[[2]]
pp_testlls3 <- five_testlls3[[3]]
cp_testlls3 <- five_testlls3[[4]]
lse_testlls3 <- five_testlls3[[5]]

#Subplot 4
five_testlls4 <- readRDS(paste0(wd,"five_testlls4_2.rds"))
chp_testlls4 <- five_testlls4[[1]]
hp_testlls4 <- five_testlls4[[2]]
pp_testlls4 <- five_testlls4[[3]]
cp_testlls4 <- five_testlls4[[4]]
lse_testlls4 <- five_testlls4[[5]]

#Find limits for plots (test)
lims1 <- unlist(five_testlls1)
lims2 <- unlist(five_testlls2)
lims3 <- unlist(five_testlls3)
lims4 <- unlist(five_testlls4)

miny1 <- min(lims1,lims2)
maxy1 <- max(lims1,lims2)

miny2 <- min(lims3,lims4)
maxy2 <- max(lims3,lims4)

beta_seq <- seq(0.15,0.5,0.025)

#Plot: Trajectories of all five model comparisons
pdf(file=paste0(wd,"five_comp_final2.pdf"), width = 6, height = 6) #png: filename
par(mfrow = c(2,2), mar = c(3,3,1.5,0.5), oma = c(2, 2, 2.5, 0))

plot(beta_seq,chp_testlls1,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=50, ", phi, " = 0.2"))),
     ylim=c(miny1,maxy1))
lines(beta_seq,hp_testlls1,col="blue")
lines(beta_seq,pp_testlls1,col="purple")
lines(beta_seq,cp_testlls1,col="orange")
lines(beta_seq,lse_testlls1,col="green")
axis(2, at = seq(0, 0.2, by = 0.1), labels = seq(0, 0.2, by = 0.1))  # Custom Y-ticks

plot(beta_seq,chp_testlls2,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=500, ", phi, " = 0.2"))),
     ylim=c(miny1,maxy1))
lines(beta_seq,hp_testlls2,col="blue")
lines(beta_seq,pp_testlls2,col="purple")
lines(beta_seq,cp_testlls2,col="orange")
lines(beta_seq,lse_testlls2,col="green")

plot(beta_seq,chp_testlls3,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=50, ", phi, " = -2"))),
     ylim=c(miny2,maxy2))
lines(beta_seq,hp_testlls3,col="blue")
lines(beta_seq,pp_testlls3,col="purple")
lines(beta_seq,cp_testlls3,col="orange")
lines(beta_seq,lse_testlls3,col="green")

plot(beta_seq,chp_testlls4,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=500, ", phi, " = -2"))),
     ylim=c(miny2,maxy2))
lines(beta_seq,hp_testlls4,col="blue")
lines(beta_seq,pp_testlls4,col="purple")
lines(beta_seq,cp_testlls4,col="orange")
lines(beta_seq,lse_testlls4,col="green")

legend("bottomright",legend=c("CHP","Hawkes","Poisson","Circadian","CISE"),
       fill=c("red","blue","purple","orange","green"),cex = 0.7)


mtext("Comparison of Test Log-Likelihoods", side = 3, outer = TRUE, line = 0.5, cex = 1.4, font = 2)
mtext(bquote(beta), side = 1, outer = TRUE, line = 0.2, cex = 1)
mtext("Test Log-Likelihoods", side = 2, outer = TRUE, line = 0.2, cex = 1)
dev.off()

