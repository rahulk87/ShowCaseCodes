# this script fits linear model and calculates standardized residuals
# from CRISPR screen data at different time points (i.e. T0 and T1)
# It also generates necessary plots

library(MASS)

setwd("./")

source("./functions.R")

T0 <- read.csv("T0.txt",header=TRUE,sep="\t")

T1 <- read.csv("T1.txt",header=TRUE,sep="\t")

hits <- lm(log2(T1$total.hits.1+0.5) ~   log2(T0$total.hits.0+0.5))

st.res.hits <- studres(hits)  ## Studentized residuals or standardized residuals

hits.pptm_psuedo <- lm(log2(T1$pptm_psuedo.1) ~   log2(T0$pptm_psuedo.0))

st.res.hits.pptm_psuedo <- studres(hits.pptm_psuedo) ## Studentized residuals or standardized residuals

total_hits <- hits$residuals

total_hits_pptm <- hits.pptm_psuedo$residuals

gene<-gsub("_.+","",perl=TRUE,T0$shrna.id)  ## Gene names

residuals_lm_totalhits <- cbind(
        gene,
	T0,
        T1,
        total_hits,
        total_hits_pptm,
        st.res.hits,
	st.res.hits.pptm_psuedo
	)

write.table(
       	residuals_lm_totalhits, 
        file="T0_T1.txt",
        col.names=TRUE,
        row.names=FALSE,
        sep="\t",
        quote=FALSE
        )

pdf("plots.pdf")

with(residuals_lm_totalhits,plot(log2(total.hits.0),log2(total.hits.1),pch=20,xlab="SUM149, T=0 (log2 total read counts)",ylab="SUM149 + Talazoparib (log2 total read counts)",cex.lab=1.5,cex=1.3,col="grey",cex.axis=1.5))

abline(lm(log2(T1$total.hits.1+0.5) ~ log2(T0$total.hits.0+0.5)))

with(residuals_lm_totalhits,plot(log2(pptm_psuedo.0),log2(pptm_psuedo.1),pch=20,xlab="SUM149, T=0 (log2 PPTM read counts)",ylab="SUM149 + Talazoparib (log2 PPTM read counts)",cex.lab=1.5,cex=1.3,col="grey",cex.axis=1.5))

abline(lm(log2(T1$pptm_psuedo.1) ~ log2(T0$pptm_psuedo.0)))

dev.off()

t0<-read.csv("T0_T1.txt",sep="\t",header=TRUE)
subset<-subset(t0,pptm_psuedo.1 > 1)
pdf("D48_residual_moddd.pdf")
par(oma=c(2,2,2,2))
with(subset,plot(log2(pptm_psuedo.1),st.res.hits.pptm_psuedo,pch=1,xlab="SUM149 + Talazoparib (log2 PPTM read counts)",ylab="Standardized residuals (SR)",cex.lab=1.5,cex=1.3,cex.axis=1.5,col="grey"))
with(subset(subset, st.res.hits.pptm_psuedo > 2), points(log2(pptm_psuedo.1),st.res.hits.pptm_psuedo, pch=1, col="black",cex=1.41))
with(subset(subset, gene=="PARP1" & st.res.hits.pptm_psuedo > 2), points(log2(pptm_psuedo.1),st.res.hits.pptm_psuedo, pch=16, bg=24, col="purple",cex=1.3))
with(subset(subset, gene=="TP53BP1" & st.res.hits.pptm_psuedo > 2), points(log2(pptm_psuedo.1),st.res.hits.pptm_psuedo, pch=16, bg=24, col="red",cex=1.25))
legend("topleft",legend=c("PARP1 gRNAs", "TP53BP1 gRNAs", "gRNAs with SR > 2"),cex=1,col=c("purple","red","black"),pch=c(16,16,1),pt.bg=c("black","black","black"))
abline(h=2,col = 4,lty=2)
dev.off()
