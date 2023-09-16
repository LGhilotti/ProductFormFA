########################################
#### Analysis individual-based data ####
########################################

source("Chao_code/Subroutine_for_iNEXT.r")	#Include subroutine for iNEXT

## Step 1 (Figure 3a), sample-size-based rarefaction and extrapolation of species diversity based on the Hill numbers of order q = 0, 1, 2.
out1 <- iNEXT.Ind(Girdled)
plot.iNEXT(out1$'q=0', xlab="Number of individuals", ylab="Effective species richness", xlim=c(1, 400), ylim=c(1, 44), main="Girdled")
par(new=TRUE)
plot.iNEXT(out1$'q=1', xlab="", ylab="", xlim=c(1, 400), ylim=c(1, 44), axes=FALSE)
par(new=TRUE)
plot.iNEXT(out1$'q=2', xlab="", ylab="", xlim=c(1, 400), ylim=c(1, 44), axes=FALSE)
text(1.1 * max(out1$"q=0"[,1]), max(out1$"q=0"[,2]), "q = 0")
text(1.1 * max(out1$"q=1"[,1]), max(out1$"q=1"[,2]), "q = 1")
text(1.1 * max(out1$"q=2"[,1]), max(out1$"q=2"[,2]), "q = 2")

out2 <- iNEXT.Ind(Logged)
plot.iNEXT(out2 $'q=0', xlab="Number of individuals", ylab="Effective species richness", xlim=c(1, 600), ylim=c(1, 59), main="Logged", col=2)
par(new=TRUE)
plot.iNEXT(out2 $'q=1', xlab="", ylab="", xlim=c(1, 600), ylim=c(1, 59), col=2, axes=FALSE)
par(new=TRUE)
plot.iNEXT(out2 $'q=2', xlab="", ylab="", xlim=c(1, 600), ylim=c(1, 59), col=2, axes=FALSE)
text(1.1 * max(out2 $"q=0"[,1]), max(out2 $"q=0"[,2]), "q = 0")
text(1.1 * max(out2 $"q=1"[,1]), max(out2 $"q=1"[,2]), "q = 1")
text(1.1 * max(out2 $"q=2"[,1]), max(out2 $"q=2"[,2]), "q = 2")

## Step 1 (Figure 3b), comparison of sample-size-based rarefaction and extrapolation of two species diversity for Hill numbers of order q = 0, 1, 2.
endpoint <- 2 * min(sum(Girdled), sum(Logged))
out1 <- iNEXT.Ind(Girdled, endpoint=endpoint)
out2 <- iNEXT.Ind(Logged, endpoint=endpoint)
plot.iNEXT(out1 $'q=0', xlab="Number of individuals", ylab="Effective species richness", xlim=c(1, 400), ylim=c(1, 59), main="q=0", col=1)
par(new=TRUE)
plot.iNEXT(out2 $'q=0', xlab="", ylab="", xlim=c(1, 400), ylim=c(1, 59), main="", col=2, axes=FALSE)
legend("bottomright", c("Girdled", "Logged"), col=1:2, lty=1, pch=19, bty="n")

plot.iNEXT(out1 $'q=1', xlab="Number of individuals", ylab="Effective species richness", xlim=c(1, 400), ylim=c(1, 18), main="q=1", col=1)
par(new=TRUE)
plot.iNEXT(out2 $'q=1', xlab="", ylab="", xlim=c(1, 400), ylim=c(1, 18), main="", col=2, axes=FALSE)
legend("bottomright", c("Girdled", "Logged"), col=1:2, lty=1, pch=19, bty="n")

plot.iNEXT(out1 $'q=2', xlab="Number of individuals", ylab="Effective species richness", xlim=c(1, 400), ylim=c(1, 10), main="q=2", col=1)
par(new=TRUE)
plot.iNEXT(out2 $'q=2', xlab="", ylab="", xlim=c(1, 400), ylim=c(1, 10), main="", col=2, axes=FALSE)
legend("bottomright", c("Girdled", "Logged"), col=1:2, lty=1, pch=19, bty="n")

## Step 2 (Figure 4), sample coverage for rarefied samples and extrapolated samples.
out1 <- iNEXT.Ind(Girdled, Knots=120, rd=0)
out2 <- iNEXT.Ind(Logged, Knots=120, rd=0)
plot.iNEXT(out2$'Cov', xlab="Number of individuals", ylab="Sample coverage", xlim=c(1, 500), ylim=c(0.6, 1), col=2)
par(new=TRUE)
plot.iNEXT(out1$'Cov', xlab="", ylab="", xlim=c(1, 500), ylim=c(0.6, 1), col=1)

## Step 3 (Figure 5a), coverage-based rarefaction and extrapolation species diversity based on Hill numbers of order q = 0, 1, 2.
out1 <- iNEXT.Ind(Girdled)
plot.iNEXT(out1$'q=0'[,c(5,2,3,4)], xlab="Sample coverage", ylab="Effective species richness", xlim=c(0, 1.15), ylim=c(1, 44), main="Girdled")
par(new=TRUE)
plot.iNEXT(out1$'q=1'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1.15), ylim=c(1, 44), axes=FALSE)
par(new=TRUE)
plot.iNEXT(out1$'q=2'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1.15), ylim=c(1, 44), axes=FALSE)
text(1.1 * max(out1$"q=0"[,5]), max(out1$"q=0"[,2]), "q = 0")
text(1.1 * max(out1$"q=1"[,5]), max(out1$"q=1"[,2]), "q = 1")
text(1.1 * max(out1$"q=2"[,5]), max(out1$"q=2"[,2]), "q = 2")

out2 <- iNEXT.Ind(Logged)
plot.iNEXT(out2 $'q=0'[,c(5,2,3,4)], xlab="Sample coverage", ylab="Effective species richness", xlim=c(0, 1.15), ylim=c(1, 59), main="Logged", col=2)
par(new=TRUE)
plot.iNEXT(out2$'q=1'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1.15), ylim=c(1, 59), col=2, axes=FALSE)
par(new=TRUE)
plot.iNEXT(out2$'q=2'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1.15), ylim=c(1, 59), col=2, axes=FALSE)
text(1.1 * max(out2$"q=0"[,5]), max(out2$"q=0"[,2]), "q = 0")
text(1.1 * max(out2$"q=1"[,5]), max(out2$"q=1"[,2]), "q = 1")
text(1.1 * max(out2$"q=2"[,5]), max(out2$"q=2"[,2]), "q = 2")

## Step 3 (Figure 5b), comparison of coverage-based rarefaction and extrapolation of two species diversity for Hill numbers of order q = 0, 1, 2.
endpoint <- 2 * min(sum(Girdled), sum(Logged))
out1 <- iNEXT.Ind(Girdled, endpoint=endpoint, Knots=120, rd=0)
out2 <- iNEXT.Ind(Logged, endpoint=endpoint, Knots=120, rd=0)
plot.iNEXT(out1 $'q=0'[,c(5,2,3,4)], xlab="Sample coverage", ylab="Effective species richness", xlim=c(0, 1), ylim=c(1, 59), main="q=0", col=1)
par(new=TRUE)
plot.iNEXT(out2 $'q=0'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1), ylim=c(1, 59), main="", col=2, axes=FALSE)
legend("topleft", c("Girdled", "Logged"), col=1:2, lty=1, pch=19, bty="n")

plot.iNEXT(out1 $'q=1'[,c(5,2,3,4)], xlab="Sample coverage", ylab="Effective species richness", xlim=c(0, 1), ylim=c(1, 18), main="q=1", col=1)
par(new=TRUE)
plot.iNEXT(out2 $'q=1'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1), ylim=c(1, 18), main="", col=2, axes=FALSE)
legend("topleft", c("Girdled", "Logged"), col=1:2, lty=1, pch=19, bty="n")

plot.iNEXT(out1 $'q=2'[,c(5,2,3,4)], xlab="Sample coverage", ylab="Effective species richness", xlim=c(0, 1), ylim=c(1, 10), main="q=2", col=1)
par(new=TRUE)
plot.iNEXT(out2 $'q=2'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1), ylim=c(1, 10), main="", col=2, axes=FALSE)
legend("topleft", c("Girdled", "Logged"), col=1:2, lty=1, pch=19, bty="n")

####################################
#### Analysis sample-based data ####
#################################### 

## Since the analysis procedure of sample-based data is parallel to individual-based data, we only demo a simple example here.

#Step1, sample-size-based rarefaction and extrapolation.
out.y50 <- iNEXT.Sam(Spec=y50[-1], T=y50[1])
plot.iNEXT(out.y50$'q=0', xlab="Number of samples", ylab="Effective species richness", xlim=c(1, 1350), ylim=c(1, 300), main="Tropical ant at 50m elevation")
par(new=TRUE)
plot.iNEXT(out.y50$'q=1', xlab="", ylab="", xlim=c(1, 1350), ylim=c(1, 300), axes=FALSE)
par(new=TRUE)
plot.iNEXT(out.y50$'q=2', xlab="", ylab="", xlim=c(1, 1350), ylim=c(1, 300), axes=FALSE)
text(1.1 * max(out.y50$"q=0"[,1]), max(out.y50$"q=0"[,2]), "q = 0")
text(1.1 * max(out.y50$"q=1"[,1]), max(out.y50$"q=1"[,2]), "q = 1")
text(1.1 * max(out.y50$"q=2"[,1]), max(out.y50$"q=2"[,2]), "q = 2")

#Step2, sample coverage for rarefied samples and extrapolated samples.
plot.iNEXT(out.y50$'Cov', xlab="Number of samples", ylab="Sample coverage", ylim=c(0.8, 1), main="Tropical ant at 50m elevation")

#Step3, coverage-based rarefaction and extrapolation species diversity based on Hill numbers
plot.iNEXT(out.y50$'q=0'[,c(5,2,3,4)], xlab="Sample coverage", ylab="Effective species richness", xlim=c(0, 1.15), ylim=c(1, 300), main="Tropical ant at 50m elevation")
par(new=TRUE)
plot.iNEXT(out.y50$'q=1'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1.15), ylim=c(1, 300), axes=FALSE)
par(new=TRUE)
plot.iNEXT(out.y50$'q=2'[,c(5,2,3,4)], xlab="", ylab="", xlim=c(0, 1.15), ylim=c(1, 300), axes=FALSE)
text(1.1 * max(out.y50$"q=0"[,5]), max(out.y50$"q=0"[,2]), "q = 0")
text(1.1 * max(out.y50$"q=1"[,5]), max(out.y50$"q=1"[,2]), "q = 1")
text(1.1 * max(out.y50$"q=2"[,5]), max(out.y50$"q=2"[,2]), "q = 2")