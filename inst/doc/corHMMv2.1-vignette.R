## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=FALSE)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=75),tidy=TRUE)

## ---- message=FALSE, warning=FALSE---------------------------------------
set.seed(1985)
require(ape)
require(expm)
require(corHMM)
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]

## ---- fig.align="center", fig.width=6, fig.height=6----------------------
plot(phy, show.tip.label = FALSE)
data.sort <- data.frame(data[, 2], data[, 3], row.names = data[,1])
data.sort <- data.sort[phy$tip.label, ]
tiplabels(pch = 16, col = data.sort[,1]+1, cex = 0.5)
tiplabels(pch = 16, col = data.sort[,2]+3, cex = 0.5, offset = 0.5)

## ---- eval=-1------------------------------------------------------------
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
load("corHMMResults.Rsave")
MK_3state

## ------------------------------------------------------------------------
head(MK_3state$data.legend)

## ------------------------------------------------------------------------
getStateMat4Dat(data)

## ---- fig.align="center", fig.width=6, fig.height=4----------------------
plotMKmodel(MK_3state)

## ------------------------------------------------------------------------
phy = MK_3state$phy
data = MK_3state$data
model  = MK_3state$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
# run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=1, nCores=1)
# we import phytools plotSimmap for plotting
phytools::plotSimmap(simmap[[1]], fsize=.5)

## ---- eval=-1------------------------------------------------------------
HMM_3state <- corHMM(phy = phy, data = data, rate.cat = 2, model = "SYM", get.tip.states = TRUE)
HMM_3state

## ---- fig.align = "center", fig.width=10, fig.height=3-------------------
plotMKmodel(HMM_3state, display = "row")

## ------------------------------------------------------------------------
#get simmap inputs from corhmm outputs
phy = HMM_3state$phy
data = HMM_3state$data
model  = HMM_3state$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)

# run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat=2, nSim=1, nCores=1)

# we import phytools plotSimmap for plotting
phytools::plotSimmap(simmap[[1]], fsize=.5)

## ------------------------------------------------------------------------
LegendAndRateMat <- getStateMat4Dat(data)
RateMat <- LegendAndRateMat$rate.mat
RateMat

## ------------------------------------------------------------------------
pars2equal <- list(c(1,2), c(3,4))
StateMatA_constrained <- equateStateMatPars(RateMat, pars2equal)
StateMatA_constrained

## ---- eval=-1------------------------------------------------------------
MK_3state_customSYM <- corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = StateMatA_constrained)
MK_3state_customSYM

## ------------------------------------------------------------------------
RateCat1 <- getStateMat4Dat(data)$rate.mat # R1
RateCat1 <- equateStateMatPars(RateCat1, c(1:4))
RateCat1
RateCat2 <- getStateMat4Dat(data)$rate.mat # R2
RateCat2 <- dropStateMatPars(RateCat2, 3)
RateCat2

## ------------------------------------------------------------------------
RateClassMat <- getRateCatMat(2) #
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
RateClassMat

## ------------------------------------------------------------------------
StateMats <- list(RateCat1, RateCat2)
StateMats

## ------------------------------------------------------------------------
FullMat <- getFullMat(StateMats, RateClassMat)
FullMat

## ---- eval=FALSE---------------------------------------------------------
#  plotMKmodel(FullMat, rate.cat = 2, display = "row", text.scale = 0.7)

## ---- eval=-1------------------------------------------------------------
HMM_3state_custom <- corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = FullMat, node.states = "none")
HMM_3state_custom

## ---- fig.width=10, fig.height=3-----------------------------------------
plotMKmodel(HMM_3state_custom, display = "row", text.scale = 0.7)

## ------------------------------------------------------------------------
phy <- read.tree("randomBD.tree")
load("simulatedData.Rsave")
head(MFT_dat)
summary(as.factor(MFT_dat[,2])) # how many of each state do we have?

## ------------------------------------------------------------------------
MFT_LegendAndRate <- getStateMat4Dat(MFT_dat)
MFT_LegendAndRate

## ------------------------------------------------------------------------
MFT_R1 <- dropStateMatPars(MFT_LegendAndRate$rate.mat, c(2,4))
MFT_R1

## ------------------------------------------------------------------------
MFT_R2 <- dropStateMatPars(MFT_LegendAndRate$rate.mat, c(4,6))
MFT_R2

## ------------------------------------------------------------------------
MFT_R3 <- MFT_LegendAndRate$rate.mat
MFT_R3

## ------------------------------------------------------------------------
MFT_ObsStateClasses <- list(MFT_R1, MFT_R2, MFT_R3)

## ------------------------------------------------------------------------
MFT_RateClassMat <- getRateCatMat(3) # we have 3 rate classes
MFT_RateClassMat <- equateStateMatPars(MFT_RateClassMat, 1:6)

## ------------------------------------------------------------------------
MFT_FullMat <- getFullMat(MFT_ObsStateClasses, MFT_RateClassMat)
MFT_FullMat

## ---- eval=FALSE---------------------------------------------------------
#  plotMKmodel(MFT_FullMat, rate.cat = 3, display = "square", text.scale = 0.9)

## ---- eval=-1------------------------------------------------------------
MFT_res.corHMM <- corHMM(phy = phy, data = MFT_dat, rate.cat = 3, rate.mat = MFT_FullMat, node.states = "none", root.p = "maddfitz")
MFT_res.corHMM

## ------------------------------------------------------------------------
head(Precur_Dat)

## ------------------------------------------------------------------------
Precur_LegendAndMat <- getStateMat4Dat(Precur_Dat)
Precur_LegendAndMat

## ------------------------------------------------------------------------
Precur_R1 <- Precur_LegendAndMat$rate.mat
Precur_R1 <- dropStateMatPars(Precur_R1, c(1,2))
Precur_R1

## ------------------------------------------------------------------------
Precur_R2 <- Precur_LegendAndMat$rate.mat
Precur_R2 <- equateStateMatPars(Precur_R2, c(1,2))
Precur_R2

## ------------------------------------------------------------------------
RateClassMat <- getRateCatMat(2) #
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
RateClassMat

## ------------------------------------------------------------------------
Precur_FullMat <- getFullMat(list(Precur_R1, Precur_R2), RateClassMat)
Precur_FullMat[c(4,2), c(2,4)] <- 0
Precur_FullMat

## ------------------------------------------------------------------------
Precur_res.corHMM <- corHMM(phy = phy, data = Precur_Dat, rate.cat = 2, rate.mat = Precur_FullMat, root.p = "maddfitz")
Precur_res.corHMM

## ------------------------------------------------------------------------
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
Limbs <- c("Limbs", "noLimbs")[data[,2]+1]
Fings <- vector("numeric", length(phy$tip.label))
Fings[which(Limbs == "Limbs")] <- round(runif(length(which(Limbs == "Limbs")), 0, 1))
Corpo <- rep("corporeal", length(phy$tip.label))
Ont_Dat <- data.frame(sp = phy$tip.label, limbs = Limbs, fings = Fings, corp = Corpo)
head(Ont_Dat)

## ------------------------------------------------------------------------
Ont_LegendAndMat <- getStateMat4Dat(Ont_Dat)
Ont_LegendAndMat

## ---- eval=-1------------------------------------------------------------
Ont_res.corHMM <- corHMM(phy = phy, data = Ont_Dat, rate.cat = 1, rate.mat = Ont_LegendAndMat$rate.mat, node.states = "none")

## ---- eval=TRUE----------------------------------------------------------
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
getStateMat4Dat(data[,c(1,2)])

## ---- eval=TRUE----------------------------------------------------------
label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
homo_gorilla <- getMRCA(phy,tip=c("Homo_sapiens", "Gorilla_gorilla"))
homo_gorilla

## ---- eval=TRUE----------------------------------------------------------
label.vector[homo_gorilla] <- 2
phy$node.label <- label.vector[-c(1:Ntip(phy))]

## ---- fig.align="center", fig.width=6, fig.height=6----------------------
plot(phy, cex=.5)
nodelabels(phy$node.label)

## ---- echo=FALSE---------------------------------------------------------
fix.node64.estrus <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=TRUE)

## ---- eval=TRUE----------------------------------------------------------
label.vector[homo_gorilla] <- 1
phy$node.label <- label.vector[-c(1:Ntip(phy))]
fix.node64.noestrus <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=TRUE)
fix.node64.noestrus

## ---- eval=TRUE, message=FALSE-------------------------------------------
library(phangorn)
data.sort <- data.frame(data[,2], row.names=data[,1])
data.sort <- data.sort[phy$tip.label,]
dat<-as.matrix(data.sort)
rownames(dat) <- phy$tip.label
dat<-phyDat(dat,type="USER", levels=c("0","1"))
mpr.recon <- ancestral.pars(phy, dat, type = c("MPR"))
mpr.recon.converted <- ConvertPhangornReconstructions(mpr.recon)
phy$node.label <- mpr.recon.converted[(Ntip(phy)+1):length(mpr.recon.converted)]

## ---- fig.align="center", fig.width=6, fig.height=6----------------------
plot(phy, cex=.5)
nodelabels(phy$node.label)

## ---- eval=TRUE----------------------------------------------------------
fixed.parsimony.recon <- corHMM(phy, data[,c(1,2)], model="ER", rate.cat=1, fixed.nodes=TRUE)
fixed.parsimony.recon

## ---- eval=TRUE----------------------------------------------------------
label.vector <- rep(NA, Ntip(phy) + Nnode(phy))
homo_gorilla <- getMRCA(phy,tip=c("Homo_sapiens", "Gorilla_gorilla"))
label.vector[homo_gorilla] <- 1
phy$node.label <- label.vector[-c(1:Ntip(phy))]
fix.node64.noestrus <- corHMM(phy, data[,c(1,2)], model="ARD", rate.cat=2, fixed.nodes=TRUE)

## ---- eval=TRUE----------------------------------------------------------
fix.node64.noestrus$states[homo_gorilla-Ntip(phy),] 

## ---- eval=TRUE----------------------------------------------------------
sum(fix.node64.noestrus$states[homo_gorilla-Ntip(phy),])

