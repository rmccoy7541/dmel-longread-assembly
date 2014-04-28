# TE analysis
library(lme4)
library(MuMIn)
library(ggplot2)

##### load data and functions #####

standardize <- function(x){
	(x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

te.results <- read.table("~/Dropbox/fly assembly ms - genome research/analysis/TE-GLM/TE2/TE.results", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
te.data <- read.table("~/Dropbox/fly assembly ms - genome research/analysis/TE-GLM/TE2/TE-data2.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

names(te.results) <- c("name", "haveit")
te.merge <- merge(te.results, te.data, by = "name", all.x = TRUE)
te.merge$seqID <- as.numeric(gsub("%", "", te.merge$seqID))/100

te.merge$family[te.merge$family=="DM88"]<-"Dm88"
te.merge$family[te.merge$family=="Doc2-element"]<-"Doc2"
te.merge$family[te.merge$family=="Doc3-element"]<-"Doc3"


families <- data.frame(table(te.merge$family))
names(families) <- c("family", "copies")
te.merge <- merge(te.merge, families, all.x = TRUE)
te.merge$family<-factor(te.merge$family)


##### GLMM: logistic response, and all predictors excluding copies #####

#te.merge<-te.merge[te.merge$family!="INE-1",]
te.merge$divergence[te.merge$family=="Fw3"]<-NA

highIDcopies<-data.frame(table(te.merge[te.merge$divergence<0.01,]$family))
names(highIDcopies)<-c("family","highIDcopies")
te.merge<-merge(te.merge, highIDcopies, "family")

# standarize predictors to unit variance and mean of zero
te.merge$length.s <- standardize(te.merge$length)
te.merge$divergence.s <- standardize(te.merge$divergence)
te.merge$gc.s <- standardize(te.merge$gc)
te.merge$copies.s <- standardize(te.merge$copies)
te.merge$highIDcopies.s <- standardize(te.merge$highIDcopies)

# the full model with interaction
g1 <- glmer(haveit ~ length.s + divergence.s + highIDcopies.s*divergence.s + gc.s + (1|family), data = te.merge, family=binomial, verbose=TRUE)
summary(g1)

# the model with overall copy number (not significant)
g2 <- glmer(haveit ~ length.s + divergence.s + copies.s + gc.s + (1|family), data = te.merge, family=binomial, verbose=TRUE)
summary(g2)

# the model without the interaction term
g3 <- glmer(haveit ~ length.s + divergence.s + highIDcopies.s + gc.s + (1|family), data = te.merge, family=binomial, verbose=TRUE)
summary(g3)

# calculate pseudo-R^2 value designed for GLMM
r.squaredGLMM(g1) #variance explained by fixed factors, variance explained by both fixed and random factors


predict_g1_mean <- function(x, Coef){
  mdat <- g1@frame #model name before @frame
  newdat <- mdat[1:length(x), ]
  newdat$divergence <- mean(mdat$divergence)
  newdat$gc <- mean(mdat$gc)
  newdat$length <- mean(mdat$length)
  newdat$highIDcopies <-mean(mdat$highIDcopies)
  newdat[, Coef] <- x
  pred_vals <- predict(g1, newdat, type = 'response')
  out <- data.frame(haveit = pred_vals, x = x)
  names(out)[2] <- Coef
  return(out)
}

predict_g1_max <- function(x, Coef){
  mdat <- g1@frame #model name before @frame
  newdat <- mdat[1:length(x), ]
  newdat$divergence <- mean(mdat$divergence)
  newdat$gc <- mean(mdat$gc)
  newdat$length <- mean(mdat$length)
  newdat$highIDcopies <-max(mdat$highIDcopies)
  newdat[, Coef] <- x
  pred_vals <- predict(g1, newdat, type = 'response')
  out <- data.frame(haveit = pred_vals, x = x)
  names(out)[2] <- Coef
  return(out)
}

predict_g1_min <- function(x, Coef){
  mdat <- g1@frame #model name before @frame
  newdat <- mdat[1:length(x), ]
  newdat$divergence <- mean(mdat$divergence)
  newdat$gc <- mean(mdat$gc)
  newdat$length <- mean(mdat$length)
  newdat$highIDcopies <-min(mdat$highIDcopies)
  newdat[, Coef] <- x
  pred_vals <- predict(g1, newdat, type = 'response')
  out <- data.frame(haveit = pred_vals, x = x)
  names(out)[2] <- Coef
  return(out)
}

g1 <- glmer(haveit ~ length + divergence + gc + highIDcopies + (1|family), data = te.merge, family=binomial, verbose=TRUE)

length.p<-predict_g1_mean(seq(1, 60000, 100), "length")
a<-ggplot(te.merge, aes(x=length, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=length.p, aes(x=length, y=haveit), size=1.5) + xlab("Length (bp)") + ylab("Prob. accurately assembled") + xlim(0,15000)

divergence.p.max<-predict_g1_max(seq(0, 1, .005), "divergence")
divergence.p.min<-predict_g1_min(seq(0, 1, .005), "divergence")
b<-ggplot(te.merge, aes(x=divergence, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=divergence.p.max, aes(x=divergence, y=haveit), size=1.5) + geom_line(data=divergence.p.min, aes(x=divergence, y=haveit), linetype="21", size=1.5) + xlab("Divergence") + ylab("Prob. accurately assembled") + xlim(0,1)

gc.p<-predict_g1_mean(seq(20, 60, .5), "gc")
c<-ggplot(te.merge, aes(x=gc, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=gc.p, aes(x=gc, y=haveit), size=1.5) + xlab("GC content (%)") + ylab("Prob. accurately assembled")

highIDcopies.p<-predict_g1_mean(seq(0, 150, 1), "highIDcopies")
d<-ggplot(te.merge, aes(x=highIDcopies, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=highIDcopies.p, aes(x=highIDcopies, y=haveit), size=1.5) + xlab("No. high identity copies") + ylab("Prob. accurately assembled")

source('multiplot.R')
multiplot(a,b,c,d, cols=2)
