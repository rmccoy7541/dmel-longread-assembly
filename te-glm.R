# TE analysis
library(lme4)
library(MuMIn)
library(ggplot2)

##### load data and functions #####

standardize <- function(x, by = NULL){
	# Function to standardize traits by mean and std. deviation
	if(is.null(by)){by <- rep("NoGroups", length(x))}
	if(is.factor(by)){by <- gdata::drop.levels(by)}  # remove unused levels if by is a factor
	groups <- unique(by)  # grouping factors
	mean.sd <- cbind(tapply(x, by, mean, na.rm = TRUE), tapply(x, by, sd, na.rm = TRUE)) ## get mean and sd for each group
	out <- rep(NA, length(x))  # dummy vector that will be filled with standaradized values
	for(i in 1:length(groups)){  # For each group
		out[by == rownames(mean.sd)[i]] <- (x[by == rownames(mean.sd)[i]] - mean.sd[i,1]) / mean.sd[i,2]  # Find all the values from that group then subtract the mean and finally divide by the sd
	}
	return(out)
}

te.results <- read.table("~/Dropbox/fly assembly ms - genome research/analysis/TE-GLM/TE.results", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
te.data <- read.table("~/Dropbox/fly assembly ms - genome research/analysis/TE-GLM/TE-data.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

names(te.results) <- c("name", "haveit")
te.merge <- merge(te.results, te.data, by = "name", all.x = TRUE)
te.merge$seqID <- as.numeric(gsub("%", "", te.merge$seqID))/100

families <- data.frame(table(te.merge$family))
names(families) <- c("family", "copies")
te.merge <- merge(te.merge, families, all.x = TRUE)

te.merge$coverage[te.merge$chromosome=="2L"]<-21.3039
te.merge$coverage[te.merge$chromosome=="2R"]<-21.126
te.merge$coverage[te.merge$chromosome=="3L"]<-20.945
te.merge$coverage[te.merge$chromosome=="3R"]<-21.4801
te.merge$coverage[te.merge$chromosome=="X"]<-12.8373
te.merge$coverage[te.merge$chromosome=="4"]<-15.8183

te.merge$family<-factor(te.merge$family)

##### GLM: start calculating family stats to get a p-value for families, since "copies" are a property of families #####

haveit.f <- aggregate( formula = haveit~family, data = te.merge, FUN = mean )
length.f <- aggregate( formula = length~family, data = te.merge, FUN = mean )
copies.f <- aggregate( formula = copies~family, data = te.merge, FUN = mean )
divergence.f <- aggregate( formula = divergence~family, data = te.merge, FUN = mean )
divergence.f$divergence[divergence.f$family=="Fw3"]<-NA
coverage.f <- aggregate( formula = coverage~family, data = te.merge, FUN = mean )

te.merge.f <- merge(haveit.f, length.f, by="family", all.x=TRUE)
te.merge.f <- merge(te.merge.f, copies.f, by="family", all.x=TRUE)
te.merge.f <- merge(te.merge.f, divergence.f, by="family", all.x=TRUE)
te.merge.f <- merge(te.merge.f, coverage.f, by="family", all.x=TRUE)

#INE-1 has disproportionate leverage, test removing from analysis
te.merge.f<-te.merge.f[te.merge.f$family!="INE-1",]

te.merge.f$length.s <- standardize(te.merge.f$length)
te.merge.f$copies.s <- standardize(te.merge.f$copies)
te.merge.f$divergence.s <- standardize(te.merge.f$divergence)
te.merge.f$coverage.s <- standardize(te.merge.f$coverage)

f1 <- glm(haveit ~ length.s + divergence.s + copies.s + coverage.s, data = te.merge.f, family=binomial, weights=te.merge.f$copies)
summary(f1)

##### GLMM: logistic response, and all predictors excluding copies #####

#te.merge<-te.merge[te.merge$family!="INE-1",]
te.merge$divergence[te.merge$family=="Fw3"]<-NA

# standarize predictors to unit variance and mean of zero
te.merge$length.s <- standardize(te.merge$length)
te.merge$divergence.s <- standardize(te.merge$divergence)
te.merge$coverage.s <- standardize(te.merge$coverage)

g1 <- glmer(haveit ~ length.s + divergence.s + coverage.s + (1|family), data = te.merge, family=binomial, verbose=TRUE)
summary(g1)
# calculate pseudo-R^2 value designed for GLMM
r.squaredGLMM(g1) #variance explained by fixed factors, variance explained by both fixed and random factors

predict_g1 <- function(x, Coef){
  mdat <- g1@frame #model name before @frame
  newdat <- mdat[1:length(x), ]
  newdat$divergence <- mean(mdat$divergence)
  newdat$coverage <- mean(mdat$coverage)
  newdat$length <- mean(mdat$length)
  newdat[, Coef] <- x
  pred_vals <- predict(g1, newdat, ReForm = NA, type = 'response')
  out <- data.frame(haveit = pred_vals, x = x)
  names(out)[2] <- Coef
  return(out)
}

g1 <- glmer(haveit ~ length + divergence + coverage + (1|family), data = te.merge, family=binomial, verbose=TRUE)

length.p<-predict_g1(seq(1, 60000, 100), "length")
a<-ggplot(te.merge, aes(x=length, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=length.p, aes(x=length, y=haveit), size=1.5) + xlab("Length (bp)") + ylab("Prob. accurately assembled") + xlim(0,15000)

divergence.p<-predict_g1(seq(0, 1, .01), "divergence")
b<-ggplot(te.merge, aes(x=divergence, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=divergence.p, aes(x=divergence, y=haveit), size=1.5) + xlab("Divergence") + ylab("Prob. accurately assembled") + xlim(0,1)

coverage.p<-predict_g1(seq(12, 22, .5), "coverage")
c<-ggplot(te.merge, aes(x=coverage, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=coverage.p, aes(x=coverage, y=haveit), size=1.5) + xlab("Depth of coverage") + ylab("Prob. accurately assembled")

source('~/Dropbox/fly assembly ms - genome research/Analysis/TE-GLM/multiplot.R', chdir = TRUE)

multiplot(a,b,c, cols=3)