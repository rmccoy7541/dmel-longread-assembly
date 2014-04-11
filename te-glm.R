# TE analysis
library(lme4)
library(MuMIn)
library(ggplot2)

##### load data and functions #####

standardize <- function(x){
	(x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

te.results <- read.table("TE.results", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
te.data <- read.table("TE-data.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

names(te.results) <- c("name", "haveit")
te.merge <- merge(te.results, te.data, by = "name", all.x = TRUE)
te.merge$seqID <- as.numeric(gsub("%", "", te.merge$seqID))/100

# remove synonyms for family names
te.merge$family[te.merge$family=="DM88"]<-"Dm88"
te.merge$family[te.merge$family=="Doc2-element"]<-"Doc2"
te.merge$family[te.merge$family=="Doc3-element"]<-"Doc3"

# calculate copy number per family
families <- data.frame(table(te.merge$family))
names(families) <- c("family", "copies")
te.merge <- merge(te.merge, families, all.x = TRUE)
te.merge$family<-factor(te.merge$family)


##### GLM: start calculating family stats to get a p-value for families, since "copies" are a property of families #####

haveit.f <- aggregate( formula = haveit~family, data = te.merge, FUN = mean )
length.f <- aggregate( formula = length~family, data = te.merge, FUN = mean )
copies.f <- aggregate( formula = copies~family, data = te.merge, FUN = mean )
divergence.f <- aggregate( formula = divergence~family, data = te.merge, FUN = mean )
divergence.f$divergence[divergence.f$family=="Fw3"]<-NA
gc.f <- aggregate( formula = gc~family, data = te.merge, FUN = mean )

te.merge.f <- merge(haveit.f, length.f, by="family", all.x=TRUE)
te.merge.f <- merge(te.merge.f, copies.f, by="family", all.x=TRUE)
te.merge.f <- merge(te.merge.f, divergence.f, by="family", all.x=TRUE)
te.merge.f <- merge(te.merge.f, gc.f, by="family", all.x=TRUE)

#INE-1 has disproportionate leverage, test removing from analysis
te.merge.f<-te.merge.f[te.merge.f$family!="INE-1",]

te.merge.f$length.s <- standardize(te.merge.f$length)
te.merge.f$copies.s <- standardize(te.merge.f$copies)
te.merge.f$divergence.s <- standardize(te.merge.f$divergence)
te.merge.f$gc.s <- standardize(te.merge.f$gc)

f1 <- glm(haveit ~ length.s + divergence.s + copies.s + gc.s, data = te.merge.f, family=binomial, weights=te.merge.f$copies)
summary(f1)



##### GLMM: logistic response, and all predictors excluding copies #####

#te.merge<-te.merge[te.merge$family!="INE-1",]
te.merge$divergence[te.merge$family=="Fw3"]<-NA

# standarize predictors to unit variance and mean of zero
te.merge$length.s <- standardize(te.merge$length)
te.merge$divergence.s <- standardize(te.merge$divergence)
te.merge$gc.s <- standardize(te.merge$gc)


g1 <- glmer(haveit ~ length.s + divergence.s + gc.s + (1|family), data = te.merge, family=binomial, verbose=TRUE)
summary(g1)
# calculate pseudo-R^2 value designed for GLMM
r.squaredGLMM(g1) #variance explained by fixed factors, variance explained by both fixed and random factors


predict_g1 <- function(x, Coef){
  mdat <- g1@frame #model name before @frame
  newdat <- mdat[1:length(x), ]
  newdat$divergence <- mean(mdat$divergence)
  newdat$gc <- mean(mdat$gc)
  newdat$length <- mean(mdat$length)
  newdat[, Coef] <- x
  pred_vals <- predict(g1, newdat, ReForm = NA, type = 'response')
  out <- data.frame(haveit = pred_vals, x = x)
  names(out)[2] <- Coef
  return(out)
}

g1 <- glmer(haveit ~ length + divergence + gc + (1|family), data = te.merge, family=binomial, verbose=TRUE)

length.p<-predict_g1(seq(1, 60000, 100), "length")
a<-ggplot(te.merge, aes(x=length, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=length.p, aes(x=length, y=haveit), size=1.5) + xlab("Length (bp)") + ylab("Prob. accurately assembled") + xlim(0,15000)

divergence.p<-predict_g1(seq(0, 1, .01), "divergence")
b<-ggplot(te.merge, aes(x=divergence, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=divergence.p, aes(x=divergence, y=haveit), size=1.5) + xlab("Divergence") + ylab("Prob. accurately assembled") + xlim(0,1)

gc.p<-predict_g1(seq(20, 60, .5), "gc")
c<-ggplot(te.merge, aes(x=gc, y=haveit)) + geom_jitter(aes(color=factor(family)), size=1.5) + theme(legend.position = "none") + geom_line(data=gc.p, aes(x=gc, y=haveit), size=1.5) + xlab("GC content (%)") + ylab("Prob. accurately assembled")

multiplot(a,b,c, cols=3)
