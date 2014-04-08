library(ggplot2)

mismatches<-read.table("mismatch.R.txt")
names(mismatches)<-c("id", "read_index", "read_base", "ref_base", "read_length")
a<-hist(mismatches$read_index, breaks=seq(0, max(mismatches$read_length, 1)))

lengthhist<-c()
for (i in 1:max(mismatches$read_length)) {
	print(i)
	counts<-sum(mismatches$read_length>=i)
	lengthhist<-c(lengthhist, counts)
}

error.frame<-as.data.frame(cbind(seq(1, max(mismatches$read_length), 1), a$counts/lengthhist))
names(error.frame)<-c("positions", "mm_counts")


error_profile <- ggplot(data=error.frame, aes(x=positions, y=mm_counts, colour="mismatches")) + stat_smooth()




b<-hist(mismatches$read_index/mismatches$read_length, breaks="FD")
#each of 100 bins in the histogram contains (the number of total aligned bases/100) bases
relative_error.frame<-as.data.frame(cbind(seq(0.01,1,0.01), b$counts/(4127039238/100)))
names(relative_error.frame)<-c("relative_positions", "relative_mm_counts")

relative_error_profile <- ggplot(data=relative_error.frame, aes(x=relative_positions, y=relative_mm_counts, colour=factor(""))) + geom_bar(stat="identity") + xlab("Relative position in read") + ylab("Mismatch rate") + theme(legend.position="none") + ylim(0, 1.5e-03) + scale_colour_manual(values = c("#00BADF"))


relative_error_profile




###################################

insertions<-read.table("insertion.R.txt")
names(insertions)<-c("id", "read_index", "read_base", "read_length")

c<-hist(insertions$read_index/insertions$read_length, breaks="FD")
relative_insertion.frame<-as.data.frame(cbind(seq(0.01,1,0.01), c$counts/(4127039238/100)))
names(relative_insertion.frame)<-c("relative_positions", "relative_insertion_counts")


relative_insertion_profile <- ggplot(data=relative_insertion.frame, aes(x=relative_positions, y=relative_insertion_counts, colour=factor(""))) + geom_bar(stat="identity") + xlab("Relative position in read") + ylab("Insertion rate") + theme(legend.position="none") + ylim(0, 1.5e-03) + scale_colour_manual(values = c("#FF62BE"))

relative_insertion_profile


###################################

deletions<-read.table("deletion.R.txt")
names(deletions)<-c("id", "read_index", "read_base", "read_length")

c<-hist(deletions$read_index/deletions$read_length, breaks="FD")
relative_deletion.frame<-as.data.frame(cbind(seq(0.01,1,0.01), c$counts/(4127039238/100)))
names(relative_deletion.frame)<-c("relative_positions", "relative_deletion_counts")


relative_deletion_profile <- ggplot(data=relative_deletion.frame, aes(x=relative_positions, y=relative_deletion_counts, colour=factor(""))) + geom_bar(stat="identity") + xlab("Relative position in read") + ylab("Deletion rate") + theme(legend.position="none") + ylim(0, 1.5e-03) + scale_colour_manual(values = c("#00BF74"))

relative_deletion_profile


###################################


data<-read.table("mol-combined-6.length", stringsAsFactors=F)
readlengths<-data$V1
rm(data)

lenhist<-qplot(readlengths/1000, geom="histogram", binwidth=.5, xlim=c(0,13), xlab="Moleculo read length (Kbp)", ylab="Number of reads")

###################################


# plot it all together
source('multiplot.R', chdir = TRUE)

multiplot(lenhist, relative_insertion_profile, relative_error_profile, relative_deletion_profile, cols=2)
