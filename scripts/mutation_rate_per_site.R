setwd("/Users/as/Documents/ResearchProjects/vls/")

library(ggplot2)
library(tidyverse)

rate.b31 <- read.table("data-analysis/r4s-B31-ref16.txt", header = T)
ggplot(data = rate.b31, aes(x=POS, y=SCORE)) + 
  geom_rect(aes(xmin=11, xmax=31, ymin=0, ymax=8), fill='yellow', alpha=0.5) + 
  geom_rect(aes(xmin=42, xmax=53, ymin=0, ymax=8), fill='yellow', alpha=0.5) +
  geom_point() + geom_line() + theme_classic()

ggplot(data = rate.b31[95:185,], aes(x=POS+117, y=SCORE)) + 
  geom_point() + geom_line() + theme_classic() +
  scale_x_continuous(breaks = seq(220,300,20))e

rate.vlse <- read.table("data-analysis/r4s-vlsE.txt", header = T)
ggplot() + 
  geom_rect(aes(xmin=1, xmax=18, ymin=-1.2, ymax=3), fill="red", alpha=0.3)+
  geom_rect(aes(xmin=19, xmax=115, ymin=-1.2, ymax=3), fill='purple', alpha=0.3)+
  geom_rect(aes(xmin=116, xmax=120, ymin=-1.2, ymax=3), fill='green', alpha=0.3)+
  geom_rect(aes(xmin=121, xmax=150, ymin=-1.2, ymax=3), fill='blue', alpha=0.3)+
  geom_rect(aes(xmin=151, xmax=165, ymin=-1.2, ymax=3), fill='yellow', alpha=0.3)+
  geom_rect(aes(xmin=166, xmax=178, ymin=-1.2, ymax=3), fill='blue', alpha=0.3)+
  geom_point(data = rate.vlse[1:178,], aes(x=POS, y=SCORE)) + 
  geom_line(data = rate.vlse[1:178,], aes(x=POS, y=SCORE), alpha=0.3) + 
  theme_classic() +
  scale_x_continuous(breaks = seq(1,180,20))

ggplot() + 
  geom_rect(aes(xmin=179, xmax=197, ymin=-1.2, ymax=3), fill="yellow", alpha=0.2)+
  geom_rect(aes(xmin=198, xmax=204, ymin=-1.2, ymax=3), fill='blue', alpha=0.2)+
  geom_rect(aes(xmin=205, xmax=210, ymin=-1.2, ymax=3), fill='yellow', alpha=0.2)+
  geom_rect(aes(xmin=211, xmax=236, ymin=-1.2, ymax=3), fill='blue', alpha=0.2)+
  geom_rect(aes(xmin=237, xmax=253, ymin=-1.2, ymax=3), fill='yellow', alpha=0.2)+
  geom_rect(aes(xmin=254, xmax=261, ymin=-1.2, ymax=3), fill='blue', alpha=0.2)+
  geom_rect(aes(xmin=262, xmax=273, ymin=-1.2, ymax=3), fill='yellow', alpha=0.2)+
  geom_rect(aes(xmin=274, xmax=298, ymin=-1.2, ymax=3), fill='blue', alpha=0.2)+
  geom_rect(aes(xmin=299, xmax=303, ymin=-1.2, ymax=3), fill='yellow', alpha=0.2)+
  geom_rect(aes(xmin=306, xmax=310, ymin=-1.2, ymax=3), fill='green', alpha=0.2)+
  geom_rect(aes(xmin=311, xmax=356, ymin=-1.2, ymax=3), fill='purple', alpha=0.2)+
  geom_point(data = rate.vlse[179:356,], aes(x=POS, y=SCORE)) + 
  geom_line(data = rate.vlse[179:356,], aes(x=POS, y=SCORE), alpha=0.2) + 
  theme_classic() +
  scale_x_continuous(breaks = seq(180,360,20))


# mean pairwise distance between evolved seq and naive seq

diff <- read.table("aln-and-trees/second-consensus/pair-diff-2nd-evolve.tsv", header = F)
diff <- diff %>% filter(str_detect(V2, "sp")) #remove extra
diff <- diff[grep("consensus|C[0-9]", diff$V1),]
diff$V1 <- gsub("primary_consensus_30","consensus",diff$V1)
diff$V1 <- gsub("secondary_consensus_30","secondary",diff$V1)
diff$Category <- rep(c("Centroid","Consensus", "2nd_Consensus"), c(640,128,128))

ggplot(data = diff, aes(x=reorder(V1, V7, FUN = mean), y=V7, color=Category)) + 
  geom_boxplot() + theme_bw() +
  labs(x="Evolved seq", y="Mean distance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))

diff_summary <- diff %>% group_by(V1) %>% summarise(mean_diff = mean(V7), sd = sd(V7))
t.test(diff_summary$mean_diff[1:5], diff_summary$mean_diff[7:11])


#Ancestor likelihood
like <- read.table("data-analysis/ancestor-test/ancestor_likelihood_n21_VR2.txt", header = F)
#like2 <- like %>% group_by(V1) %>% summarise(mean=mean(V2))
like$V1 <- gsub("Consensus_30","Cons",like$V1)
like <- like %>% arrange(V1)
like$group <- rep(c("Centroid","Consensus","Random"), c(15,3,15))
ggplot(data=like, aes(x=reorder(V1, -V2, FUN = mean), y=V2, color=group)) + 
  geom_boxplot() + theme_bw() +
  labs(x="Evolved seq", y="Ancestral likelihood") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))

var <- read.table("aln-and-trees/second-consensus/variability-score.txt",header = T)
ggplot(data = var, aes(x=POS, y=SCORE)) + geom_line() + theme_bw() +
  labs(title = "variablity of 7 primary consensus", x="AA position", y="variability score")
