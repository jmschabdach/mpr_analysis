library(psych)
library(ggplot2)
library(patchwork)

gradeFn <- "/Users/youngjm/Data/slip/images/qc/2022-12-08_jms_mg_ratings.csv"
gradeDf <- read.csv(gradeFn)

slipDf <- read.csv("/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_fs.csv")
gradeDf <- gradeDf[gradeDf$subject %in% slipDf$patient_id, ]
gradeDf <- gradeDf[gradeDf$session %in% slipDf$sess_id, ]


# On a .PNG level
kap <- cohen.kappa(x=cbind(gradeDf$jenna_schabdach_grades, gradeDf$margaret_gardner_grades))

rater1 <- c(-1:2)
rater2 <- c(-1:2)
X <- c()
Y <- c()
Z <- c()
for (i in rater1){
  print(i)
  for (j in rater2){
    print(j)
    
    X <- append(X, i)
    Y <- append(Y, j)
    Z <- append(Z, dim(gradeDf[gradeDf$jenna_schabdach_grades == i & gradeDf$margaret_gardner_grades == j, ])[1])
  }
}
tmp <- log(Z, base=10)
tmp[!is.finite(tmp)] <- NA
df <- data.frame(x=X, y=Y, z=tmp)

p <- ggplot(df, aes(x, y, fill=z)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(min(df$z), max(df$z)), name="log(count(Grade 1 = Grade 2))") +
  xlab("Rater 1 QC Grades") +
  ylab("Rater 2 QC Grades") +
  labs(title=paste0("Grading Consensus Across All .PNGs (Cohen's Kappa =", format(kap$weighted.kappa, digits=3), ")")) +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

patch <- wrap_plots(p)
png(file="/Users/youngjm/Data/slip/figures/2022-12-09_qc_grading_consensus.png",
    width=800, height=600)
print(patch)
dev.off()

# Across MRI scans
# Get the average grade for each rater
rater1Avgs <- aggregate(jenna_schabdach_grades ~ subject, gradeDf, FUN=mean)
rater2Avgs <- aggregate(margaret_gardner_grades ~ subject, gradeDf, FUN=mean)

kap <- cohen.kappa(x=cbind(rater1Avgs$jenna_schabdach_grades, rater2Avgs$margaret_gardner_grades))

# Join the rater averages dataframes
jointMriDf <- merge(rater1Avgs, rater2Avgs)

p <- ggplot(jointMriDf, aes(jenna_schabdach_grades, margaret_gardner_grades)) +
  geom_jitter() +
  xlab("Rater 1 QC Grades") +
  ylab("Rater 2 QC Grades") +
  labs(title=paste0("Grading Consensus Across All MRIs (Cohen's Kappa =", format(kap$weighted.kappa, digits=3), ")")) +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

print(p)

# Add new columns with a mutation for +/- 1 for both raters
jointMriDf <- jointMriDf %>%
  mutate(jennaUsable = case_when(
    jenna_schabdach_grades >= 1.0 ~ 1.0,
    jenna_schabdach_grades < 1.0 ~ 0.0
  ))

jointMriDf <- jointMriDf %>%
  mutate(margaretUsable = case_when(
    margaret_gardner_grades >= 1.0 ~ 1.0,
    margaret_gardner_grades < 1.0 ~ 0.0
  ))

kap <- cohen.kappa(x=cbind(jointMriDf$jennaUsable, jointMriDf$margaretUsable))

p <- ggplot(jointMriDf, aes(jennaUsable, margaretUsable)) +
  geom_jitter() +
  xlab("Rater 1 QC Grades") +
  ylab("Rater 2 QC Grades") +
  labs(title=paste0("Usability Consensus Across All MRIs (Cohen's Kappa =", format(kap$weighted.kappa, digits=3), ")")) +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

print(p)

rater1 <- c(0:1)
rater2 <- c(0:1)
agree <- 0
X <- c()
Y <- c()
Z <- c()
for (i in rater1){
  print(i)
  for (j in rater2){
    print(j)
    measure <- dim(jointMriDf[jointMriDf$jennaUsable == i & jointMriDf$margaretUsable == j, ])[1]
    if (i == j) {
      agree <- agree + measure
    }
    X <- append(X, i)
    Y <- append(Y, j)
    Z <- append(Z, measure)
  }
}

df <- data.frame(x=X, y=Y, z=Z)

p <- ggplot(df, aes(x, y, fill=z)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(min(df$z), max(df$z)), name="count") +
  xlab("Rater 1 QC Grades") +
  ylab("Rater 2 QC Grades") +
  labs(title=paste0("Grading Consensus Across All .PNGs (Cohen's Kappa =", format(kap$weighted.kappa, digits=3), ")")) +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))
print(p)
