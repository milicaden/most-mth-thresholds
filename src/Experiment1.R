source("./Most-MTH-functions.R")
##########################################################################################
# Set the folder
##########################################################################################
Folder = "../data/"
theme_set(theme_bw())
##########################################################################################
#
# Load and organize the data
#
##########################################################################################

rm(d)
d = load.data(Folder, "experiment-one")

# Load the design files
design1 <- read.csv(paste0(Folder, "Experiment1-most-design.csv"), header = TRUE, sep = ",")
design2 <- read.csv(paste0(Folder, "Experiment1-mth-design.csv"), header = TRUE, sep = ",")

# Merge the results with the design file
d1 <- merge(subset(d, (Exp == 1| Exp == 2)), design1, by="Rtype")
d2 <- merge(subset(d, (Exp == 3| Exp == 4)), design2, by="Rtype")
d <- rbind(d1, d2)

# Check out the number of participants in each experiment
part_num = part_n(d)
##########################################################################################
#
# Exclusions and data organization
#
##########################################################################################
### Exclude non-native speakers of English :
d = english_excl(d)

###Exclude based on performance in the control conditions
d = control_excl(d)

### Exclude responses whose RT takes more/less than 2 SD from the mean logRT within each experiment
d = RT_excl(d)

# Check out the number of participants in each experiment post exclusions
part_num = part_n(d)

##########################################################################################
#
# Determine the threshold for most for each participant: method 1 (logistic function fitting)
#
##########################################################################################
# We fit the logistic function separately for upward and downward versions (they are mirror images of each other)

#simple logistic function fitting
fitted1 <- dlply(subset(d, Exp == 1 & Quantifier == "Most"), .(ProlificID), model1)
fitted1success=fitted1[(which(sapply(fitted1,is.numeric),arr.ind=TRUE))]

fitted2 <- dlply(subset(d, Exp == 2 & Quantifier == "Most"), .(ProlificID), model2)
fitted2success=fitted2[(which(sapply(fitted2,is.numeric), arr.ind=TRUE))]

fitted3 <- dlply(subset(d, Exp == 3 & Quantifier == "More-than-half"), .(ProlificID), model1)
fitted3success=fitted3[(which(sapply(fitted3,is.numeric),arr.ind=TRUE))]

fitted4 <- dlply(subset(d, Exp == 4 & Quantifier == "More-than-half"), .(ProlificID), model2)
fitted4success=fitted4[(which(sapply(fitted4,is.numeric), arr.ind=TRUE))]

thresholds <- data.frame(Reduce(rbind, c(fitted1success, fitted2success, fitted3success, fitted4success)))
thresholds$ProlificID=names(c(fitted1success, fitted2success, fitted3success, fitted4success))

dth <- merge(d, thresholds, by=c("ProlificID"))


##########################################################################################
#
# Evaluate thresholds for those participants for whom the method 1 failed: method 2
#
##########################################################################################
##### Experiment 1 and 3 (positive versions) ####
mostdfexp1 <- subset(d, ((Exp == 1 & Quantifier == "Most")|(Exp == 3 & Quantifier == "More-than-half")) & !(ProlificID %in% dth$ProlificID))
##### Experiment 2 and 4 (negative versions) ####
mostdfexp2 <- subset(d, ((Exp == 2 & Quantifier == "Most")|(Exp == 4 & Quantifier == "More-than-half")) & !(ProlificID %in% dth$ProlificID))




# Estimate thresholds for experimental versions when monotonicity is upward
thresholds1 = estimate_thresh(mostdfexp1, "upward")

# Estimate thresholds for experimental versions when monotonicity is upward
thresholds2 = estimate_thresh(mostdfexp2, "downward")

# bind and sort over experiments
allthresholds <- rbind(thresholds1, thresholds2)
allthresholds <- merge(allthresholds, unique(d[c("ProlificID", "Exp")]), by=c("ProlificID"))

##########################################################################################
#
# Put together thresholds from methods 1 and 2
#
##########################################################################################
#merge with thresholds estimated by logistic function fitting
finalthresholds <- rbind(unique(dth[,c("ProlificID", "xmid", "Exp")]), allthresholds)

# Exclude those people for whom smth went wrong (threshold 0 or lower or 100 or higher)
finalthresholds <-subset(finalthresholds, xmid >0 & xmid < 100)
part_num = part_n(finalthresholds)

finalthresholds$Exp <- as.factor(finalthresholds$Exp)

# Rename experiments
finalthresholds$Condition <- ifelse(finalthresholds$Exp == "1", "most-pos", ifelse(finalthresholds$Exp == "2", "most-neg", ifelse(finalthresholds$Exp == "3", "mth-pos", "mth-neg")))

##########################################################################################
#
# Plot final thresholds
#
##########################################################################################

# Mean threshold and SD per version
mu <- ddply(finalthresholds, "Condition", summarise, grp.mean=mean(xmid), grp.sd = sd(xmid))

#Plot threshold distribution
finalthresholds$Quantifier <- ifelse(finalthresholds$Exp == "1"|finalthresholds$Exp == "2", "Most", "More_than_half")
finalthresholds$Monotonicity <- ifelse(finalthresholds$Exp == "1"|finalthresholds$Exp == "3", "Upward", "Downward")

min(finalthresholds$xmid)
max(finalthresholds$xmid)

#png("Experiment1-thresholds.png", width = 160, height = 100, units='mm', res = 300)
ggplot(data=finalthresholds, aes(finalthresholds$xmid)) + 
  geom_histogram(position = "identity", breaks=seq(25, 85, by =10), 
                 alpha = .4)+
  facet_grid(Monotonicity~Quantifier)+
  labs(x="Thresholds in Experiment 1")
#dev.off()


#T-tests
t.test(subset(finalthresholds, Exp == "1")$xmid, subset(finalthresholds, Exp =="3")$xmid)
t.test(subset(finalthresholds, Exp == "2")$xmid, subset(finalthresholds, Exp =="4")$xmid)

#Two-way anova
anovadf <- finalthresholds
anovadf$Monotonicity = as.factor(anovadf$Monotonicity)
anovadf$Quantifier = as.factor(anovadf$Quantifier)
fit = aov(xmid ~ Monotonicity*Quantifier, data = anovadf)
summary(fit)

#Levene's test
leveneTest(xmid ~ Monotonicity*Quantifier, data = anovadf)# Levene's test significant.

#Two-way aligned rank anova
aranova = art(xmid ~ Monotonicity*Quantifier, data = anovadf)
anova(aranova)

#Power analysis for anova (done online with Shiny app http://shiny.ieis.tue.nl/anova_power/):
#Already with 40 participants per group we should be fine.

##########################################################################################
#
# Vagueness: is there a higher uncertainty around the threshold for most than for mth?
#
##########################################################################################
dvague = vague_prep(finalthresholds, d)


mod1 <- glmer(correctansw ~ distance*Quantifier + (1 | ProlificID), data = dvague, family = binomial)
summary(mod1)
mod2 <- glmer(correctansw ~ distance+Quantifier + (1 | ProlificID), data = dvague, family = binomial)
summary(mod2)
anova(mod1, mod2)

#Compute probabilities for the plot
probabilities_correct = compute_prob(mod1)

#Plot#
#png("Experiment1-borderline.png", width = 160, height = 100, units='mm', res = 300)
ggplot(probabilities_correct, aes(x=x, y=prob, color=Quantifier)) + 
  geom_line(lwd=2) + 
  labs(x="Absolute distance from threshold", y="P(response is correct)")
#dev.off()
