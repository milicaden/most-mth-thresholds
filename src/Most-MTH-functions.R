library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(sciplot)
library(nlme) 
library(lme4)
library(purrr)
library("lsr")
library(tidyr)
library("lsr")
library("DescTools")
library("rstatix")
library("ARTool")
library(car)
options(scipen = 999)

length.unique <- function(x){length(unique(x))}

# Load the four results files whose names start with 'experiment' and are in 'Folder'
load.data = function(Folder, experiment){
  
  for (filenumber in c(1,2,3,4)) {
    
    
    # Load the result files
    filename <- paste(Folder, experiment, filenumber, ".txt", sep="")
    all <- readLines(filename)
    
    # Extract the questionnaire data
    questionnaire <- subset(all, grepl(",form,", all))
    questionnaire <- read.csv(textConnection(questionnaire), header = F)
    names(questionnaire) <- c("timeresults", "ip", "controller", "itemnumber", "elementnumber", "type", "group", "fieldname", "fieldvalue")
    questionnaire <- dcast(questionnaire, timeresults + ip + controller + itemnumber + elementnumber + type + group ~ fieldname, value.var="fieldvalue")
    
    # Extract the responses
    results <- subset(all, grepl("Question2", all))
    results <- read.csv(textConnection(results), header = F)	
    names(results) <- c("timeresults", "ip", "Rcontroller", "Ritemnumber", "Relementnumber", "Rtype", "Rgroup", "question", "answer", "correct", "RT")
    results <- subset(results, !is.na(RT))
    
    # Merge the data from the questionnaire with responses
    E <- merge(questionnaire, results)
    E$Exp <- filenumber
    E$comments <- NULL
    E$answer <- factor(E$answer)
    
    if (exists("d")) {
      d <- rbind(d, E)
    }
    else
    { d <- E }
    
  }
  
  d$ProlificID <- as.factor(d$ProlificID)
  d$timeresults <- as.factor(as.character(d$timeresults))
  
  #Recode responses numerically
  d$answernum <- ifelse(d$answer== "True", 1, 0)
  
  
  return(d)
}

#Check the number of participant in each experimental version
part_n = function(df){
  part_num = df %>%
    group_by(Exp) %>%
    summarise(part_n = length.unique(ProlificID))
  return(part_num)
}

# Exclude non-native speakers of English
english_excl = function(df){
  df$NativeEnglish <- toupper(df$language) %in% c("ENGLISH", "ENG", "EN", "ENGLISH ", "UNITED STATES", "ENLISH", "ENGLISH (US)", "ENLGISH", "UNITED KINGDOM", "BRITISH", "EMGLISH", " ENGLISH", "ENGLIS", "CANADA", "ENLGLISH", "NGLISH", "EN_US", "ENGLISH LANGUAGE ")
  print(length.unique(df$ProlificID[!(df$NativeEnglish)])) #print how many participants were excluded for not being native speakers
  df <- subset(df, NativeEnglish)
  return(df)
}

# Exclude based on performance on uncontroversial controls

control_excl = function(d){
  ###Exclude based on responses in the control conditions
  #For positive sentences, we expect people to say true with 'some' at least for 0 < X < 50 
  excl1 <- subset(d, ((Quantifier == "Some" & Percentage < 50)|(Quantifier == "All")))
  excl1 <- mutate(excl1, expected = ifelse((Quantifier == "Some" & Percentage > 0)|(Quantifier == "All" & Percentage == 100), 1, 0))
  excl1 <- mutate(excl1, correct = ifelse(excl1$answernum == excl1$expected, 1, 0))

  #For negative sentences, we expect people to say true with ' not all' at least for 50 > X > 100 
  excl2 <- subset(d, ((Quantifier == "All-Neg" & Percentage > 50)|(Quantifier == "Some-Neg")))
  excl2 <- mutate(excl2, expected = ifelse((Quantifier == "Some-Neg" & Percentage == 0)|(Quantifier == "All-Neg" & Percentage < 100), 1, 0))
  excl2 <- mutate(excl2, correct = ifelse(excl2$answernum == excl2$expected, 1, 0))

  excl <- rbind(excl1, excl2)

  # Calculate for each participant his score on controls
  excl.score <- ddply(excl, c("ProlificID"), 
                      function(df)c(CTLmean=mean(df$correct)))

  # Keep only the participants who scored at least 80% on controls
  d <- subset(d, ProlificID %in% subset(excl.score, excl.score$CTLmean>= 0.8)$ProlificID)
  return(d)
}

### Exclude responses whose RT takes more/less than 2 SD from the mean logRT within each experiment

RT_excl = function(d){
  d$logRT <-log(d$RT) #Take the log from response time.
  dRT = d %>% group_by(Exp) %>%  summarise(
    RTbound1 = mean(logRT)+sd(logRT)*c(-2),
    RTbound2 = mean(logRT)+sd(logRT)*c(2)
  )

  d = merge(d, dRT, by = c("Exp"))
  d=subset(d, logRT>RTbound1 & logRT<RTbound2)
  return(d)
}

# Function for upward versions, start scale = 4 (simple logistic function)
model1 <- function(x){
  try(coef(nls(answernum ~ 1/(1+exp((xmid-Percentage)/scal)), data=x, 
               start=list(xmid=50,scal=4))))
}

# Function for downward versions, start scale = -4 (simple logistic function)
model2 <- function(x){
  try(coef(nls(answernum ~ 1/(1+exp((Percentage-xmid)/scal)), data=x, 
               start=list(xmid=50,scal=4))))
}

# Estimate thresholds: method 2

# Calculate per participant per percentage tested:
#mean response on smaller percentages  
#mean response or bigger percentages 

estimate_thresh = function(df, monotonicity){
  participants <- c()
  percentages <- c()
  downwardmeans <- c()
  upwardmeans <- c()
  for(i in unique(df$ProlificID)){
    for(k in unique(df$Percentage)){
      j <- subset(df, ProlificID == i & Percentage < k)
      u <- subset(df, ProlificID == i & Percentage > k)
      down <- mean(j$answernum)
      up <- mean(u$answernum)
      participants <- c(participants, i)
      percentages <- c(percentages, k)
      downwardmeans <- c(downwardmeans, down)
      upwardmeans <- c(upwardmeans, up)
    }
  }

  experiment <- data.frame(participants, percentages, downwardmeans, upwardmeans)
  names(experiment)[names(experiment) == "participants"] <- "ProlificID"
  names(experiment)[names(experiment) == "percentages"] <- "Percentage"

  experiment<- merge(experiment, df, by= c("ProlificID", "Percentage"))
  experiment[is.na(experiment)] <- 0

  #estimate thresholds based on monotonicity
  thresholds <- c()
  participants <- c()
  for(i in unique(experiment$ProlificID)){
    temp <- subset(experiment, ProlificID == i)
    if(monotonicity == "upward"){
      threshold <- (max(temp$Percentage[which(temp$answernum == 0 & temp$downwardmeans <0.2)]) + min(temp$Percentage[which(temp$answernum == 1 & temp$upwardmeans >0.8)]))/2 
    } else if(monotonicity == "downward"){
      threshold <- (max(temp$Percentage[which(temp$answernum == 1 & temp$downwardmeans > 0.8)]) + min(temp$Percentage[which(temp$answernum == 0 & temp$upwardmeans <0.2)]))/2 
    }
    thresholds <- c(thresholds, threshold)
    participants <- c(participants, i)
  }

  thresholds <- data.frame(thresholds, participants)
  names(thresholds)[names(thresholds) == "thresholds"] <- "xmid"
  names(thresholds)[names(thresholds) == "participants"] <- "ProlificID"
  return(thresholds)
}

# Prepare the data frame for vagueness analysis

vague_prep = function(finalthresholds, d){
  finalthresholds$Quantifier = NULL
  finalthresholds$Monotonicity = NULL
  dvague = merge(finalthresholds, subset(d, Quantifier == "Most"|Quantifier == "More-than-half"), by = c("ProlificID", "Exp"))
  dvague = dvague[c("ProlificID", "Exp", "xmid", "Condition", "answer", "Quantifier", "Percentage", "answernum")]

  #Subset to cases +/- 25 around threshold 
  dvague = subset(dvague, Percentage >= xmid-25 & Percentage <= xmid + 25)

  #Add the measure of distance
  dvague$distance = abs(dvague$xmid-dvague$Percentage)

  #Code for correct vs. not for positive
  dvaguepos = subset(dvague, Condition == 'mth-pos'|Condition == 'most-pos')
  dvaguepos$correctansw = ifelse((dvaguepos$Percentage > dvaguepos$xmid & dvaguepos$answernum == 1)|(dvaguepos$Percentage < dvaguepos$xmid & dvaguepos$answernum == 0), 1, 0)

  #Code for correct vs. not for positive
  dvagueneg = subset(dvague, Condition == 'mth-neg'|Condition == 'most-neg')
  dvagueneg$correctansw = ifelse((dvagueneg$Percentage > dvagueneg$xmid & dvagueneg$answernum == 0)|(dvagueneg$Percentage < dvagueneg$xmid & dvagueneg$answernum == 1), 1, 0)

  #Positive and negative
  dvaguetogether = rbind(dvaguepos, dvagueneg)
  dvaguetogether$version = ifelse(dvaguetogether$Condition == "most-pos"|dvaguetogether$Condition == "mth-pos", "positive", "negative")
  dvaguetogether$version = as.factor(dvaguetogether$version)

  # Model
  dvaguetogether[] <- lapply(dvaguetogether, function(x) if(is.factor(x)) factor(x) else x)

  #Sum contrast coding of the relevant variables
  contrasts(dvaguetogether$Quantifier) <- contr.sum(2)
  contrasts(dvaguetogether$version) <- contr.sum(2)
  return(dvaguetogether)
}


#compute probabilities for the interaction plot
compute_prob = function(model){
  ###Plot
  #Extract model's coefficients
  b0 <- fixef(model)[1] # intercept
  distance <- fixef(model)[2]
  quant1 <- fixef(model)[3]
  distance.quant1 <- fixef(model)[4]
  
  #X-range
  x_range <- seq(from=0, to=25, by=.01)
  
  #Logits (quantifier is contrast coded, most = 1, more than half = -1)
  mth_logits <- b0 + 
    distance*x_range + 
    quant1*(-1) +
    distance.quant1*x_range*(-1) 
  
  most_logits <- b0 + 
    distance*x_range + 
    quant1*1+
    distance.quant1*x_range*1
  
  #Probabilities
  mth_probs <- exp(mth_logits)/(1 + exp(mth_logits))
  most_probs <- exp(most_logits)/(1 + exp(most_logits))
  
  plot.data <- data.frame(More_than_half=mth_probs, Most=most_probs, x=x_range)
  plot.data <- gather(plot.data, key=Quantifier, value=prob, c(More_than_half, Most))
  return(plot.data)
}