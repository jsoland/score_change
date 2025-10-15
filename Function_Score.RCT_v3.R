
#UPDATED ON 10/15/2025


library(mirt)
library(psych)

library(ggplot2)
library(dplyr)
library(reshape2)


#https://combine-australia.github.io/r-pkg-dev/functions.html

##################################################################################
##################################################################################
#CREATE A DATASET WITH 5 CATEGORIES, 4 ITEMS PER TIMEPOINT
##################################################################################
##################################################################################

## Simulation:
a <- matrix(c(0.75, 0, # specify slope parameters 
              1.30, 0, 
              1.80, 0, 
              1.65, 0,
              0, 0.75,
              0, 1.30, 
              0, 1.80, 
              0, 1.65),
            nrow=8, ncol=2, byrow=T) 
d <- matrix(c(2.13, 0.37, -0.86, -1.88, # specify intercepts 
              2.67, 1.33, 0.31, -0.67,
              3.57, 1.70, 0.07, -1.90,
              3.77, 2.32, 1.09, -0.31,
              2.13, 0.37, -0.86, -1.88,  
              2.67, 1.33, 0.31, -0.67,
              3.57, 1.70, 0.07, -1.90,
              3.77, 2.32, 1.09, -0.31),
            nrow=8, ncol=4, byrow=T) 
mu.g1 <- c(0, 0.1); mu.g2 <- c(0, 0.3) # specify group means
s.g1 <- matrix(c(1,0.9,0.9,1.4),2,2)    # specify group 1 variances 
s.g2 <- matrix(c(1,0.7,0.7,1.2),2,2)    # specify group 2 variances 

N = 1000 # sample size
set.seed(22903)

dat <- data.frame(rbind(                                         # create data frame of item responses
  simdata(a, d, N, itemtype='graded', mu = mu.g1, sigma = s.g1),  # simulate data for group 1 
  simdata(a, d, N, itemtype='graded', mu = mu.g2, sigma = s.g2))) # simulate data for group 2

dat$group <- c(rep('C', N), rep('T', N))                       # create grouping variable

dat3<-dat




##################################################################################
##################################################################################
#CREATE A DATASET WITH 4 CATEGORIES, 5 ITEMS PER TIMEPOINT
##################################################################################
##################################################################################

## Simulation:
a <- matrix(c(0.75, 0, # specify slope parameters 
              1.30, 0, 
              1.80, 0, 
              1.65, 0,
              1.65, 0,
              0, 0.75,
              0, 1.30, 
              0, 1.80, 
              0, 1.80, 
              0, 1.65),
            nrow=10, ncol=2, byrow=T) 
d <- matrix(c(2.13, 0.37, -1.88, # specify intercepts 
              2.67, 1.33, -0.67,
              3.57, 1.70, -1.90,
              2.18, 1.32, -0.11,
              3.77, 2.32, -0.31,
              2.13, 0.37, -0.86, 
              2.57, 1.78, -0.57,
              2.67, 1.33, -0.67,
              3.57, 1.70, 0.07,
              2.32, 1.09, -0.31),
            nrow=10, ncol=3, byrow=T) 
mu.g1 <- c(0, 0.1); mu.g2 <- c(0, 0.3) # specify group means
s.g1 <- matrix(c(1,0.9,0.9,1.4),2,2)    # specify group 1 variances 
s.g2 <- matrix(c(1,0.7,0.7,1.2),2,2)    # specify group 2 variances 

N = 1000 # sample size
set.seed(22903)

dat <- data.frame(rbind(                                         # create data frame of item responses
  simdata(a, d, N, itemtype='graded', mu = mu.g1, sigma = s.g1),  # simulate data for group 1 
  simdata(a, d, N, itemtype='graded', mu = mu.g2, sigma = s.g2))) # simulate data for group 2

dat$group <- c(rep('C', N), rep('T', N))                       # create grouping variable

dat4<-dat


##################################################################################
##################################################################################
#LOAD IN A DATASET
##################################################################################
##################################################################################


G1<-read.table("C:\\Users\\jgs8e\\Dropbox\\papersNWEA\\IRTpower\\z_package\\G1.csv",sep=" ")
G1$group<-"C"
G2<-read.table("C:\\Users\\jgs8e\\Dropbox\\papersNWEA\\IRTpower\\z_package\\G2.csv",sep=" ")
G2$group<-"T"

dat2<-rbind(G1,G2)


##################################################################################
##################################################################################
#FUNCTION
##################################################################################
##################################################################################

#ITEMS PER TIMEPOINT
#TIMEPOINTS
#CATEGORIES
#GROUPING

nc<-5
N<-4
sampN<-2000

score.RCT<-function(nc,N,sampN,dat) {

#########################################################
#MODEL
#########################################################

L1<- paste0(" T1 = ","1-",N) 
L2<- paste0("T2 = ",1+N,"-",N+N) 
L3<- "COV = T1*T2, T2*T2"
L4<- "CONSTRAIN = "


#----------------------------
#SLOPES
#----------------------------
vec <- rep(1:N, each = 1)
vec2 <- rep(((N+1):(2*N)), each = 1)

# Create the pattern
pattern <- paste("(", vec, ",", vec2, ", a1, a2),", sep = "")

# Combine all the patterns into a single string
L5 <- paste(pattern, collapse = " ")


#----------------------------
#INTERCEPTS
#----------------------------

# Create the pattern

for(cats in 1:(nc-1)) {

assign(paste0("LL",cats+5),paste("(", vec, ",", vec2, ",d",cats,"),", sep = ""))  
  
}


assign(paste0("LLL",5+1+nc-1),"MEAN = T2 ")



#----------------------------
#COMBINE
#----------------------------

object_names <- ls(pattern="L")


# Initialize a list to store object values
object_values <- list()

# Loop through each object name
for (obj_name in object_names) {
  # Get the value of the object and add it to the list
  object_values <- c(object_values, get(obj_name))
}

# Concatenate all values with a line break separator
mod <- paste(object_values, collapse = "\n")


#CALIBRATE MIRT
mod.est <- multipleGroup(
  dat[, 1:(N*2)],
  model = mod,
  group = dat$group,
  invariance = c(
    "free_means", "free_variances", 
    "slopes", "intercepts"), 
  itemtype='graded', method='MHRM', SE = TRUE) 
coef(mod.est, simplify=TRUE) # get item parameters

## CLEAN UNIDIM
dat.T1<-dat[,c(1:N)]
dat.T2<-dat[,c((N+1):(N*2))]
colnames(dat.T1) <- seq_along(dat.T1)
colnames(dat.T2) <- seq_along(dat.T2)
dat.U<-rbind(dat.T1,dat.T2)

## CALIBRATE UNIDIM
uni.mod <- mirt(dat.U[,(1:N)], 1, itemtype = "graded")


## EAP MIRT SCORING 
scores.EAP <- data.frame(fscores(mod.est, method='EAP'))
summary(scores.EAP)
dat$EAP1 <- scores.EAP$T1
dat$EAP2 <- scores.EAP$T2

## PV MIRT SCORING
scores.PV <- data.frame(fscores(mod.est, plausible.draws=10))
colnames(scores.PV)<-c("T1.PV1","T2.PV1",
                        "T1.PV2","T2.PV2",
                        "T1.PV3","T2.PV3",
                        "T1.PV4","T2.PV4",
                        "T1.PV5","T2.PV5",
                        "T1.PV6","T2.PV6",
                        "T1.PV7","T2.PV7",
                        "T1.PV8","T2.PV8",
                        "T1.PV9","T2.PV9",
                        "T1.PV10","T2.PV10"
                        )

#pvmods <- lapply(pv, function(x, covdata) lm(x ~ covdata$X1 + covdata$X2),
 #                covdata=covdata)

## SUM SCORING
dat$Sum1<-rowMeans(dat[,(1:N)])
dat$Sum2<-rowMeans(dat[,((N+1):(2*N))])

ms <- dat %>%
  filter(group == "C") %>%
  summarise(mean_Sum1 = mean(Sum1, na.rm = TRUE)) %>%
  pull(mean_Sum1)

vs <- dat %>%
  filter(group == "C") %>%
  summarise(variance_Sum1 = var(Sum1, na.rm = TRUE)) %>%
  pull(variance_Sum1)

dat$Sum1stand <- (dat$Sum1 - ms) / vs
dat$Sum2stand <- (dat$Sum2 - ms) / vs

#UNI EAP SCORING
scores.EAPu <- data.frame(fscores(uni.mod, method='EAP'))
EAPu.T1<-scores.EAPu[(1:sampN),]
EAPu.T2<-scores.EAPu[(sampN+1):(2*sampN),]
scores.EAPunidim<-cbind(EAPu.T1,EAPu.T2)
dat <- cbind(dat,scores.EAPunidim)


#UNI ML SCORING
scores.MLu <- data.frame(fscores(uni.mod, method='ML',max_theta=4))
MLu.T1<-scores.MLu[(1:sampN),]
MLu.T2<-scores.MLu[(sampN+1):(2*sampN),]
scores.MLunidim<-cbind(MLu.T1,MLu.T2)
dat <- cbind(dat,scores.MLunidim)
dat$MLu.T1[dat$Sum1==nc-1] <- 4
dat$MLu.T2[dat$Sum2==nc-1] <- 4
dat$MLu.T1[dat$Sum1==0] <- -4
dat$MLu.T2[dat$Sum2==0] <- -4


## PV UNI SCORING
scores.PVu <- data.frame(fscores(uni.mod, plausible.draws=10))
scores.PVuA<-scores.PVu[1:sampN,]
scores.PVuB<-scores.PVu[(sampN+1):(2*sampN),]
scores.PVuni<-cbind(scores.PVuA,scores.PVuB)

colnames(scores.PVuni)<-c("T1.PVu1",
                       "T1.PVu2",
                       "T1.PVu3",
                       "T1.PVu4",
                       "T1.PVu5",
                       "T1.PVu6",
                       "T1.PVu7",
                       "T1.PVu8",
                       "T1.PVu9",
                       "T1.PVu10",
                       "T2.PVu1",
                       "T2.PVu2",
                       "T2.PVu3",
                       "T2.PVu4",
                       "T2.PVu5",
                       "T2.PVu6",
                       "T2.PVu7",
                       "T2.PVu8",
                       "T2.PVu9",
                       "T2.PVu10"
)


dat<-cbind(dat,scores.PV,scores.PVuni)


## ML scoring 
#scores.ML <- data.frame(fscores(mod.est, method='ML',max_theta=4))
#summary(scores.ML)
#dat$ML1 <- scores.ML$T1
#dat$ML2 <- scores.ML$T2
#describeBy(dat[,c("ML1","ML2")],dat$group)


#dens1<-density(scores.EAP$T1)
#dens2<-density(scores.EAP$T2)
#plot(dens1,dens2)

##############################################
#DESCRIBE RESULTS
##############################################

print(describeBy(dat[,c("Sum1stand","Sum2stand")],dat$group))
print(describeBy(dat[,c("EAP1","EAP2")],dat$group))
print(describeBy(dat[,c("EAPu.T1","EAPu.T2")],dat$group))
print(describeBy(dat[,c("MLu.T1","MLu.T2")],dat$group))

write.table(dat,"C:\\Users\\jgs8e\\Dropbox\\papersNWEA\\IRTpower\\z_package\\dat.score.csv",
            row.names=FALSE, col.names=TRUE, sep = ",",na = "",quote = FALSE)

dat.scores<<-dat



}


score.RCT(4,5,2000,dat4)
score.RCT(5,4,2000,dat3)
score.RCT(3,5,4000,dat2)