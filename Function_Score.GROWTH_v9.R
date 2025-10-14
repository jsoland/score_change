





#https://combine-australia.github.io/r-pkg-dev/functions.html




##################################################################################
##################################################################################
#LOAD IN A DATASET
##################################################################################
##################################################################################


dat2<-read.table("C:\\Users\\jgs8e\\Dropbox\\papersNWEA\\IRTpower\\z_package\\gm_mirt_v2.csv",sep=",")


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
sampN<-2774
dat2->dat
Tpoints<-4

score.GROWTH<-function(nc,N,sampN,Tpoints,dat) {
  
#########################################################  
# List of packages to be installed and loaded
#########################################################
  packages <- c("mirt", "psych", "ggplot2", "dplyr", "reshape2", "parallel")
  
  # Function to check if a package is installed
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Loop through each package and install/load it
  for (pkg in packages) {
    install_if_missing(pkg)
  }  
  
  
  cl <- makeCluster(detectCores()) #use more cores
  
  Ti<-Tpoints

  dat <- dat[!apply(dat, 1, function(x) all(is.na(x))), ]
  sampN<-dim(dat)[1]
  
#########################################################
#MODEL
#########################################################

L1<- paste0(" T1 = ","1-",N) 
L2a<- paste0("T2 = ",1+N,"-",N+N) 
L2b<- paste0("T3 = ",1+N*2,"-",N*3) 
if(Ti==4) {
  L2c<- paste0("T4 = ",1+N*3,"-",N*4)   
} else if(Ti==5) {
  L2c<- paste0("T4 = ",1+N*3,"-",N*4)  
  L2d<- paste0("T5 = ",1+N*4,"-",N*5)  
}

L3<- "COV = T1*T2, T2*T2,T3*T3,T1*T3,T2*T3"
if(Ti==4) {
  L3b<- ",T4*T4, T1*T4,T2*T4,T3*T4" 
} else if(Ti==5) {
  L3b<- ",T4*T4, T1*T4,T2*T4,T3*T4," 
  L3c<- "T5*T5, T1*T5,T2*T5,T3*T5,T4*T5" 
}

L4<- "CONSTRAIN = "


#----------------------------
#SLOPES
#----------------------------
vec <- rep(1:N, each = 1)
vec2 <- rep(((N+1):(2*N)), each = 1)
vec3 <- rep(((2*N+1):(3*N)), each = 1)

# Create the pattern
if(Ti==3) {
pattern <- paste("(", vec, ",", vec2, ",",vec3, ", a1, a2, a3),", sep = "")
} else if(Ti==4){
  vec4 <- rep(((3*N+1):(4*N)), each = 1)  
  pattern <- paste("(", vec, ",", vec2, ",",vec3,",",vec4, ", a1, a2, a3, a4),", sep = "")
} else if(Ti==5){
  vec4 <- rep(((3*N+1):(4*N)), each = 1) 
  vec5 <- rep(((4*N+1):(5*N)), each = 1) 
  pattern <- paste("(", vec, ",", vec2, ",",vec3,",",vec4,",",vec5, ", a1, a2, a3, a4, a5),", sep = "")
}

# Combine all the patterns into a single string
L5 <- paste(pattern, collapse = " ")


#----------------------------
#INTERCEPTS
#----------------------------

# Create the pattern

for(cats in 1:(nc-1)) {

if(Ti==3) {
assign(paste0("LL",cats+5),paste("(", vec, ",", vec2, ",",vec3, ",d",cats,"),", sep = ""))  
} else if(Ti==4){
  assign(paste0("LL",cats+5),paste("(", vec, ",", vec2, ",",vec3,",",vec4, ",d",cats,"),", sep = ""))    
}  else if(Ti==5){
  assign(paste0("LL",cats+5),paste("(", vec, ",", vec2, ",",vec3,",",vec4,",",vec5, ",d",cats,"),", sep = ""))  
}
  
}



if(Ti==3){
assign(paste0("LLL",6+1+nc-1),"MEAN = T2, T3 ")
} else if(Ti==4){
  assign(paste0("LLL",6+1+nc-1),"MEAN = T2, T3, T4 ")
} else if (Ti==5){
  assign(paste0("LLL",6+1+nc-1),"MEAN = T2, T3, T4, T5 ")  
}


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


#----------------------------------------------------
#CALIBRATE MIRT
#----------------------------------------------------
mod.est <- mirt(
  dat[, 1:(N*Ti)],
  model = mod,
  itemtype='graded', method='MHRM', SE = TRUE) 
coef(mod.est, simplify=TRUE) # get item parameters

#----------------------------------------------------
## CLEAN UNIDIM
#----------------------------------------------------
dat.T1<-dat[,c(1:N)]
dat.T2<-dat[,c((N+1):(N*2))]
dat.T3<-dat[,c((2*N+1):(N*3))]
if(Ti==4){
  dat.T4<-dat[,c((3*N+1):(N*4))]  
} else if (Ti==5) {
  dat.T4<-dat[,c((3*N+1):(N*4))] 
  dat.T5<-dat[,c((4*N+1):(N*5))] 
  
}

#get list of time-specific datasets
all_objects <- ls()
dat.T.names <- grep("^dat\\.T", all_objects, value = TRUE)
dat.T.list <- mget(dat.T.names)

# Change column names of each dataframe in the list
for(tl in 1:Ti) {
  colnames(dat.T.list[[tl]]) <- seq_along(dat.T.list[[tl]])
}


dat.U<-do.call(rbind, dat.T.list)

## CALIBRATE UNIDIM
dat.Unm <- dat.U[!apply(dat.U, 1, function(x) all(is.na(x))), ]
uni.mod <- mirt(dat.U[,(1:N)], 1, itemtype = "graded")



#########################################################
#SCORING
#########################################################


#----------------------------------------------------
## EAP MIRT SCORING 
#----------------------------------------------------
scores.EAP <- data.frame(fscores(mod.est, method='EAP',QMC=TRUE))
summary(scores.EAP)
colnames(scores.EAP)<-c("EAP.T1","EAP.T2","EAP.T3","EAP.T4")
dat<-cbind(dat,scores.EAP)

#----------------------------------------------------
## SUM SCORING
#----------------------------------------------------

dat$Sum1<-rowMeans(dat[,(1:N)])
for (Su in 2:Ti) {
  dat <- dat %>%
    mutate(!!paste0("Sum", Su) := rowMeans(dat[, ((Su-1)*N+1):(Su*N)]))
}

for (Su in 1:Ti) {
  # Calculate the start and end columns for the current time point
  start_col <- (Su - 1) * N + 1
  end_col <- Su * N
  
  # Create a new column with the row means
  dat[[paste0("Sum", Su)]] <- rowMeans(dat[, start_col:end_col])
}


#STANDARDIZE
ms <- mean(dat$Sum1,na.rm = TRUE)
vs <- sd(dat$Sum1,na.rm = TRUE)


for (Su in 1:Ti) {
  # Create a new standardized column
  dat[[paste0("Sum", Su, "stand")]] <- (dat[[paste0("Sum", Su)]] - ms) / vs
}



#----------------------------------------------------
#UNI EAP SCORING
#----------------------------------------------------
scores.EAPu <- data.frame(fscores(uni.mod, method='EAP'))

EAPu.T1<-scores.EAPu[(1:sampN),]
for(TT in 2:Ti){
assign(paste0("EAPu.T",TT),scores.EAPu[((TT-1)*sampN+1):(TT*sampN),])  
}


for (TT in 1:Ti) {
  col_name <- paste0("EAPu.T", TT)
  dat[[col_name]] <- get(col_name)  # Fetch the object using get()
}

#----------------------------------------------------
#UNI ML SCORING
#----------------------------------------------------
scores.MLu <- data.frame(fscores(uni.mod, method='ML',max_theta=4))

# Initialize a list to store MLu.TT dataframes
MLu_list <- list()

# Assign MLu.T1 to MLu.Ti with correct column names
for (TT in 1:Ti) {
  start_idx <- (TT - 1) * sampN + 1
  end_idx <- TT * sampN
  MLu_list[[TT]] <- as.data.frame(scores.MLu[start_idx:end_idx, ])
  colnames(MLu_list[[TT]]) <- paste0("MLu.T", TT)
}

# Assign scores to dataframe dat with proper column names
for (TT in 1:Ti) {
  col_name <- paste0("MLu.T", TT)
  dat[[col_name]] <- unlist(MLu_list[[TT]])  # Unlist the dataframe before assigning
}

ml_cols <- grep("ML", names(dat), value = TRUE)

# Change values to 4 if Inf, -4 if -Inf
for (col_name in ml_cols) {
  dat[[col_name]][dat[[col_name]] == Inf] <- 4
  dat[[col_name]][dat[[col_name]] == -Inf] <- -4
}

#----------------------------------------------------
## PV MIRT SCORING
#----------------------------------------------------
scores.PV <- data.frame(fscores(mod.est, plausible.draws=5,QMC=TRUE))

#RENAME PV COLUMNS
new_col_names <- c()
for (pv in 1:5) {
  for (t in 1:Ti) {
    new_col_names <- c(new_col_names, paste0("T", t, ".pv", pv))
  }
}


# Assign new column names to the dataframe
colnames(scores.PV) <- new_col_names
dat<-cbind(dat,scores.PV)


#----------------------------------------------------
## PV UNI SCORING
#----------------------------------------------------
scores.PVu <- data.frame(fscores(uni.mod, plausible.draws=5))

PVu.T1<-scores.PVu[(1:sampN),]
for(TT in 2:Ti){
  assign(paste0("PVu.T",TT),scores.PVu[((TT-1)*sampN+1):(TT*sampN),])  
}



# Create a list of dataframes
df_list <- lapply(1:Ti, function(i) get(paste0("PVu.T", i)))

# Combine the dataframes using cbind
combined_df <- do.call(cbind, df_list)

# Function to generate the original combined column names
generate_colnames_original <- function(Ti, ncol_each) {
  colnames_combined <- c()
  for (i in 1:Ti) {
    colnames_combined <- c(colnames_combined, paste0("T", i, ".PVu", 1:ncol_each))
  }
  return(colnames_combined)
}

# Assume all dataframes have the same number of columns
ncol_each <- ncol(df_list[[1]])

# Generate original column names for the combined dataframe
original_colnames <- generate_colnames_original(Ti, ncol_each)

# Assign these column names to the combined dataframe
colnames(combined_df) <- original_colnames

# Print the combined dataframe to verify
#print(combined_df)

# Function to generate the new column names in the desired order
generate_colnames_by_pv <- function(Ti, ncol_each) {
  colnames_combined <- c()
  for (pv in 1:ncol_each) {
    colnames_combined <- c(colnames_combined, paste0("T", 1:Ti, ".PVu", pv))
  }
  return(colnames_combined)
}

# Generate new column names in the desired order
new_colnames_by_pv <- generate_colnames_by_pv(Ti, ncol_each)

# Verify the generated column names match the existing ones before reordering
#print(new_colnames_by_pv %in% colnames(combined_df))

# Reorder the columns in the combined dataframe according to the new column names
combined_df <- combined_df[, new_colnames_by_pv]

dat<-cbind(dat,combined_df)





##############################################
#DESCRIBE RESULTS
##############################################


#--------------------------------------------------
#SUM
#--------------------------------------------------
# Find columns with "sum" in their names
sum_cols <- grep("Sum", colnames(dat), value = TRUE)

# Subset the dataframe to include only these columns
sum_vars <- dat[, sum_cols]

# Summary or further operations on sum_vars if needed
print("-----------------SUM SCORE SUMMARY STATS-----------------")
print(summary(sum_vars))


#--------------------------------------------------
#EAP
#--------------------------------------------------
# Find columns with "EAP" in their names
eap_cols <- grep("EAP", colnames(dat), value = TRUE)

# Subset the dataframe to include only these columns
eap_vars <- dat[, eap_cols]

# Summary or further operations on eap_vars if needed
print("-----------------EAP (MIRT & UNIDIM.) SCORE SUMMARY STATS-----------------")
print(summary(eap_vars))

#--------------------------------------------------
#MLE
#--------------------------------------------------
# Example to troubleshoot and calculate column means
ml_cols <- grep("ML", colnames(dat), value = TRUE)

# Subset the dataframe to include only these columns
ml_vars <- dat[, ml_cols]

# Summary or further operations on eap_vars if needed
print("-----------------MLE UNIDIM. SCORE SUMMARY STATS-----------------")
print(summary(ml_vars))


#--------------------------------------------------
#PV UNI
#--------------------------------------------------

# Initialize an empty list to store the dataframes
PVM_list <- list()

# Loop through each value of T
for (i in 1:N) {
  T_value <- paste0("T", i)
  
  # Select columns that include both "pv" and the current T value
  selected_cols <- grep(paste0("PVu.*", T_value, "|", T_value, ".*PVu"), names(dat), value = TRUE)
  
  # Create a new dataframe for the current T value
  PVM_list[[T_value]] <- dat[, selected_cols]
}

# You can now access each dataframe with PVM_list[["T1"]], PVM_list[["T2"]], etc.
# For example, to view PVM1:
#head(PVM_list[["T1"]])



# Initialize a new dataframe to store the results with the correct number of rows
# Initialize with NA values
first_df <- PVM_list[[1]]
PVM_summary <- data.frame(matrix(NA, nrow = nrow(first_df), ncol = 0))

# Loop through each dataset in PVM_list
for (T_value in names(PVM_list)) {
  # Extract the current dataset
  current_df <- PVM_list[[T_value]]
  
  # Check if the number of rows in the current dataframe matches
  if (nrow(current_df) != nrow(PVM_summary)) {
    stop("Row count mismatch between datasets")
  }
  
  # Calculate row-wise mean and standard deviation
  row_mean <- rowMeans(current_df, na.rm = TRUE)
  row_sd <- apply(current_df, 1, sd, na.rm = TRUE)
  
  # Add these as new columns to the summary dataframe
  PVM_summary[[paste0(T_value, "_mean")]] <- row_mean
  PVM_summary[[paste0(T_value, "_sd")]] <- row_sd
}

# View the resulting summary dataframe
#head(PVM_summary)

print("-----------------PV UNIDIM. SCORE SUMMARY STATS-----------------")
print(summary(PVM_summary))



#--------------------------------------------------
#PV MIRT
#--------------------------------------------------

# Initialize an empty list to store the dataframes
PVM_list <- list()

# Loop through each value of T
for (i in 1:N) {
  T_value <- paste0("T", i)
  
  # Select columns that include both "pv" and the current T value
  selected_cols <- grep(paste0("pv.*", T_value, "|", T_value, ".*pv"), names(dat), value = TRUE)
  
  # Create a new dataframe for the current T value
  PVM_list[[T_value]] <- dat[, selected_cols]
}

# You can now access each dataframe with PVM_list[["T1"]], PVM_list[["T2"]], etc.
# For example, to view PVM1:
#head(PVM_list[["T1"]])



# Initialize a new dataframe to store the results with the correct number of rows
# Initialize with NA values
first_df <- PVM_list[[1]]
PVM_summary <- data.frame(matrix(NA, nrow = nrow(first_df), ncol = 0))

# Loop through each dataset in PVM_list
for (T_value in names(PVM_list)) {
  # Extract the current dataset
  current_df <- PVM_list[[T_value]]
  
  # Check if the number of rows in the current dataframe matches
  if (nrow(current_df) != nrow(PVM_summary)) {
    stop("Row count mismatch between datasets")
  }
  
  # Calculate row-wise mean and standard deviation
  row_mean <- rowMeans(current_df, na.rm = TRUE)
  row_sd <- apply(current_df, 1, sd, na.rm = TRUE)
  
  # Add these as new columns to the summary dataframe
  PVM_summary[[paste0(T_value, "_mean")]] <- row_mean
  PVM_summary[[paste0(T_value, "_sd")]] <- row_sd
}

# View the resulting summary dataframe
#head(PVM_summary)

print("-----------------PV MIRT SCORE SUMMARY STATS-----------------")
print(summary(PVM_summary))


#--------------------------------------------------
#SAVE OUT DATA
#--------------------------------------------------


dat.scores<<-dat



}

nc<-5
N<-4
sampN<-2774
dat2->dat
Ti<-4

score.GROWTH(5,4,2774,4,dat2)
