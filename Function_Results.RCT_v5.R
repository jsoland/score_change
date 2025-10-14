




dat.s<-dat.scores
plots<-"True"
ests<-"True"


results.RCT<-function(score.file,plots="True",ests="True") {
  
  # Vector of required packages
  packages <- c("tidyverse", "lmtest", "sandwich","mice","tidyr","dplyr","mitools")
  
  # Install any missing packages
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
  
  if (length(to_install) > 0) {
    install.packages(to_install)
  }
  
  # Load the packages
  lapply(packages, library, character.only = TRUE)

#DATA FILE    
dat.s<-score.file  

#####################################################################
# PLOTS (full replacement)
#####################################################################

if (plots == "True") {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(stringr)
    library(tidyr)
  })
  
  # ---------- helper: robust column finder ----------
  find_cols <- function(x, patterns) {
    unique(unlist(lapply(patterns, function(p)
      grep(p, names(x), value = TRUE, ignore.case = TRUE)
    )))
  }
  
  # ---------- detect PV columns (matches your screenshot) ----------
  # MIRT PVs: T1.PV1 ... T1.PV10 (exclude PVu); same for T2
  mirt_t1 <- grep("^T1[._]PV(?!u)\\d+$", names(dat.s), value = TRUE, perl = TRUE)
  mirt_t2 <- grep("^T2[._]PV(?!u)\\d+$", names(dat.s), value = TRUE, perl = TRUE)
  
  # UNIDIM PVs: T1.PVu1 ... T1.PVu10; same for T2
  uni_t1  <- grep("^T1[._]PVu\\d+$", names(dat.s), value = TRUE)
  uni_t2  <- grep("^T2[._]PVu\\d+$", names(dat.s), value = TRUE)
  
  # (Optional) diagnostics
  message(sprintf("MIRT T1: %d -> %s", length(mirt_t1), paste(head(mirt_t1, 10), collapse = ", ")))
  message(sprintf("MIRT T2: %d -> %s", length(mirt_t2), paste(head(mirt_t2, 10), collapse = ", ")))
  message(sprintf("UNI  T1: %d -> %s", length(uni_t1),  paste(head(uni_t1,  10), collapse = ", ")))
  message(sprintf("UNI  T2: %d -> %s", length(uni_t2),  paste(head(uni_t2,  10), collapse = ", ")))
  
  
  # ---------- compute rowwise PV means and gains ----------
  dat.s <- dat.s %>%
    mutate(
      MIRT_PVmean_T1 = if (length(mirt_t1) > 0) rowMeans(across(all_of(mirt_t1)), na.rm = TRUE) else NA_real_,
      MIRT_PVmean_T2 = if (length(mirt_t2) > 0) rowMeans(across(all_of(mirt_t2)), na.rm = TRUE) else NA_real_,
      UNIDIM_PVmean_T1 = if (length(uni_t1) > 0) rowMeans(across(all_of(uni_t1)), na.rm = TRUE) else NA_real_,
      UNIDIM_PVmean_T2 = if (length(uni_t2) > 0) rowMeans(across(all_of(uni_t2)), na.rm = TRUE) else NA_real_,
      MIRT_PVmean_gain   = MIRT_PVmean_T2   - MIRT_PVmean_T1,
      UNIDIM_PVmean_gain = UNIDIM_PVmean_T2 - UNIDIM_PVmean_T1
    )
  
  # ---------- your existing gains ----------
  dat.s$MIRT_EAP    <- dat.s$EAP2      - dat.s$EAP1
  dat.s$UNIDIM_EAP  <- dat.s$EAPu.T2   - dat.s$EAPu.T1
  dat.s$SUM_SCORE   <- dat.s$Sum2stand - dat.s$Sum1stand
  dat.s$UNIDIM_ML   <- dat.s$MLu.T2    - dat.s$MLu.T1
  
  # ---------- long data for boxplot ----------
  data_long <- melt(
    dat.s,
    id.vars = "group",
    measure.vars = c("SUM_SCORE","UNIDIM_ML","UNIDIM_EAP","MIRT_EAP",
                     "MIRT_PVmean_gain","UNIDIM_PVmean_gain"),
    variable.name = "variable",
    value.name    = "value"
  )
  
  # ---------- BOXLOT with MEAN overlays (points + labels) ----------
  p1 <- ggplot(data_long, aes(x = variable, y = value, fill = group)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.alpha = 0.25) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 0.75),
                 shape = 23, size = 2) +
    stat_summary(fun = mean, geom = "text",
                 aes(label = round(..y.., 2)),
                 position = position_dodge(width = 0.75),
                 vjust = -0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Gain Scores by Group",
         x = "Score Type", y = "Gain Scores") +
    theme_minimal()
  
  print(p1); plot.1 <<- p1
  
  # ---------- Difference-in-Differences (T - C) ----------
  difference_summary <- dat.s %>%
    group_by(group) %>%
    summarise(
      mean_diffEAP    = mean(MIRT_EAP,           na.rm = TRUE),
      mean_diffEAPu   = mean(UNIDIM_EAP,         na.rm = TRUE),
      mean_diffML     = mean(UNIDIM_ML,          na.rm = TRUE),
      mean_diffSUM    = mean(SUM_SCORE,          na.rm = TRUE),
      mean_diffMIRTpv = mean(MIRT_PVmean_gain,   na.rm = TRUE),
      mean_diffUNIpv  = mean(UNIDIM_PVmean_gain, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Safe extractors for T and C (return NA if missing)
  getTC <- function(vec, grp, g) if (any(grp == g)) vec[grp == g][1] else NA_real_
  T_EAP    <- getTC(difference_summary$mean_diffEAP,    difference_summary$group, "T")
  C_EAP    <- getTC(difference_summary$mean_diffEAP,    difference_summary$group, "C")
  T_EAPu   <- getTC(difference_summary$mean_diffEAPu,   difference_summary$group, "T")
  C_EAPu   <- getTC(difference_summary$mean_diffEAPu,   difference_summary$group, "C")
  T_ML     <- getTC(difference_summary$mean_diffML,     difference_summary$group, "T")
  C_ML     <- getTC(difference_summary$mean_diffML,     difference_summary$group, "C")
  T_SUM    <- getTC(difference_summary$mean_diffSUM,    difference_summary$group, "T")
  C_SUM    <- getTC(difference_summary$mean_diffSUM,    difference_summary$group, "C")
  T_MPV    <- getTC(difference_summary$mean_diffMIRTpv, difference_summary$group, "T")
  C_MPV    <- getTC(difference_summary$mean_diffMIRTpv, difference_summary$group, "C")
  T_UPV    <- getTC(difference_summary$mean_diffUNIpv,  difference_summary$group, "T")
  C_UPV    <- getTC(difference_summary$mean_diffUNIpv,  difference_summary$group, "C")
  
  MIRT_EAP_DiD    <- T_EAP - C_EAP
  UNI_EAP_DiD     <- T_EAPu - C_EAPu
  UNI_ML_DiD      <- T_ML - C_ML
  SUM_SCORE_DiD   <- T_SUM - C_SUM
  MIRT_PVmean_DiD <- T_MPV - C_MPV
  UNI_PVmean_DiD  <- T_UPV - C_UPV
  
  plot_data <- data.frame(
    variable = c("SUM_SCORE_DiD", "UNI_ML_DiD", "UNI_EAP_DiD", "MIRT_EAP_DiD",
                 "MIRT_PVmean_DiD", "UNI_PVmean_DiD"),
    value = c(SUM_SCORE_DiD, UNI_ML_DiD, UNI_EAP_DiD, MIRT_EAP_DiD,
              MIRT_PVmean_DiD, UNI_PVmean_DiD)
  )
  plot_data$variable <- factor(
    plot_data$variable,
    levels = c("SUM_SCORE_DiD","UNI_ML_DiD","UNI_EAP_DiD","MIRT_EAP_DiD",
               "MIRT_PVmean_DiD","UNI_PVmean_DiD")
  )
  
  p2 <- ggplot(plot_data, aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(value, 2)), vjust = -0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Difference in Differences Between Groups (T - C)",
         x = "Difference Type",
         y = "DiD (Treat - Control)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p2); plot.2 <<- p2
  
} else {
  
  print("No Plots Requested")
  
}  # END PLOT LOOP


#####################################################################
#ESTS
#####################################################################

if(ests=="True") {
  
  
#--------------------------------------
#MIRT
#--------------------------------------
  
  dat.sub <- dat.s %>%
    select(group, EAP1, EAP2)
  
  df_long <- dat.sub %>%
    pivot_longer(cols = starts_with("EAP"), 
                 names_to = "EAP", 
                 values_to = "Score") %>%
    # Create the new variables
    mutate(
      Trt = if_else(group == "T", 1, 0),
      Post = if_else(EAP == "EAP2", 1, 0)
    )
  
df_long$TP<-df_long$Trt*df_long$Post
  
  
  RG<-lm(formula = Score ~ Trt + Post + TP, data = df_long)
  RG2<-summary(RG)
  RG2$coefficients
  
  # newly adjusted SEs
  RG3<-coeftest(RG, vcov = vcovHC(RG, type = "HC0", cluster = "ID"))
  
  #est
  est<-RG3[4,1]
  #est<-RG2$coefficients[4,1]
  est
  
  #SE
  SE<-RG3[4,2]
  #SE<-RG2$coefficients[4,2]
  SE
  
  #PV
  PVal<-RG3[4,4]
  #PV<-RG2$coefficients[4,4]
  PVal
  
  res.MIRT<-cbind("MIRT",est,SE,PVal)
  
  
  #--------------------------------------
  #SUM
  #--------------------------------------
  
  dat.sub <- dat.s %>%
    select(group, Sum1stand, Sum2stand)
  
  df_long <- dat.sub %>%
    pivot_longer(cols = starts_with("Sum"), 
                 names_to = "Sum", 
                 values_to = "Score") %>%
    # Create the new variables
    mutate(
      Trt = if_else(group == "T", 1, 0),
      Post = if_else(Sum == "Sum2stand", 1, 0)
    )
  
  df_long$TP<-df_long$Trt*df_long$Post
  
  
  RG<-lm(formula = Score ~ Trt + Post + TP, data = df_long)
  RG2<-summary(RG)
  RG2$coefficients
  
  # newly adjusted SEs
  RG3<-coeftest(RG, vcov = vcovHC(RG, type = "HC0", cluster = "ID"))
  
  #est
  est<-RG3[4,1]
  #est<-RG2$coefficients[4,1]
  est
  
  #SE
  SE<-RG3[4,2]
  #SE<-RG2$coefficients[4,2]
  SE
  
  #PV
  PVal<-RG3[4,4]
  #PV<-RG2$coefficients[4,4]
  PVal  
  
  res.SUM<-cbind("Sum",est,SE,PVal)
  
  #--------------------------------------
  #EAP UNI
  #--------------------------------------
  
  
  dat.sub <- dat.s %>%
    select(group, EAPu.T1, EAPu.T2)
  
  df_long <- dat.sub %>%
    pivot_longer(cols = starts_with("EAP"), 
                 names_to = "EAP", 
                 values_to = "Score") %>%
    # Create the new variables
    mutate(
      Trt = if_else(group == "T", 1, 0),
      Post = if_else(EAP == "EAPu.T2", 1, 0)
    )
  
  df_long$TP<-df_long$Trt*df_long$Post
  
  
  RG<-lm(formula = Score ~ Trt + Post + TP, data = df_long)
  RG2<-summary(RG)
  RG2$coefficients
  
  # newly adjusted SEs
  RG3<-coeftest(RG, vcov = vcovHC(RG, type = "HC0", cluster = "ID"))
  
  #est
  est<-RG3[4,1]
  #est<-RG2$coefficients[4,1]
  est
  
  #SE
  SE<-RG3[4,2]
  #SE<-RG2$coefficients[4,2]
  SE
  
  #PV
  PVal<-RG3[4,4]
  #PV<-RG2$coefficients[4,4]
  PVal
  
  res.UNIeap<-cbind("UNIDIM. EAP",est,SE,PVal)
  
  #--------------------------------------
  #MLE UNI
  #--------------------------------------
  
  
  dat.sub <- dat.s %>%
    select(group, MLu.T1, MLu.T2)
  
  df_long <- dat.sub %>%
    pivot_longer(cols = starts_with("ML"), 
                 names_to = "ML", 
                 values_to = "Score") %>%
    # Create the new variables
    mutate(
      Trt = if_else(group == "T", 1, 0),
      Post = if_else(ML == "MLu.T2", 1, 0)
    )
  
  df_long$TP<-df_long$Trt*df_long$Post
  
  
  RG<-lm(formula = Score ~ Trt + Post + TP, data = df_long)
  RG2<-summary(RG)
  RG2$coefficients
  
  # newly adjusted SEs
  RG3<-coeftest(RG, vcov = vcovHC(RG, type = "HC0", cluster = "ID"))
  
  #est
  est<-RG3[4,1]
  #est<-RG2$coefficients[4,1]
  est
  
  #SE
  SE<-RG3[4,2]
  #SE<-RG2$coefficients[4,2]
  SE
  
  #PV
  PVal<-RG3[4,4]
  #PV<-RG2$coefficients[4,4]
  PVal
  
  res.UNImle<-cbind("UNIDIM. MLE",est,SE,PVal)
  
  
  
  #--------------------------------------
  #MIRT.PV
  #--------------------------------------
  
  # 0) Add a unique person ID to your base data
  dat.s <- dat.s %>%
    mutate(ID = dplyr::row_number())
  
  # 1) Keep group, ID, and plausible values; reshape to long
  dat.pv <- dat.s %>%
    select(ID, group, starts_with("T1.PV"), starts_with("T2.PV"))
  
  df_long <- dat.pv %>%
    pivot_longer(
      cols = matches("^T[12]\\.PV\\d+$"),
      names_to = c("Time", "PV"),
      names_pattern = "^T(\\d)\\.PV(\\d+)$",
      values_to = "Score"
    ) %>%
    mutate(
      Trt  = if_else(group == "T", 1, 0),
      Post = if_else(Time == "2", 1, 0),
      TP   = Trt * Post,
      imp  = as.integer(PV)  # 1..10
    )
  
  # 2) Fit DiD model within each PV with cluster-robust SEs by ID
  per_imp <- lapply(split(df_long, df_long$imp), function(d) {
    mod <- lm(Score ~ Trt + Post + TP, data = d, na.action = na.omit)
    # Clustered (by ID) HC0 covariance
    Vcl <- vcovCL(mod, cluster = d$ID, type = "HC0")
    rg  <- coeftest(mod, vcov. = Vcl)
    c(est = rg["TP", 1], se = rg["TP", 2])
  })
  
  per_imp_df <- do.call(rbind, per_imp)
  
  # 3) Pool across plausible values (Rubin's rules)
  mi_res <- MIcombine(
    results   = as.list(per_imp_df[, "est"]),
    variances = as.list(per_imp_df[, "se"]^2)
  )
  mi_sum <- summary(mi_res)
  
  mis <- summary(mi_res)
  
  # mi: result of MIcombine(...)
  # mis: summary(mi)
  
  # 1) Pull pooled estimate and SE as scalars
  est <- as.numeric(drop(mis$results[1]))
  SE  <- as.numeric(drop(mis$se[1]))
  
  # 2) Compute t-stat and p-value
  tval <- est / SE
  
  # Try to use MI degrees of freedom if available; otherwise normal approx
  df <- tryCatch(as.numeric(mi$df)[1], error = function(e) NA_real_)
  if (is.na(df) || !is.finite(df)) {
    PVal <- 2 * pnorm(abs(tval), lower.tail = FALSE)  # z-approx
  } else {
    PVal <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)  # t with MI df
  }
  
  # 3) Return a clean one-row data frame
  res.pv <- data.frame(
    V1 = "MIRT PV",
    est = est,
    SE  = SE,
    PVal = PVal,
    stringsAsFactors = FALSE
  )
  res.pv
  
  
  
  
  #--------------------------------------
  # UNI.PV (with T1.PVu*, T2.PVu* variables)
  #--------------------------------------
  

  
  # 0) Add a unique person ID to your base data
  dat.s <- dat.s %>%
    mutate(ID = dplyr::row_number())
  
  # 1) Keep group, ID, and plausible values; reshape to long
  dat.pv <- dat.s %>%
    select(ID, group, starts_with("T1.PVu"), starts_with("T2.PVu"))
  
  df_long <- dat.pv %>%
    pivot_longer(
      cols = matches("^T[12]\\.PVu\\d+$"),
      names_to = c("Time", "PV"),
      names_pattern = "^T(\\d)\\.PVu(\\d+)$",
      values_to = "Score"
    ) %>%
    mutate(
      Trt  = if_else(group == "T", 1, 0),
      Post = if_else(Time == "2", 1, 0),
      TP   = Trt * Post,
      imp  = as.integer(PV)  # 1..m (e.g. 10 PVs)
    )
  
  # 2) Fit DiD model within each PV with cluster-robust SEs by ID
  per_imp <- lapply(split(df_long, df_long$imp), function(d) {
    mod <- lm(Score ~ Trt + Post + TP, data = d, na.action = na.omit)
    Vcl <- vcovCL(mod, cluster = d$ID, type = "HC0")  # cluster-robust
    rg  <- coeftest(mod, vcov. = Vcl)
    c(est = rg["TP", 1], se = rg["TP", 2])
  })
  
  per_imp_df <- do.call(rbind, per_imp)
  
  # 3) Pool across plausible values (Rubin's rules)
  mi_res <- MIcombine(
    results   = as.list(per_imp_df[, "est"]),
    variances = as.list(per_imp_df[, "se"]^2)
  )
  mis <- summary(mi_res)
  
  # 4) Pull pooled estimate and SE as scalars
  est <- as.numeric(drop(mis$results[1]))
  SE  <- as.numeric(drop(mis$se[1]))
  
  # 5) Compute t-stat and p-value
  tval <- est / SE
  df   <- tryCatch(as.numeric(mi$df)[1], error = function(e) NA_real_)
  if (is.na(df) || !is.finite(df)) {
    PVal <- 2 * pnorm(abs(tval), lower.tail = FALSE)  # z-approx
  } else {
    PVal <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)  # t with MI df
  }
  
  # 6) Return a clean one-row data frame
  res.pv.uni <- data.frame(
    V1 = "Unidim. PV",
    est = est,
    SE  = SE,
    PVal = PVal,
    stringsAsFactors = FALSE
  )
  res.pv.uni
  
  
  
  
  
  #--------------------------------------
  #COMBINE
  #--------------------------------------
  
ests.all<-rbind(res.SUM,res.UNImle,res.UNIeap,res.MIRT,res.pv,res.pv.uni)  
  print(ests.all)
  
  ests.all<<-ests.all
  
  
  #--------------------------------------
  #PLOT ESTS
  #--------------------------------------
  
  
  ests_plotdat <- ests.all %>%
    rename(score_approach = V1) %>%
    mutate(
      est = as.numeric(est),
      SE  = as.numeric(SE),
      ci  = 1.96 * SE,
      lower = est - ci,
      upper = est + ci,
      score_approach = reorder(score_approach, est)
    )
  
  p3 <- ggplot(ests_plotdat, aes(x = score_approach, y = est)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
    geom_point(size = 2) +
    labs(
      title   = "Estimates with 95% Intervals (±1.96·SE)",
      x       = "Score approach",
      y       = "Estimate",
      caption = "Each point shows an estimate for a scoring approach with ±1.96×SE error bars (approx. 95% CI). 
The dashed line marks zero, so you can quickly see which estimates are statistically different from zero."
    ) +
    theme_minimal(base_size = 12)
  

  
  
  
  print(p3)
  plot.3<<-p3
  
  
} else {
  
print("No Estimates Requested")  

} #END EST LOOPS

  
  
}




results.RCT(dat.scores)

results.RCT(dat.scores,"False","True")



