
#https://chatgpt.com/c/68cd6bfa-3510-8321-a1aa-11f42e67b729

results.GROWTH <- function(
    score.file,
    plots = TRUE,
    ests  = TRUE,
    timepoints,
    quad  = FALSE
) {
  # ---- dependencies ----
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Please install 'stringr'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  
  # ---- coerce flags (accept TRUE/FALSE or "True"/"False") ----
  as_logical <- function(x) {
    if (is.logical(x)) return(x)
    if (is.character(x)) return(tolower(x) %in% c("true", "t", "1", "yes", "y"))
    as.logical(x)
  }
  plots <- as_logical(plots)
  ests  <- as_logical(ests)
  quad  <- as_logical(quad)
  
  # ---- read data or accept data frame ----
  dat <- if (is.data.frame(score.file)) {
    score.file
  } else if (is.character(score.file) && length(score.file) == 1 && file.exists(score.file)) {
    read.csv(score.file, check.names = FALSE)
  } else {
    stop("`score.file` must be a data.frame or a valid path to a .csv file.")
  }
  
  # ---- helper: safely compute rowMeans over a regex selection ----
  row_mean_by_regex <- function(data, pattern) {
    cols <- grep(pattern, names(data), value = TRUE)
    if (length(cols) == 0) return(rep(NA_real_, nrow(data)))
    rowMeans(data[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  # ---- build per-timepoint person-level PV means (MIRT PV and unidim PVu) ----
  # creates new columns: T1.PVmean, ..., Tt.PVmean, and T1.PVu_mean, ..., Tt.PVu_mean
  for (t in seq_len(timepoints)) {
    # MIRT PV columns like T1.pv1, T1.pv2, ...
    pv_mean_vec  <- row_mean_by_regex(dat, paste0("^T", t, "\\.pv\\d+$"))
    if (!all(is.na(pv_mean_vec))) {
      dat[[paste0("T", t, ".PVmean")]] <- pv_mean_vec
    }
    
    # Unidim PV columns like T1.PVu1, T1.PVu2, ...
    pvu_mean_vec <- row_mean_by_regex(dat, paste0("^T", t, "\\.PVu\\d+$"))
    if (!all(is.na(pvu_mean_vec))) {
      dat[[paste0("T", t, ".PVu_mean")]] <- pvu_mean_vec
    }
  }
  
  # ---- define score families & their timepoint-specific name templates ----
  # Each entry: label = user-facing name; template = function(t) -> column name
  families <- list(
    EAP      = function(t) paste0("EAP.T", t),
    EAPu     = function(t) paste0("EAPu.T", t),
    MLu      = function(t) paste0("MLu.T", t),
    Sum_std  = function(t) paste0("Sum", t, "stand"),
    PVmean   = function(t) paste0("T", t, ".PVmean"),
    PVu_mean = function(t) paste0("T", t, ".PVu_mean")
  )
  
  # ---- gather means by timepoint & score family (column means over people) ----
  summary_list <- list()
  for (fam in names(families)) {
    for (t in seq_len(timepoints)) {
      colname <- families[[fam]](t)
      if (colname %in% names(dat)) {
        mu <- mean(dat[[colname]], na.rm = TRUE)
        n  <- sum(!is.na(dat[[colname]]))
        summary_list[[length(summary_list) + 1]] <- data.frame(
          Timepoint  = t,
          Score_Type = fam,
          Mean       = mu,
          N          = n,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(summary_list) == 0) {
    stop("No matching score columns found for the specified `timepoints`.")
  }
  
  summary_df <- do.call(rbind, summary_list) %>%
    arrange(Score_Type, Timepoint)
  
  
  # ---- relabel Score_Type for legend ----
  summary_df <- summary_df %>%
    mutate(Score_Type = dplyr::recode(
      Score_Type,
      "EAPu"     = "EAP.unidim",
      "MLu"      = "ML.unidim",
      "PVu_mean" = "PV.unidim",
      "PVmean"   = "PV.mirt",
      "EAP"      = "EAP.mirt",
      "Sum_std"  = "Sum_stand"
    ))
  
  # ---- optional plot ----
  p <- NULL
  if (plots) {
    p <- ggplot(summary_df, aes(x = Timepoint, y = Mean, group = Score_Type, linetype = Score_Type)) +
      geom_line() +
      geom_point(aes(shape = Score_Type)) +
      labs(
        x = "Timepoint",
        y = "Mean score",
        linetype = "Score type",
        shape = "Score type",
        title = "Growth of Mean Scores Across Timepoints"
      ) +
      scale_x_continuous(breaks = seq_len(timepoints)) +
      theme_minimal(base_size = 12)
    
    if (quad) {
      # add quadratic trends per series
      p <- p + geom_smooth(
        method = "lm",
        formula = y ~ poly(x, 2),
        se = FALSE,
        aes(group = Score_Type)
      )
    }
    
    print(p)
    plot.1<<-p
  }
  
  # ---- return ----
  out <- list()
  if (ests) out$estimates <- summary_df
  if (plots) out$plot <- p
  invisible(out)
  
  
  # ---- latent growth curve models (lavaan), incl. MI pooling for PVs ----
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("Please install 'lavaan'.")
  
  
    build_lgcm_syntax <- function(varnames, quadratic = FALSE) {
      # time index 0,1,2,... for linear slope loadings
      t_idx <- seq_along(varnames) - 1L
      
      syn <- paste0(
        # loadings
        "i =~ ", paste0("1*", varnames, collapse = " + "), "\n",
        "s =~ ", paste0(t_idx, "*", varnames, collapse = " + "), "\n",
        # freely estimate latent means (otherwise lavaan fixes them to 0)
        "i ~ 1\n",
        "s ~ 1\n"
      )
      
      if (quadratic) {
        syn <- paste0(
          syn,
          "q =~ ", paste0((t_idx^2), "*", varnames, collapse = " + "), "\n",
          "q ~ 1\n"
        )
      }
      syn
    }
    
  
  extract_core_params <- function(fit, quadratic = FALSE) {
    pe <- lavaan::parameterEstimates(fit, standardized = FALSE)
    need_means <- c("i", "s", if (quadratic) "q" else NULL)
    need_vars  <- c("i", "s", if (quadratic) "q" else NULL)
    
    means_df <- pe[pe$op == "~1" & pe$lhs %in% need_means, c("lhs", "est", "se")]
    vars_df  <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% need_vars, c("lhs", "est", "se")]
    
    means_df$param <- paste0(means_df$lhs, "_mean")
    vars_df$param  <- paste0(vars_df$lhs, "_var")
    
    out <- rbind(
      data.frame(param = means_df$param, est = means_df$est, se = means_df$se),
      data.frame(param = vars_df$param,  est = vars_df$est,  se = vars_df$se)
    )
    out
  }
  
  rubin_pool <- function(est_vec, se_vec) {
    m <- length(est_vec)
    if (m == 1L) return(c(est = est_vec, se = se_vec))
    W <- mean(se_vec^2, na.rm = TRUE)
    B <- stats::var(est_vec, na.rm = TRUE)
    Tvar <- W + (1 + 1/m) * B
    c(est = mean(est_vec, na.rm = TRUE), se = sqrt(Tvar))
  }
  
  fit_one_family <- function(df, varnames, quadratic, model_label, score_label) {
    syn <- build_lgcm_syntax(varnames, quadratic = quadratic)
    fit <- lavaan::lavaan(
      model = syn, data = df, estimator = "MLR", missing = "fiml",
      meanstructure = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE, auto.fix.first = FALSE,
      warn = FALSE
    )
    pars <- extract_core_params(fit, quadratic = quadratic)
    pars$model      <- model_label
    pars$score_type <- score_label
    pars
  }
  
  # helper to build data.frames for PV imputations k=1..K
  build_pv_imputations <- function(dat, prefix_time = "T", pv_prefix, timepoints) {
    # pv_prefix examples: "pv" for MIRT PVs => columns like T1.pv1, T2.pv1, ...
    #                     "PVu" for unidim PVs => T1.PVu1, T2.PVu1, ...
    # find how many imputations K by scanning T1
    k_cols <- grep(paste0("^", prefix_time, "1\\.", pv_prefix, "\\d+$"), names(dat), value = TRUE)
    if (!length(k_cols)) return(list())
    K <- max(as.integer(gsub(paste0("^", prefix_time, "1\\.", pv_prefix, "(\\d+)$"), "\\1", k_cols)))
    
    imps <- vector("list", K)
    for (k in seq_len(K)) {
      vn <- paste0(prefix_time, seq_len(timepoints), ".", pv_prefix, k)
      keep <- vn[vn %in% names(dat)]
      if (length(keep) != timepoints) next
      dfk <- dat[, keep, drop = FALSE]
      names(dfk) <- paste0("y", seq_len(timepoints))
      imps[[k]] <- dfk
    }
    # remove any NULLs if missing timepoints broke a slot
    Filter(Negate(is.null), imps)
  }
  
  lgcm_rows <- list()
  
  # 1) Non-PV single-dataset families
  single_fams <- list(
    "EAP.mirt"   = function(t) paste0("EAP.T", t),
    "EAP.unidim" = function(t) paste0("EAPu.T", t),
    "ML.unidim"  = function(t) paste0("MLu.T", t),
    "Sum_stand"  = function(t) paste0("Sum", t, "stand")
  )
  
  for (lab in names(single_fams)) {
    vn <- vapply(seq_len(timepoints), single_fams[[lab]], character(1))
    present <- vn %in% names(dat)
    if (!all(present)) next
    df <- dat[, vn, drop = FALSE]
    names(df) <- paste0("y", seq_len(timepoints))
    # linear
    lgcm_rows[[length(lgcm_rows) + 1]] <- fit_one_family(df, names(df), FALSE, "linear", lab)
    # quadratic
    if (timepoints >= 3) {
      lgcm_rows[[length(lgcm_rows) + 1]] <- fit_one_family(df, names(df), TRUE,  "quadratic", lab)
    }
  }
  
  # 2) PV families with MI pooling
  #    (a) MIRT PVs: Tt.pv1..K
  pv_imps_mirt <- build_pv_imputations(dat, pv_prefix = "pv",  timepoints = timepoints)
  if (length(pv_imps_mirt)) {
    # linear
    ests <- lapply(pv_imps_mirt, function(df) fit_one_family(df, names(df), FALSE, "linear", "PV.mirt"))
    # pool by param
    pool_df <- do.call(rbind, ests)
    pooled <- do.call(rbind, lapply(split(pool_df, pool_df$param), function(dd) {
      z <- rubin_pool(dd$est, dd$se)
      data.frame(param = unique(dd$param), est = z["est"], se = z["se"], model = "linear", score_type = "PV.mirt")
    }))
    lgcm_rows[[length(lgcm_rows) + 1]] <- pooled
    
    if (timepoints >= 3) {
      ests_q <- lapply(pv_imps_mirt, function(df) fit_one_family(df, names(df), TRUE, "quadratic", "PV.mirt"))
      pool_dfq <- do.call(rbind, ests_q)
      pooled_q <- do.call(rbind, lapply(split(pool_dfq, pool_dfq$param), function(dd) {
        z <- rubin_pool(dd$est, dd$se)
        data.frame(param = unique(dd$param), est = z["est"], se = z["se"], model = "quadratic", score_type = "PV.mirt")
      }))
      lgcm_rows[[length(lgcm_rows) + 1]] <- pooled_q
    }
  }
  
  #    (b) Unidim PVs: Tt.PVu1..K
  pv_imps_uni <- build_pv_imputations(dat, pv_prefix = "PVu", timepoints = timepoints)
  if (length(pv_imps_uni)) {
    ests <- lapply(pv_imps_uni, function(df) fit_one_family(df, names(df), FALSE, "linear", "PV.unidim"))
    pool_df <- do.call(rbind, ests)
    pooled <- do.call(rbind, lapply(split(pool_df, pool_df$param), function(dd) {
      z <- rubin_pool(dd$est, dd$se)
      data.frame(param = unique(dd$param), est = z["est"], se = z["se"], model = "linear", score_type = "PV.unidim")
    }))
    lgcm_rows[[length(lgcm_rows) + 1]] <- pooled
    
    if (timepoints >= 3) {
      ests_q <- lapply(pv_imps_uni, function(df) fit_one_family(df, names(df), TRUE, "quadratic", "PV.unidim"))
      pool_dfq <- do.call(rbind, ests_q)
      pooled_q <- do.call(rbind, lapply(split(pool_dfq, pool_dfq$param), function(dd) {
        z <- rubin_pool(dd$est, dd$se)
        data.frame(param = unique(dd$param), est = z["est"], se = z["se"], model = "quadratic", score_type = "PV.unidim")
      }))
      lgcm_rows[[length(lgcm_rows) + 1]] <- pooled_q
    }
  }
  
  lgcm_estimates <- if (length(lgcm_rows)) {
    do.call(rbind, lgcm_rows) %>%
      dplyr::select(score_type, model, param, est, se) %>%
      dplyr::arrange(score_type, model, param)
  } else {
    warning("No LGCM estimates were produced; check available columns/timepoints.")
    NULL
  }
  
  # save into function return
  out$lgcm_estimates <- lgcm_estimates
  
  lgcm_results<<-lgcm_estimates
  
  
  # ---- Model-based trends for ALL score types (linear & quadratic),
  #      LRT computed ONLY for EAP.mirt (quadratic vs linear) ----
  
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("Please install 'lavaan'.")
  
  get_latent_means <- function(fit, quadratic = FALSE) {
    pe <- lavaan::parameterEstimates(fit)
    i  <- pe$est[pe$op == "~1" & pe$lhs == "i"]
    s  <- pe$est[pe$op == "~1" & pe$lhs == "s"]
    q  <- if (quadratic) pe$est[pe$op == "~1" & pe$lhs == "q"] else 0
    c(i = i, s = s, q = q)
  }
  
  # families that exist as single datasets (no MI pooling needed)
  single_fams <- list(
    "EAP.mirt"   = function(t) paste0("EAP.T", t),
    "EAP.unidim" = function(t) paste0("EAPu.T", t),
    "ML.unidim"  = function(t) paste0("MLu.T", t),
    "Sum_stand"  = function(t) paste0("Sum", t, "stand")
  )
  
  pred_rows <- list()
  t_idx <- 0:(timepoints - 1L)
  
  # Keep handles to EAP.mirt fits for LRT
  eap_lin_fit  <- NULL
  eap_quad_fit <- NULL
  
  for (lab in names(single_fams)) {
    vn <- vapply(seq_len(timepoints), single_fams[[lab]], character(1))
    if (!all(vn %in% names(dat))) next
    df <- dat[, vn, drop = FALSE]
    names(df) <- paste0("y", seq_len(timepoints))
    
    # linear
    syn_lin <- build_lgcm_syntax(names(df), quadratic = FALSE)
    fit_lin <- lavaan::lavaan(
      model = syn_lin, data = df, estimator = "MLR", missing = "fiml",
      meanstructure = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE, auto.fix.first = FALSE,
      warn = FALSE
    )
    mu <- get_latent_means(fit_lin, quadratic = FALSE)
    pred_rows[[length(pred_rows) + 1]] <-
      data.frame(Timepoint = seq_len(timepoints),
                 yhat = as.numeric(mu["i"] + mu["s"] * t_idx),
                 Score_Type = lab, Model = "Linear")
    
    if (lab == "EAP.mirt") eap_lin_fit <- fit_lin
    
    # quadratic
    if (timepoints >= 3) {
      syn_q <- build_lgcm_syntax(names(df), quadratic = TRUE)
      fit_q <- lavaan::lavaan(
        model = syn_q, data = df, estimator = "MLR", missing = "fiml",
        meanstructure = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE, auto.fix.first = FALSE,
        warn = FALSE
      )
      muq <- get_latent_means(fit_q, quadratic = TRUE)
      pred_rows[[length(pred_rows) + 1]] <-
        data.frame(Timepoint = seq_len(timepoints),
                   yhat = as.numeric(muq["i"] + muq["s"] * t_idx + muq["q"] * (t_idx^2)),
                   Score_Type = lab, Model = "Quadratic")
      
      if (lab == "EAP.mirt") eap_quad_fit <- fit_q
    }
  }
  
  # PV families with MI pooling → use pooled latent means for predictions
  mk_pv_preds <- function(pv_prefix, label_out) {
    imps <- build_pv_imputations(dat, pv_prefix = pv_prefix, timepoints = timepoints)
    if (!length(imps)) return(NULL)
    
    # linear
    est_lin <- lapply(imps, function(dfk) {
      syn <- build_lgcm_syntax(names(dfk), quadratic = FALSE)
      fit <- lavaan::lavaan(
        model = syn, data = dfk, estimator = "MLR", missing = "fiml",
        meanstructure = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE, auto.fix.first = FALSE,
        warn = FALSE
      )
      pe <- lavaan::parameterEstimates(fit)
      c(i = pe$est[pe$op == "~1" & pe$lhs == "i"],
        s = pe$est[pe$op == "~1" & pe$lhs == "s"])
    })
    est_mat <- do.call(rbind, est_lin)
    i_lin <- mean(est_mat[, "i"], na.rm = TRUE)
    s_lin <- mean(est_mat[, "s"], na.rm = TRUE)
    preds <- list(
      data.frame(Timepoint = seq_len(timepoints),
                 yhat = as.numeric(i_lin + s_lin * t_idx),
                 Score_Type = label_out, Model = "Linear")
    )
    
    # quadratic
    if (timepoints >= 3) {
      est_q <- lapply(imps, function(dfk) {
        syn <- build_lgcm_syntax(names(dfk), quadratic = TRUE)
        fit <- lavaan::lavaan(
          model = syn, data = dfk, estimator = "MLR", missing = "fiml",
          meanstructure = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE, auto.fix.first = FALSE,
          warn = FALSE
        )
        pe <- lavaan::parameterEstimates(fit)
        c(i = pe$est[pe$op == "~1" & pe$lhs == "i"],
          s = pe$est[pe$op == "~1" & pe$lhs == "s"],
          q = pe$est[pe$op == "~1" & pe$lhs == "q"])
      })
      est_mat_q <- do.call(rbind, est_q)
      i_q <- mean(est_mat_q[, "i"], na.rm = TRUE)
      s_q <- mean(est_mat_q[, "s"], na.rm = TRUE)
      q_q <- mean(est_mat_q[, "q"], na.rm = TRUE)
      preds[[length(preds) + 1]] <-
        data.frame(Timepoint = seq_len(timepoints),
                   yhat = as.numeric(i_q + s_q * t_idx + q_q * (t_idx^2)),
                   Score_Type = label_out, Model = "Quadratic")
    }
    do.call(rbind, preds)
  }
  
  # PV.mirt (Tt.pv#) and PV.unidim (Tt.PVu#)
  pv_mirt_preds <- mk_pv_preds("pv",  "PV.mirt")
  pv_uni_preds  <- mk_pv_preds("PVu", "PV.unidim")
  if (!is.null(pv_mirt_preds)) pred_rows[[length(pred_rows) + 1]] <- pv_mirt_preds
  if (!is.null(pv_uni_preds))  pred_rows[[length(pred_rows) + 1]] <- pv_uni_preds
  
  pred_df_all <- if (length(pred_rows)) do.call(rbind, pred_rows) else {
    stop("No model-based predictions could be computed; check columns/timepoints.")
  }
  
  # LRT ONLY for EAP.mirt
  lrt_label <- "LRT (Quad vs Linear) — EAP.mirt: not available"
  if (!is.null(eap_lin_fit) && !is.null(eap_quad_fit)) {
    lrt_tab <- as.data.frame(lavaan::lavTestLRT(eap_lin_fit, eap_quad_fit))
    chi_col <- intersect(c("Chisq.diff.scaled", "Chisq diff.scaled", "Chisq.diff", "Chisq diff"), names(lrt_tab))[1]
    p_col   <- intersect(c("Pr(>Chisq).scaled", "Pr(>Chisq)"), names(lrt_tab))[1]
    df_col  <- intersect(c("Df.diff", "Df diff"), names(lrt_tab))[1]
    if (!is.na(chi_col) && !is.na(p_col) && !is.na(df_col)) {
      delta_chi <- lrt_tab[[chi_col]][nrow(lrt_tab)]
      delta_df  <- lrt_tab[[df_col]][nrow(lrt_tab)]
      pval      <- lrt_tab[[p_col]][nrow(lrt_tab)]
      lrt_label <- sprintf("LRT (Quadratic vs Linear) — EAP.mirt: Δχ²(df=%d)=%.2f, p=%.3g",
                           as.integer(delta_df), as.numeric(delta_chi), as.numeric(pval))
    }
  }
  
  # Two plots: one for Linear, one for Quadratic (if available), showing ALL score types
  mk_trend_plot <- function(df_pred, title_txt, caption_txt) {
    ggplot(df_pred, aes(x = Timepoint, y = yhat,
                        group = Score_Type, linetype = Score_Type, shape = Score_Type)) +
      geom_line() +
      geom_point() +
      labs(x = "Timepoint", y = "Model-based mean", title = title_txt, caption = caption_txt) +
      scale_x_continuous(breaks = seq_len(timepoints)) +
      theme_minimal(base_size = 12) +
      theme(plot.caption = element_text(hjust = 0.5))
  }
  
  p_all_linear  <- mk_trend_plot(subset(pred_df_all, Model == "Linear"),
                                 "LGCM Model-Based Trends (Linear)",  lrt_label)
  p_all_quadratic <- NULL
  if (timepoints >= 3 && any(pred_df_all$Model == "Quadratic")) {
    p_all_quadratic <- mk_trend_plot(subset(pred_df_all, Model == "Quadratic"),
                                     "LGCM Model-Based Trends (Quadratic)", lrt_label)
  }
  
  if (plots) {
    print(p_all_linear)
    if (!is.null(p_all_quadratic)) print(p_all_quadratic)
  }
  
  out$trend_plots <- c(out$trend_plots,
                       list(AllScores_linear = p_all_linear,
                            AllScores_quadratic = p_all_quadratic))
  
  
plot.2<<-p_all_linear
plot.3<<-p_all_quadratic
  
}


results.GROWTH(dat.scores,"True","True",4,"False")
