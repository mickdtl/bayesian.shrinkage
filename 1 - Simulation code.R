# ==========================#
# SIMULATION POWER ANALYSIS #
# ==========================#


# --- libraries --- #
library(tidyverse)
library(truncnorm)
library(glmnet)
library(bayesreg)
library(car)
library(BoomSpikeSlab)
library(future)
library(furrr)
library(progressr)
library(pheatmap)
library(patchwork)
library(ggforce)


# =================== #
# PARAMETRI GLOBALI   #
# =================== #

# --- Extended Run  ------------------------------------------------------------
SAMPLE_SIZES <- seq(20, 100, by = 10)   # 9 livelli: 20 30 40 50 60 70 80 90 100
N_SIM        <- 500                      # simulazioni per sample size
CHUNK_SIZE   <- 10                       # checkpoint ogni 10 sim

# --- Workers -----------------------------------------------------------------
N_WORKERS <- min(7, parallelly::availableCores() - 2)
message(sprintf("Worker paralleli: %d", N_WORKERS))

# --- Checkpoint dir ----------------------------------------------------------
CHECKPOINT_DIR <- "sim_checkpoints_v11"
if (!dir.exists(CHECKPOINT_DIR)) dir.create(CHECKPOINT_DIR)


# ========================== #
# ----- FUNZIONE CORE -----  #
# ========================== #

do_sim <- function(n) {
  
  # --- Helpers ---------------------------------------------------------------
  quiet_model <- function(expr) {
    suppressWarnings(suppressMessages(capture.output(res <- eval(expr))))
    return(res)
  }
  
  extract_sign <- function(mod_obj) {
    suppressWarnings(suppressMessages(capture.output(res <- summary(mod_obj))))
    ci <- res$CI.coef
    ci <- ci[rownames(ci) != "(Intercept)", , drop = FALSE]
    as.vector(sign(ci[, 1]) == sign(ci[, 2]))
  }
  
  # make_var: estrae media/sd dal gruppo e genera dati truncated-normal.
  # Fallback su media/sd globali di colonna se il gruppo ha <3 osservazioni.
  make_var <- function(outcome_val, col, lower = 0, upper = Inf) {
    sub  <- database_per_simulazione[[col]][
      database_per_simulazione$outcome_early == outcome_val
    ]
    full <- database_per_simulazione[[col]]
    m <- if (sum(!is.na(sub)) >= 3) mean(sub, na.rm = TRUE) else mean(full, na.rm = TRUE)
    s <- if (sum(!is.na(sub)) >= 3) sd(sub,   na.rm = TRUE) else sd(full,   na.rm = TRUE)
    # sd=0 o NA -> piccola perturbazione per evitare degenerate truncnorm
    if (is.na(s) || s == 0) s <- 0.01
    rtruncnorm(n, a = lower, b = upper, mean = m, sd = s)
  }
  
  # --- Data Generation ------------------------------------------------------
  nprs   <- c(make_var("Non Responder", "baseline_NPRS",       0, 10),
              make_var("Responder",     "baseline_NPRS",       0, 10))
  il9    <- c(make_var("Non Responder", "baseline_apc_il9",    0),
              make_var("Responder",     "baseline_apc_il9",    0))
  mip1b  <- c(make_var("Non Responder", "baseline_apc_mip1b",  0),
              make_var("Responder",     "baseline_apc_mip1b",  0))
  rantes <- c(make_var("Non Responder", "baseline_apc_rantes", 0),
              make_var("Responder",     "baseline_apc_rantes", 0))
  vegf   <- c(make_var("Non Responder", "baseline_apc_vegf",   0),
              make_var("Responder",     "baseline_apc_vegf",   0))
  
  x1 <- rbinom(2 * n, 1, 0.6)
  x2 <- rgamma(2 * n, shape = 2, scale = 1)
  x3 <- rbeta(2 * n, 1, 2)
  x4 <- rnorm(2 * n)
  
  responder <- c(rep(0, n), rep(1, n))
  dat <- data.frame(responder, nprs, il9, mip1b, rantes, vegf, x1, x2, x3, x4)
  
  var_names <- c("nprs", "il9", "mip1b", "rantes", "vegf", "x1", "x2", "x3", "x4")
  true_vars <- c("nprs", "il9", "mip1b", "rantes", "vegf")
  n_vars    <- length(var_names)
  
  # --- Models  ---------------------------------------------------------------
  
  # GLM — tryCatch per separazione completa (frequente a n=20)
  vector_glm <- tryCatch({
    mod        <- glm(responder ~ ., family = "binomial", data = dat)
    mod_sum    <- summary(mod)
    coef_table <- mod_sum$coefficients
    vif_values <- suppressWarnings(vif(mod))
    setNames(
      sapply(var_names, function(var) {
        if (!(var %in% rownames(coef_table))) return(FALSE)
        if (!(var %in% names(vif_values)))    return(FALSE)
        (coef_table[var, "Pr(>|z|)"] < 0.05) & (vif_values[var] < 5)
      }),
      var_names
    )
  }, error = function(e) {
    # In caso di non-convergenza totale: nessuna variabile selezionata
    setNames(rep(FALSE, length(var_names)), var_names)
  })
  
  # LASSO
  mod_lasso <- cv.glmnet(as.matrix(dat[, -1]), dat[, "responder"],
                         family = "binomial", alpha = 1)
  lasso.1se <- as.vector(coef(mod_lasso, s = "lambda.1se"))[-1] != 0
  
  # Bayesian Ridge
  ridge_bayesreg <- quiet_model(quote(
    bayesreg(as.factor(responder) ~ nprs + il9 + mip1b + rantes + vegf +
               x1 + x2 + x3 + x4,
             data = dat, model = "logistic", prior = "ridge")
  ))
  bayes_ridge_sign <- extract_sign(ridge_bayesreg)
  
  # Bayesian Lasso
  lasso_bayesreg <- quiet_model(quote(
    bayesreg(as.factor(responder) ~ nprs + il9 + mip1b + rantes + vegf +
               x1 + x2 + x3 + x4,
             data = dat, model = "logistic", prior = "lasso")
  ))
  bayes_lasso_sign <- extract_sign(lasso_bayesreg)
  
  # Horseshoe
  hs_bayesreg <- quiet_model(quote(
    bayesreg(as.factor(responder) ~ nprs + il9 + mip1b + rantes + vegf +
               x1 + x2 + x3 + x4,
             data = dat, model = "logistic", prior = "horseshoe")
  ))
  horseshoe_sign <- extract_sign(hs_bayesreg)
  
  # Horseshoe+
  hsplus_bayesreg <- quiet_model(quote(
    bayesreg(as.factor(responder) ~ nprs + il9 + mip1b + rantes + vegf +
               x1 + x2 + x3 + x4,
             data = dat, model = "logistic", prior = "horseshoe+")
  ))
  horseshoe_plus_sign <- extract_sign(hsplus_bayesreg)
  
  # Spike-and-Slab (PIP threshold = 0.3)
  mod_spike <- quiet_model(quote(
    logit.spike(responder ~ nprs + il9 + mip1b + rantes + vegf +
                  x1 + x2 + x3 + x4,
                data = dat, niter = 5000)
  ))
  pips_all  <- summary(mod_spike)$coefficients[, "inc.prob"]
  pips      <- pips_all[names(pips_all) != "(Intercept)"]
  spike_sel <- pips > 0.3
  
  # --- Metrics --------------------------------------------------------------
  compute_metrics <- function(selection_vector) {
    selected  <- var_names[as.logical(selection_vector)]
    tp <- sum(selected %in% true_vars)
    fp <- sum(!(selected %in% true_vars))
    fn <- sum(!(true_vars %in% selected))
    tn <- n_vars - tp - fp - fn
    
    precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall    <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    fdr       <- ifelse((tp + fp) == 0, 0, fp / (tp + fp))
    fnr       <- ifelse((tp + fn) == 0, 0, fn / (tp + fn))
    f1        <- ifelse((precision + recall) == 0, 0,
                        2 * precision * recall / (precision + recall))
    mcc_d     <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc       <- ifelse(mcc_d == 0, 0, (tp * tn - fp * fn) / mcc_d)
    
    c(TP = tp, FP = fp, FN = fn, TN = tn,
      FDR = fdr, Precision = precision, Recall = recall,
      F1 = f1, FNR = fnr, MCC = mcc,
      Selected = length(selected),
      ExactRecovery = as.integer(setequal(selected, true_vars)))
  }
  
  metrics_mat <- rbind(
    glm         = compute_metrics(vector_glm),
    lasso_1se   = compute_metrics(lasso.1se),
    bayes_lasso = compute_metrics(bayes_lasso_sign),
    bayes_ridge = compute_metrics(bayes_ridge_sign),
    hs          = compute_metrics(horseshoe_sign),
    hsplus      = compute_metrics(horseshoe_plus_sign),
    spike       = compute_metrics(spike_sel)
  )
  
  selection_mat <- rbind(
    glm         = as.integer(vector_glm),
    lasso_1se   = as.integer(lasso.1se),
    bayes_lasso = as.integer(bayes_lasso_sign),
    bayes_ridge = as.integer(bayes_ridge_sign),
    hs          = as.integer(horseshoe_sign),
    hsplus      = as.integer(horseshoe_plus_sign),
    spike       = as.integer(spike_sel)
  )
  colnames(selection_mat) <- var_names
  
  list(metrics = metrics_mat, selection = selection_mat)
  }


# =============================== #
# ------ CHECKPOINT RUNNER -----  #
# =============================== #

run_with_checkpoints <- function(sample_n, n_sim, chunk_size,
                                 checkpoint_dir, n_workers) {
  chunks   <- split(seq_len(n_sim), ceiling(seq_len(n_sim) / chunk_size))
  n_chunks <- length(chunks)
  all_chunks <- vector("list", n_chunks)
  
  for (k in seq_len(n_chunks)) {
    cp_file <- file.path(checkpoint_dir,
                         sprintf("n%03d_chunk%03d.rds", sample_n, k))
    
    if (file.exists(cp_file)) {
      message(sprintf("  [n=%d | chunk %d/%d] Resume da checkpoint.", sample_n, k, n_chunks))
      all_chunks[[k]] <- readRDS(cp_file)
      next
    }
    
    message(sprintf("  [n=%d | chunk %d/%d] Avvio %d sim con %d worker...",
                    sample_n, k, n_chunks, length(chunks[[k]]), n_workers))
    
    plan(multisession, workers = n_workers)
    handlers(global = TRUE)
    handlers("progress")
    
    chunk_result <- withCallingHandlers(
      with_progress({
        p <- progressr::progressor(steps = length(chunks[[k]]))
        future_map(chunks[[k]], function(i) { p(); do_sim(sample_n) },
                   .options = furrr_options(seed = TRUE))
      }),
      warning = function(w) {
        message("  WARN: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    plan(sequential)
    saveRDS(chunk_result, cp_file)
    message(sprintf("  [n=%d | chunk %d/%d] Salvato: %s", sample_n, k, n_chunks, cp_file))
    all_chunks[[k]] <- chunk_result
    gc()
  }
  
  do.call(c, all_chunks)
}


# ================================ #
# ---------- MAIN LOOP ----------  #
# ================================ #

results_by_n <- list()

for (sn in SAMPLE_SIZES) {
  message(sprintf("\n======= n = %d =======", sn))
  
  sim_results <- run_with_checkpoints(sn, N_SIM, CHUNK_SIZE, CHECKPOINT_DIR, N_WORKERS)
  
  metrics_n <- do.call(rbind, lapply(seq_along(sim_results), function(i) {
    df <- as.data.frame(sim_results[[i]]$metrics)
    df$model <- rownames(df)
    df$sim   <- i
    df$n     <- sn
    df
  })) %>%
    rename(TPR = Recall, Size = Selected) %>%
    mutate(
      Jaccard = ifelse((TP + FP + FN) == 0, 0, TP / (TP + FP + FN)),
      Brier   = (FP + FN) / 9
    )
  
  results_by_n[[as.character(sn)]] <- metrics_n
  message(sprintf("  n=%d: %d righe aggregate.", sn, nrow(metrics_n)))
}


# ========================= #
# ----- AGGREGATION ------  #
# ========================= #

all_metrics <- bind_rows(results_by_n)
saveRDS(all_metrics, file.path(CHECKPOINT_DIR, "all_metrics_v11_500 4di4.rds"))

power_summary <- all_metrics_2000 %>%
  group_by(n, model) %>%
  summarise(
    FDR_mean           = mean(FDR),
    FDR_se             = sd(FDR) / sqrt(dplyr::n()),
    TPR_mean           = mean(TPR),
    TPR_se             = sd(TPR) / sqrt(dplyr::n()),
    MCC_mean           = mean(MCC),
    MCC_se             = sd(MCC) / sqrt(dplyr::n()),
    Jaccard_mean       = mean(Jaccard),
    Jaccard_se         = sd(Jaccard) / sqrt(dplyr::n()),  
    Brier_mean         = mean(Brier),
    Brier_se           = sd(Brier) / sqrt(dplyr::n()),  
    F1_mean            = mean(F1),
    F1_se              = sd(F1) / sqrt(dplyr::n()),
    ExactRecovery_mean = mean(ExactRecovery),
    ExactRecovery_se   = sd(ExactRecovery) / sqrt(dplyr::n()),
    Size               = mean(Size),
    n_sim              = dplyr::n(),
    .groups            = "drop"
  )

# =============================================================================
# PALETTE COLORI — consistente tra tutti i plot
# =============================================================================

model_levels <- sort(unique(all_metrics$model))
model_colors_custom <- c(
  "bayes_lasso" = "#E07B7B",   # rosa
  "bayes_ridge" = "#C8A83C",   # oro
  "glm"         = "#5CB85C",   # verde
  "hs"          = "#4ABFB0",   # teal
  "hsplus"      = "#4A90D9",   # azzurro
  "lasso_1se"   = "#9B6FC8",   # viola
  "spike"       = "#E87CBF"    # rosa acceso
)

# Subset bayesiani + frequentist per highlight nel scatter
bayesian_models    <- c("bayes_lasso", "bayes_ridge", "hs", "hsplus", "spike")
frequentist_models <- c("glm", "lasso_1se")

power_summary <- power_summary %>%
  mutate(
    model_type = ifelse(model %in% bayesian_models, "Bayesian", "Frequentist")
  )


# =============================================================================
# PLOT 1 — Power curves composite (A+B+C+D)
# =============================================================================

ribbon_plot <- function(data, y_mean, y_se, y_lab, title_lab,
                        hline = NULL, hline_lab = NULL) {
  sym_mean <- sym(y_mean)
  sym_se   <- sym(y_se)
  p <- ggplot(data, aes(x = n, y = !!sym_mean,
                        color = model, fill = model, group = model)) +
    geom_ribbon(aes(ymin = !!sym_mean - 1.96 * !!sym_se,
                    ymax = !!sym_mean + 1.96 * !!sym_se),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_color_manual(values = model_colors_custom) +
    scale_fill_manual(values  = model_colors_custom) +
    scale_x_continuous(breaks = SAMPLE_SIZES) +
    labs(title = title_lab, x = "Sample size (n per group)", y = y_lab) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  if (!is.null(hline))
    p <- p +
    geom_hline(yintercept = hline, linetype = "dashed", color = "gray40") +
    annotate("text", x = min(SAMPLE_SIZES), y = hline + 0.02,
             label = hline_lab, hjust = 0, size = 3, color = "gray40")
  p
}

pc_tpr <- ribbon_plot(power_summary, "TPR_mean", "TPR_se",
                      "True Positive Rate (Power)",
                      "A. Power (TPR) vs sample size",
                      hline = 0.80, hline_lab = "0.80")

pc_fdr <- ribbon_plot(power_summary, "FDR_mean", "FDR_se",
                      "False Discovery Rate",
                      "B. FDR vs sample size",
                      hline = 0.05, hline_lab = "FDR = 0.05")

pc_mcc <- ribbon_plot(power_summary, "MCC_mean", "MCC_se",
                      "Matthews Correlation Coefficient",
                      "C. MCC vs sample size",
                      hline = 0, hline_lab = "MCC = 0")

pc_er <- ggplot(power_summary,
                aes(x = n, y = ExactRecovery, color = model, group = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = model_colors_custom) +
  scale_x_continuous(breaks = SAMPLE_SIZES) +
  labs(title = "D. Exact Recovery vs sample size",
       x = "Sample size (n per group)", y = "P(Exact Recovery)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

power_composite <- (pc_tpr | pc_fdr) / (pc_mcc | pc_er)
print(power_composite)


# =============================================================================
# PLOT 2 — FDR-TPR SCATTER (plot indipendente)
#
# Logica:
#   - Un punto per ogni (model, n): posizione = (FDR_mean, TPR_mean)
#   - Dimensione punto = MCC_mean (più grande = meglio globalmente)
#   - Colore = modello
#   - Forma = tipo (Bayesian / Frequentist)
#   - Frecce che collegano lo stesso modello attraverso n crescenti
#     (da n_min a n_max, per vedere la traiettoria nel piano FDR-TPR)
#   - Ellissi di concentrazione per gruppo Bayesian vs Frequentist
#   - Quadranti: riga FDR=0.05, colonna TPR=0.80
#   - Labelling automatico dell'ultimo punto (n più grande) per ogni modello
# =============================================================================

# Ordina per n per poter disegnare le frecce in sequenza
scatter_data <- power_summary %>%
  arrange(model, n)

# Crea dataset segmenti: ogni riga = segmento da (n_k) a (n_{k+1})
arrow_data <- scatter_data %>%
  group_by(model) %>%
  mutate(
    FDR_end = lead(FDR_mean),
    TPR_end = lead(TPR_mean)
  ) %>%
  filter(!is.na(FDR_end)) %>%
  ungroup()

# Labels: solo per l'ultimo n di ogni modello (evita sovraffollamento)
label_data <- scatter_data %>%
  group_by(model) %>%
  slice_max(order_by = n, n = 1) %>%
  ungroup()

fdr_tpr_scatter <- ggplot(scatter_data,
                          aes(x = FDR_mean, y = TPR_mean,
                              color = model, fill = model)) +
  
  # Quadranti di riferimento
  annotate("rect",
           xmin = -Inf, xmax = 0.05, ymin = 0.80, ymax = Inf,
           fill = "green", alpha = 0.04) +
  annotate("text", x = 0.001, y = 0.99,
           label = "Low FDR\nHigh TPR", hjust = 0, vjust = 1,
           size = 3, color = "darkgreen", fontface = "italic") +
  geom_vline(xintercept = 0.05, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  
  # Ellissi per gruppo Bayesian / Frequentist
  #stat_ellipse(aes(group = model_type, color = NULL,
  #                 fill = NULL),
  #             geom = "polygon",
  #             type = "norm", level = 0.80,
  #             alpha = 0.06, color = "gray60",
  #             linetype = "dotted", linewidth = 0.7) +
  
  
  # Punti: dimensione proporzionale a MCC
  geom_point(aes(size = MCC_mean, shape = model_type),
             alpha = 0.85) +
  
  # Frecce traiettoria (n crescente)
  geom_segment(data = arrow_data,
               aes(x = FDR_mean, y = TPR_mean,
                   xend = FDR_end, yend = TPR_end,
                   color = "black"),
               arrow = arrow(length = unit(0.15, "cm"),
                             type = "closed"),
               linewidth = 0.55) +
  
  # Labels modello sull'ultimo punto
  ggrepel::geom_text_repel(data = label_data,
                           aes(label = model),
                           size = 3, fontface = "bold",
                           box.padding = 0.4,
                           show.legend = FALSE) +
  
  scale_color_manual(values = model_colors_custom) +
  scale_fill_manual(values  = model_colors_custom) +
  scale_size_continuous(name   = "MCC (mean)",
                        range  = c(2, 8),
                        breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_shape_manual(name   = "Model type",
                     values = c(Bayesian = 16, Frequentist = 17)) +
  
  labs(
    title    = "FDR–TPR trade-off across sample sizes",
    subtitle = paste0("Each point = one (model, n) combination  |  ",
                      "Point size = MCC  |  Arrows = trajectory as n increases  |  ",
                      "n = ", paste(range(SAMPLE_SIZES), collapse = "–")),
    x        = "False Discovery Rate (FDR)",
    y        = "True Positive Rate (Power)",
    color    = "Model"
  ) +
  
  coord_cartesian(xlim = c(0, max(scatter_data$FDR_mean) * 1.15),
                  ylim = c(min(scatter_data$TPR_mean) * 0.90, 1.02)) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 9, color = "gray40")
  )

print(fdr_tpr_scatter)


# =============================================================================
# PLOT 3 — Violin distributions (MCC, FDR, TPR)
# =============================================================================

make_violin_n <- function(data, metric, title_lab, hline = NULL) {
  p <- ggplot(data,
              aes(x = factor(.data$n), y = .data[[metric]], fill = model)) +
    geom_violin(alpha = 0.5, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, outlier.shape = NA,
                 position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = model_colors_custom) +
    labs(title = title_lab,
         x = "Sample size (n per group)", y = metric) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline,
                        linetype = "dashed", color = "gray40")
  p
}

print(make_violin_n(all_metrics, "MCC",
                    "MCC distribution by n and model", hline = 0))
print(make_violin_n(all_metrics, "FDR",
                    "FDR distribution by n and model", hline = 0.05))
print(make_violin_n(all_metrics, "TPR",
                    "TPR distribution by n and model", hline = 0.80))


# =============================================================================
# PLOT 4 — Performance heatmap per ogni n
# =============================================================================

for (sn in SAMPLE_SIZES) {
  heat_n <- power_summary %>%
    filter(n == sn) %>%
    column_to_rownames("model") %>%
    select(FDR_mean, TPR_mean, MCC_mean, Jaccard_mean,
           Brier_mean, ExactRecovery, Size)
  
  colnames(heat_n) <- c("FDR", "TPR", "MCC", "Jaccard",
                        "Brier", "ExactRecovery", "Size")
  
  pheatmap(scale(as.matrix(heat_n)),
           main          = sprintf("Performance heatmap (n = %d)", sn),
           cluster_rows  = TRUE,
           cluster_cols  = FALSE,
           display_numbers = TRUE,
           number_format = "%.2f",
           color         = colorRampPalette(
             c("steelblue3", "white", "firebrick3"))(50))
}


# =============================================================================
# SUMMARY CONSOLE
# =============================================================================

cat("\n\n========== POWER SUMMARY ==========\n")
print(
  power_summary %>%
    select(n, model, TPR_mean, FDR_mean, MCC_mean, ExactRecovery, F1_mean) %>%
    arrange(n, model),
  n = Inf
)

# Evidenzia n minimo per TPR >= 0.80 per ogni modello
cat("\n\n========== n MINIMO PER TPR >= 0.80 ==========\n")
print(
  power_summary %>%
    filter(TPR_mean >= 0.80) %>%
    group_by(model) %>%
    slice_min(order_by = n, n = 1) %>%
    select(model, n, TPR_mean, FDR_mean, MCC_mean) %>%
    arrange(n)
)
