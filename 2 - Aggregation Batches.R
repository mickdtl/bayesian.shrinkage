library(dplyr)
library(tidyr)
library(gt)

# Carica i due aggregati
batch1 <- readRDS("sim_checkpoints_v11/all_metrics_v11_2000_1di4.rds")
batch2 <- readRDS("sim_checkpoints_v11/all_metrics_v11_2000_2di4.rds")
batch3 <- readRDS("sim_checkpoints_v11/all_metrics_v11_2000_3di4.rds")
batch4 <- readRDS("sim_checkpoints_v11/all_metrics_v11_2000_4di4.rds")

# Ri-indicizza sim per evitare duplicati (batch2 parte da 501)
batch2 <- batch2 %>%
  mutate(sim = sim + max(batch1$sim))
batch3 <- batch3 %>%
  mutate(sim = sim + max(batch2$sim))
batch4 <- batch4 %>%
  mutate(sim = sim + max(batch3$sim))

# Unisci
all_metrics_1000 <- bind_rows(batch1, batch2)
all_metrics_1500 <- bind_rows(batch1, batch2, batch3)
all_metrics_2000 <- bind_rows(batch1, batch2, batch3, batch4)

# Verifica: deve essere 2000 sim × 7 modelli × 9 n = 126000 righe
cat("Righe totali:", nrow(all_metrics_2000), "\n")
cat("Sim per n e modello:\n")
print(all_metrics_2000 %>% 
        group_by(n, model) %>% 
        summarise(n_sim = n(), .groups = "drop") %>% 
        pivot_wider(names_from = model, values_from = n_sim))


# =============================================================================
# ----- AGGREGAZIONE -----
# =============================================================================

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

# estrapolo dati
power_summary %>%
  # Uniamo media e SE in un'unica stringa formattata (es. 0.160 ± 0.005)
  mutate(valore = sprintf("%.3f ± %.3f", ExactRecovery_mean, ExactRecovery_se)) %>%
  # Selezioniamo solo le colonne necessarie alla nuova struttura
  select(n, model, valore) %>% arrange(model) %>%
  # Ruotiamo la tabella portando i modelli come colonne
  pivot_wider(names_from = model, values_from = valore) %>%
  gt()


# =============================================================
# MAPPING LABELS
# =============================================================

model_labels <- c(
  "bayes_lasso" = "Bayesian Lasso",
  "bayes_ridge" = "Bayesian Ridge",
  "glm"         = "GLM",
  "hs"          = "Horseshoe",
  "hsplus"      = "Horseshoe+",
  "lasso_1se"   = "LASSO (1se)",
  "spike"       = "Spike-and-Slab"
)

variable_labels <- c(
  "nprs"   = "NPRS",
  "il9"    = "IL-9",
  "mip1b"  = "MIP-1\u03b2",
  "rantes" = "RANTES",
  "vegf"   = "VEGF",
  "x1"     = "Noise 1",
  "x2"     = "Noise 2",
  "x3"     = "Noise 3",
  "x4"     = "Noise 4"
)

all_metrics_2000$model
all_metrics_2000 <- all_metrics_2000 %>%
  mutate(
    model = case_match(model,
                       "bayes_lasso" ~ "Bayesian Lasso",
                       "bayes_ridge" ~ "Bayesian Ridge",
                       "glm"         ~ "GLM",
                       "hs"          ~ "Horseshoe",
                       "hsplus"      ~ "Horseshoe+",
                       "lasso_1se"   ~ "LASSO (1se)",
                       "spike"       ~ "Spike-and-Slab"
    )
  )
# =============================================================================
# PALETTE COLORI — consistente tra tutti i plot
# =============================================================================

model_levels <- sort(unique(all_metrics_2000$model))
model_colors_custom <- c(
  "bayes_lasso" = "#79AF97FF",
  "bayes_ridge" = "#00A1d5ff",
  "glm"         = "#b24745ff",
  "hs"          = "#374E55FF",
  "hsplus"      = "#80796bff",
  "lasso_1se"   = "#df8f44ff",
  "spike"       = "#6a6599ff"
)

names(model_colors_custom) <- model_labels[names(model_colors_custom)]

# Subset bayesiani + frequentist per highlight nel scatter
bayesian_models    <- c("Bayesian Lasso", "Bayesian Ridge", "Horseshoe", "Horseshoe+", "Spike-and-Slab")
frequentist_models <- c("GLM", "Lasso (1se)")

power_summary <- power_summary %>%
  mutate(
    model_type = ifelse(model %in% bayesian_models, "Bayesian", "Frequentist")
  )


# =============================================================================
# PLOT 1 — Power curves composite (A+B+C+D)
# =============================================================================
library(tidyverse)

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
    labs(title = title_lab, x = "Sample size (n per group)", y = y_lab,
         color = "Model", fill = "Model") +
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
                      "True Positive Rate",
                      "A. TPR vs sample size",
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
                aes(x = n, y = ExactRecovery_mean, color = model, group = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = model_colors_custom) +
  scale_x_continuous(breaks = SAMPLE_SIZES) +
  labs(title = "D. Exact Recovery vs sample size",
       x = "Sample size (n per group)", y = "P(Exact Recovery)", color ="Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

power_composite <- (pc_tpr | pc_fdr) / (pc_mcc | pc_er)
print(power_composite)


# =============================================================================
# PLOT 2 — FDR-TPR SCATTER (plot indipendente)
#
# Logica:
#   - Un punto per ogni (model, n): posizione = (FDR_mean, TPR_mean)
#   - Dimensione punto = sample size
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
scatter_data
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
  
  
  # Punti: dimensione proporzionale a size
  geom_point(aes(size = n, shape = model_type),
             alpha = 0.85) +
  
  # Frecce traiettoria (n crescente)
  geom_segment(data = arrow_data,
               aes(x = FDR_mean, y = TPR_mean,
                   xend = FDR_end, yend = TPR_end,
                   color = "black"),
               arrow = arrow(length = unit(0.15, "cm"),
                             type = "closed"),
               linewidth = 0.55, show.legend = F) +
  
  # Labels modello sull'ultimo punto
  ggrepel::geom_label_repel(data = label_data,
                           aes(label = model,color = model), fill = "white",
                           size = 3, fontface = "bold",
                           box.padding = 0.4,
                           show.legend = FALSE) +
  
  scale_color_manual(values = model_colors_custom, ) +
  scale_fill_manual(values  = model_colors_custom, guide="none") +
  scale_size_continuous(name   = "Sample Size",
                        range  = c(2, 9),
                        breaks = c(20,30,40,50,60,70,80,90,100)) +
  scale_shape_manual(name   = "Model type",
                     values = c(Bayesian = 16, Frequentist = 17)) +
  
  labs(
#    title    = "FDR–TPR trade-off across sample sizes",
#    subtitle = paste0("Each point = one (model, n) combination  |  ",
#                      "Point size = MCC  |  Arrows = trajectory as n increases  |  ",
#                      "n = ", paste(range(SAMPLE_SIZES), collapse = "–")),
    x        = "False Discovery Rate (FDR)",
    y        = "True Positive Rate (TPR)",
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

make_violin_facet <- function(data, metric, title_lab, 
                              hline = NULL, ylim = NULL) {
  p <- ggplot(data,
              aes(x = factor(n), y = .data[[metric]], 
                  fill = model)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA,
                 fill = "white", alpha = 0.7) +
    facet_wrap(~ model, nrow = 2) +
    scale_fill_manual(values = model_colors_custom) +
    labs(title  = title_lab,
         x      = "Sample size (n per group)", 
         y      = metric,
         caption = if (!is.null(ylim))
           "Y-axis truncated for improve readibility" else NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",   # ridondante con facet
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text  = element_text(face = "bold"))
  
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline,
                        linetype = "dashed", color = "gray40")
  if (!is.null(ylim))
    p <- p + coord_cartesian(ylim = ylim)
  p
}

# Chiamate
make_violin_facet(all_metrics_2000, "Jaccard", 
                  "MCC distribution by model and n",
                  hline = 0, ylim = c(-1.0, 1.0))

make_violin_facet(all_metrics_2000, "Brier",
                  "FDR distribution by model and n",
                  hline = 0.05, ylim = c(-1, 1))

make_violin_facet(all_metrics_2000, "F1",
                  "TPR distribution by model and n",
                  hline = 0.80, ylim = c(0, 1.0))

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
  rownames(heat_n) <- c("Bayesian Lasso","Bayesian Ridge","GLM","Horseshoe","Horseshoe+",
                        "LASSO (1se)", "Spike-and-Slab")
  
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

# =============================================================
# SELECTION FREQUENCY DA BATCH 3
# =============================================================

# =============================================================
# DIRECTORY DEI DUE BATCH
# =============================================================

CHECKPOINT_DIRS <- c(
  "sim_checkpoints_v11/500sims - batch 3",
  "sim_checkpoints_v11/500sims - batch 4"
)

SAMPLE_SIZES <- seq(20, 100, by = 10)

# =============================================================
# CARICAMENTO E CALCOLO FREQUENZE (batch 3 + batch 4 uniti)
# =============================================================

sel_freq_by_n <- lapply(SAMPLE_SIZES, function(sn) {
  
  # Raccoglie i chunk da ENTRAMBE le cartelle
  sim_results <- list()
  
  for (dir in CHECKPOINT_DIRS) {
    cp_files <- list.files(dir,
                           pattern = sprintf("n%03d_chunk\\d{3}\\.rds", sn),
                           full.names = TRUE)
    if (length(cp_files) == 0) {
      message(sprintf("  [%s] Nessun chunk per n=%d", basename(dir), sn))
      next
    }
    chunks <- do.call(c, lapply(sort(cp_files), readRDS))
    sim_results <- c(sim_results, chunks)
    message(sprintf("  [%s] n=%d: +%d sim", basename(dir), sn, length(chunks)))
  }
  
  if (length(sim_results) == 0) return(NULL)
  
  # Verifica selection
  if (is.null(sim_results[[1]]$selection)) {
    message(sprintf("n=%d: selection non trovata", sn))
    return(NULL)
  }
  
  B_sim <- length(sim_results) # Numero effettivo di simulazioni (es. 2000)
  message(sprintf("n=%d: TOTALE %d sim", sn, B_sim))
  
  # Array 3D: [Modelli, Variabili, Simulazioni]
  sel_arr <- simplify2array(lapply(sim_results, `[[`, "selection"))
  
  # Calcoliamo media e deviazione standard lungo la terza dimensione (le simulazioni)
  mean_arr <- apply(sel_arr, c(1, 2), mean)
  sd_arr   <- apply(sel_arr, c(1, 2), sd)
  
  # Calcoliamo l'Errore Standard (SE = SD / sqrt(B))
  se_arr   <- sd_arr / sqrt(B_sim)
  
  # Trasformiamo in dataframe e formattiamo in stile "pubblicazione"
  # 1. Gestione Medie
  df_mean <- as.data.frame(mean_arr) %>%
    rownames_to_column("model") %>%
    pivot_longer(-model, names_to = "variable", values_to = "mean_val")
  
  # 2. Gestione SE
  df_se <- as.data.frame(se_arr) %>%
    rownames_to_column("model") %>%
    pivot_longer(-model, names_to = "variable", values_to = "se_val")
  
  # Uniamo i due dataframe e creiamo la stringa "Media ± SE"
  df_mean %>%
    left_join(df_se, by = c("model", "variable")) %>%
    mutate(
      n = sn,
      # Formattazione a 3 decimali (es. 0.456 ± 0.007)
      formatted_val = sprintf("%.3f ± %.3f", mean_val, se_val)
    ) %>%
    select(n, model, variable, formatted_val)
  
}) %>% bind_rows()

# Verifica: atteso lo stesso numero di righe (9 n × 7 modelli × X variabili)
glimpse(sel_freq_by_n)


# estrapolo metriche per tabelle
sel_freq_by_n %>%
  filter(variable == "x4") %>%
  select(-variable) %>% arrange(model) %>%
  pivot_wider(names_from = model, values_from = formatted_val) %>%
  gt

sel_freq_by_n
# =============================================================
# MAPPING LABELS (definisci una volta sola)
# =============================================================

model_labels <- c(
  "bayes_lasso" = "Bayesian Lasso",
  "bayes_ridge" = "Bayesian Ridge",
  "glm"         = "GLM",
  "hs"          = "Horseshoe",
  "hsplus"      = "Horseshoe+",
  "lasso_1se"   = "LASSO (1se)",
  "spike"       = "Spike-and-Slab"
)

variable_labels <- c(
  "nprs"   = "NPRS",
  "il9"    = "IL-9",
  "mip1b"  = "MIP-1\u03b2",
  "rantes" = "RANTES",
  "vegf"   = "VEGF",
  "x1"     = "Noise 1",
  "x2"     = "Noise 2",
  "x3"     = "Noise 3",
  "x4"     = "Noise 4"
)

# Palette con nomi rinominati
model_colors_custom <- c(
  "bayes_lasso" = "#79AF97FF",
  "bayes_ridge" = "#00A1d5ff",
  "glm"         = "#b24745ff",
  "hs"          = "#374E55FF",
  "hsplus"      = "#80796bff",
  "lasso_1se"   = "#df8f44ff",
  "spike"       = "#6a6599ff"
)
names(model_colors_custom) <- model_labels[names(model_colors_custom)]

# Ordini con nomi rinominati
var_order_labels   <- variable_labels[c("nprs","il9","mip1b","rantes","vegf",
                                        "x1","x2","x3","x4")]
model_order_labels <- model_labels[c("lasso_1se","glm","bayes_ridge",
                                     "bayes_lasso","hs","hsplus","spike")]


# =============================================================
# RINOMINA NEL DATAFRAME — applicata UNA volta, usata in entrambi i plot
# =============================================================

sel_freq_renamed <- sel_freq_by_n %>%
  mutate(
    model = case_match(model,
                       "bayes_lasso" ~ "Bayesian Lasso",
                       "bayes_ridge" ~ "Bayesian Ridge",
                       "glm"         ~ "GLM",
                       "hs"          ~ "Horseshoe",
                       "hsplus"      ~ "Horseshoe+",
                       "lasso_1se"   ~ "LASSO (1se)",
                       "spike"       ~ "Spike-and-Slab"
    ),
    variable = case_match(variable,
                          "nprs"   ~ "NPRS",
                          "il9"    ~ "IL-9",
                          "mip1b"  ~ "MIP-1\u03b2",
                          "rantes" ~ "RANTES",
                          "vegf"   ~ "VEGF",
                          "x1"     ~ "Noise 1",
                          "x2"     ~ "Noise 2",
                          "x3"     ~ "Noise 3",
                          "x4"     ~ "Noise 4"
    )
  )


# =============================================================
# PLOT 1 — Heatmap variabile x modello, facet per n
# =============================================================

sel_plot <- sel_freq_renamed %>%
  filter(n %in% c(20, 40, 60, 80, 100)) %>%
  mutate(
    variable = factor(variable, levels = rev(var_order_labels)),
    model    = factor(model,    levels = model_order_labels)
  )

ggplot(sel_plot, aes(x = model, y = variable, fill = sel_freq)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", sel_freq)),
            size = 2.8, color = "black") +
  geom_hline(yintercept = 4.5, color = "black",
             linewidth = 0.8, linetype = "dashed") +
  scale_fill_gradient2(low      = "white",
                       mid      = "skyblue",
                       high     = "steelblue4",
                       midpoint = 0.5,
                       limits   = c(0, 1),
                       name     = "Selection\nFrequency") +
  facet_wrap(~ n, nrow = 1, labeller = label_both) +
  labs(
#    title = "Variable selection frequency by model and sample size",
       x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_text(face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold"),
        panel.grid  = element_blank())


# =============================================================
# PLOT 2 — Power curve per variabile
# =============================================================
library(tidyverse)
models_focus <- c("Bayesian Ridge", "LASSO (1se)", "GLM", "Bayesian Lasso")

sel_curves <- sel_freq_renamed %>%
  mutate(
    var_type = ifelse(variable %in% var_order_labels[1:5],
                      "True signal", "Noise"),
    variable = factor(variable, levels = var_order_labels),
    model    = factor(model,    levels = model_order_labels)
  )
sel_curves$variable
ggplot(sel_curves %>% filter(model %in% models_focus),
       aes(x = n, y = sel_freq,
           color = model, linetype = var_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = model_colors_custom[models_focus],
                     name   = "Model") +
  scale_linetype_manual(values = c("True signal" = "solid",
                                   "Noise"        = "dashed"),
                        name   = "Variable type") +
  scale_x_continuous(breaks = SAMPLE_SIZES) +
  facet_wrap(~ variable, nrow = 2) +
  labs(
#    title = "Per-variable selection frequency across sample sizes",
       x     = "Sample size (n per group)",
       y     = "Selection frequency") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
        strip.text      = element_text(face = "bold"),
        legend.position = "bottom")


library(patchwork)
Fig2a/fig2b+plot_annotation(tag_levels = c("A","B"))
