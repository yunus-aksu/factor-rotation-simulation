# ===============================================================
# --- SIMULATION RESULTS ANALYSİS CODE --
# ===============================================================
Sys.setlocale("LC_ALL", "tr_TR.UTF-8") 

required_packages <- c("nnet", "dplyr", "ggplot2", "data.table",
                       "forcats", "MASS", "scales", "vcd", "rlang", "emmeans") 
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

data_dir <- file.path(getwd(), "simulation_results")
freq_file <- file.path(data_dir, "scenario_frequency.rds")
winner_file <- file.path(data_dir, "scenario_winner.rds")

if (file.exists(freq_file) && file.exists(winner_file)) {
  freq_df <- as.data.frame(readRDS(freq_file))
  winner_df <- as.data.frame(readRDS(winner_file))
  
  freq_df$method <- as.factor(freq_df$method)
  
  if ("none" %in% levels(freq_df$method)) {
    levels(freq_df$method)[levels(freq_df$method) == "none"] <- "unrotated"
    freq_df$method <- relevel(freq_df$method, ref = "unrotated")
    cat("Method faktörü: 'none' -> 'unrotated' olarak güncellendi.\n")
  } else {
    cat("Uyarı: 'none' yöntemi 'method' değişkeninin seviyelerinde bulunamadı. Referans kategori ayarlanmadı.\n")
  }
  
  freq_df$measure <- as.factor(freq_df$measure)
  freq_df$kk_tip <- as.factor(freq_df$kk_tip)
  freq_df$ds <- as.numeric(as.character(freq_df$ds))
  freq_df$np <- as.numeric(as.character(freq_df$np))
  
  winner_df$winner <- as.factor(winner_df$winner)
  
  if ("none" %in% levels(winner_df$winner)) {
    levels(winner_df$winner)[levels(winner_df$winner) == "none"] <- "unrotated"
    winner_df$winner <- relevel(winner_df$winner, ref = "unrotated")
    cat("'winner' faktörü: 'none' -> 'unrotated' olarak güncellendi.\n")
  }
  winner_df$kk_tip <- as.factor(winner_df$kk_tip)
  winner_df$ds <- as.numeric(as.character(winner_df$ds))
  winner_df$np <- as.numeric(as.character(winner_df$np))
  winner_df$winner <- as.factor(winner_df$winner) 
  winner_df$measure <- as.factor(winner_df$measure) 
  
  kk_tip_labels <- c(
    "sabit_kk0p2"           = "Constant correlation (0.2)",
    "sabit_kk0p5"           = "Constant correlation (0.5)",
    "sabit_kk0p8"           = "Constant correlation (0.8)"
  )
  
  freq_df$kk_tip <- factor(freq_df$kk_tip,
                           levels = names(kk_tip_labels),
                           labels = kk_tip_labels)
  
  winner_df$kk_tip <- factor(winner_df$kk_tip,
                             levels = names(kk_tip_labels),
                             labels = kk_tip_labels)
  
} else {
  stop("Hata: Gerekli simülasyon dosyaları bulunamadı. Lütfen dosya yolunu kontrol edin.")
}

gc()

# ================================================================

# ================================================================
measure_label <- function(x) {
  labs <- c(
    "max" = "Maximum loading",
    "min" = "Minimum loading",
    "range_max" = "Maximum loading range"
  )
  labs[as.character(x)]
}
# ================================================================

# ================================================================
my_colors = c( 
  "bentlerT" = "#E41A1C",    
  "entropy" = "#4DAF4A",     
  "geominT" = "#377EB8",     
  "varimax" = "#FFD700",     
  "unrotated" = "#FF7F00",   
  "quartimax" = "#FFFF33",   
  "equamax" = "#A65628",     
  "cf_varimax" = "#F781BF",  
  "infomaxT" = "#999999",    
  "bifactorT" = "#E6AB02",   
  "varimin" = "#00CED1"      
)

# ===============================================================
# --- GRAPHIC FUNCTION ---
# ===============================================================

plot_academic <- function(data,
                          x_var, y_var, group_var = NULL, 
                          facet_var1 = NULL, facet_var2 = NULL, 
                          geom_type = "line_point", 
                          title = "", x_label = "", y_label = "", legend_title = "",
                          max_groups = NULL, 
                          y_is_percent = TRUE,  
                          legend_ncol = NULL,   
                          legend_nrow = NULL,  
                          dpi = 600, 
                          width = 15, height = 10, units = "in", 
                          save_path = "academic_plot.tiff") { 
  
  if (!is.null(max_groups) && !is.null(group_var)) {
    group_var_sym <- rlang::sym(group_var)
    top_groups <- data %>%
      count(!!group_var_sym, sort = TRUE) %>%
      slice_head(n = max_groups) %>%
      pull(!!group_var_sym)
    
    data <- data %>% filter(!!group_var_sym %in% top_groups)
  }
  
  p <- ggplot(data, aes(x = !!rlang::sym(x_var), y = !!rlang::sym(y_var)))
  
  if (geom_type == "line_point") {
    if (!is.null(group_var)) {
      p <- p + aes(color = !!rlang::sym(group_var),
                   shape = !!rlang::sym(group_var),
                   linetype = !!rlang::sym(group_var)) +
        stat_summary(fun = mean, geom = "line", linewidth = 0.8) + 
        stat_summary(fun = mean, geom = "point", size = 1.8) + 
        scale_color_brewer(palette = "Dark2") + 
        scale_shape_manual(values = 1:length(unique(data[[group_var]]))) + 
        scale_linetype_manual(values = 1:length(unique(data[[group_var]]))) 
    } else {
      p <- p + geom_line(linewidth = 0.8) + geom_point(size = 1.8) 
    }
  } else if (geom_type == "bar") {
    p <- p + geom_bar(stat = "identity", position = "stack", aes(fill = !!rlang::sym(group_var))) +
      { if (isTRUE(y_is_percent)) scale_y_continuous(labels = scales::percent_format(scale = 1)) else scale_y_continuous() } +
      scale_fill_manual(values = c( 
        "bentlerT" = "#E41A1C",    
        "entropy" = "#4DAF4A",    
        "geominT" = "#377EB8",     
        "varimax" = "#FFD700",    
        "unrotated" = "#FF7F00",   
        "quartimax" = "#FFFF33",   
        "equamax" = "#A65628",     
        "cf_varimax" = "#F781BF",  
        "infomaxT" = "#999999",    
        "bifactorT" = "#E6AB02",   
        "varimin" = "#00CED1"      
      ))
  } else if (geom_type == "tile") {
    p <- p + geom_tile(color = "white", linewidth = 1, aes(fill = !!rlang::sym(group_var))) +
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + 
      scale_fill_manual(values = c( 
        "bentlerT" = "#E41A1A",    
        "entropy" = "#4DAF4A",    
        "geominT" = "#377EB8",    
        "varimax" = "#FFD700",    
        "unrotated" = "#FF7F00",   
        "quartimax" = "#FFFF33",   
        "equamax" = "#A65628",     
        "cf_varimax" = "#F781BF",  
        "infomaxT" = "#999999",    
        "bifactorT" = "#E6AB02",  
        "varimin" = "#00CED1"     
      ))
  }
  
  if (!is.null(facet_var1) || !is.null(facet_var2)) {
    rows_facet_arg <- NULL
    cols_facet_arg <- NULL
    
    if (!is.null(facet_var1)) {
      rows_facet_arg <- rlang::inject(ggplot2::vars(!!!rlang::syms(facet_var1)))
    }
    if (!is.null(facet_var2)) {
      cols_facet_arg <- rlang::inject(ggplot2::vars(!!!rlang::syms(facet_var2)))
    }
    
    if (!is.null(rows_facet_arg) && !is.null(cols_facet_arg)) {
      p <- p + ggplot2::facet_grid(rows = rows_facet_arg,
                                   cols = cols_facet_arg,
                                   scales = "free", 
                                   labeller = label_value) 
    } else if (!is.null(rows_facet_arg)) { 
      facet_formula <- as.formula(paste("~", paste(facet_var1, collapse = " + ")))
      p <- p + facet_wrap(facet_formula, scales = "fixed", labeller = label_value)
    } else if (!is.null(cols_facet_arg)) { 
      facet_formula <- as.formula(paste("~", paste(facet_var2, collapse = " + ")))
      p <- p + facet_wrap(facet_formula, scales = "fixed", labeller = label_value)
    }
  }
  
  p <- p + theme_minimal(base_size = 12) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
      axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 8, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), 
      legend.position = "bottom",
      legend.box.spacing = unit(0.01, "cm"), 
      legend.box = "horizontal",
      legend.key.width = unit(0.3, "cm"), 
      legend.key.height = unit(0.3, "cm"),
      legend.spacing.x = unit(0.2, "cm"), 
      legend.spacing.y = unit(0.2, "cm"),
      legend.text = element_text(size = 6), 
      legend.title = element_text(size = 6, face = "bold"), 
      strip.text = element_text(face = "bold", size = 6), 
      panel.spacing = unit(0.5, "cm") 
      
    ) + 
    labs(title = title, x = x_label, y = y_label,
         color = legend_title, shape = legend_title, linetype = legend_title, fill = legend_title)
  
  if (!is.null(legend_ncol) || !is.null(legend_nrow)) {
    guide_args <- list()
    if (!is.null(legend_ncol)) guide_args$ncol <- legend_ncol
    if (!is.null(legend_nrow)) guide_args$nrow <- legend_nrow
    
    my_guide <- do.call(guide_legend, guide_args)
    
    p <- p + guides(
      fill     = my_guide,
      color    = my_guide,
      shape    = my_guide,
      linetype = my_guide
    )
  }
  
  ggsave(save_path, 
         plot = p, 
         device = "tiff",
         width = 17.5, 
         height = 12, 
         units = "cm", 
         dpi = 600,
         compression = "lzw")
  
  return(p)
}

# ===============================================================
# --- ANALYSIS FUNCTIONS --- 
# ===============================================================

#### SCENARIO-BASED ANALYSIS FUNCTIONS ####

run_mnl_analysis <- function(data, measure_type, output_base_dir) {
  cat(paste0("\n--- MNL Analizleri Başlıyor (Ölçüt: ", measure_type, ") ---\n"))
  
  output_dir <- file.path(output_base_dir, "senaryo_bazli", "mnl_analiz", measure_type)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  freq_subset <- data %>% dplyr::filter(measure == measure_type)
  if (nrow(freq_subset) == 0 || length(unique(freq_subset$method)) <= 1) {
    cat(paste0("Uyarı: '", measure_type, "' için MNL modeli kurulamıyor. Yeterli veri yok.\n"))
    return(NULL)
  }
  
  set.seed(123) 
  
  mnl_model <- nnet::multinom(method ~ ds * np + ds * kk_tip + np * kk_tip,
                              data = freq_subset, weights = freq, trace = FALSE)
  
  
  summ <- summary(mnl_model)
  coefs <- summ$coefficients
  ses <- summ$standard.errors 
  
  or_df <- data.frame(
    term = rep(rownames(coefs), each = ncol(coefs)),
    reference_method = rep(colnames(coefs), times = nrow(coefs)),
    OR = as.numeric(exp(t(coefs))),
    LCL = as.numeric(exp(t(coefs - 1.96 * ses))), 
    UCL = as.numeric(exp(t(coefs + 1.96 * ses))), 
    p = as.numeric(2 * (1 - pnorm(abs(t(coefs / ses))))) 
  )
  
  file_name <- paste0("mnl_odds_ratios_for_", measure_type, ".csv")
  write.csv(or_df, file.path(output_dir, file_name), row.names = FALSE)
  
  significant_or_df <- or_df %>%
    dplyr::filter(p < 0.05)
  
  if (nrow(significant_or_df) == 0) {
    cat(paste0("Uyarı: '", measure_type, "' için istatistik olarak anlamlı Odds Oranı bulunamadı. Grafik oluşturulmayacak.\n"))
  } else {
    facet_labels <- c(
      "(Intercept)" = "Intercept",
      "ds" = "Number of variables (p)",
      "np" = "Sample-to-variable ratio (n/p)",
      "kk_tipConstant correlation (0.5)" = "Constant correlation (0.5)",
      "kk_tipConstant correlation (0.8)" = "Constant correlation (0.8)",
      "ds:np" = "p × n/p Interaction",
      "ds:kk_tipConstant correlation (0.5)" = "p × Constant correlation (0.5)",
      "ds:kk_tipConstant correlation (0.8)" = "p × Constant correlation (0.8)",
      "np:kk_tipConstant correlation (0.5)" = "n/p × Constant correlation (0.5)",
      "np:kk_tipConstant correlation (0.8)" = "n/p × Constant correlation (0.8)"
      
    )
    
    p_mnl_plot_direct <- ggplot2::ggplot(significant_or_df, ggplot2::aes(x = OR, y = term)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = LCL, xmax = UCL), height = 0.2, color = "darkblue") + 
      ggplot2::geom_point(size = 3, color = "steelblue") + 
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      ggplot2::facet_wrap(~ reference_method, scales = "free", labeller = ggplot2::as_labeller(facet_labels)) +
      ggplot2::labs(title = paste0("Statistically Significant Odds Ratios (Criterion: ", measure_label(measure_type), ")"), x = "Odds Ratio", y = "Rotation Method") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
                     axis.title = element_text(size = 12, face = "bold"),
                     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                     axis.text.y = element_text(size = 8),
                     strip.text = element_text(face = "bold", size = 9))
    
    ggplot2::ggsave(file.path(output_dir, paste0("mnl_odds_ratios_plot_significant_", measure_type, ".tiff")),
                    plot = p_mnl_plot_direct, width = 15, height = 10, units = "in", dpi = 600, compression = "lzw")
    
    cat(paste0("'", measure_type, "' için istatistiksel olarak anlamlı MNL odds oranları grafiği oluşturuldu.\n"))
  }
  cat(paste0("'", measure_type, "' için MNL analizleri tamamlandı.\n"))
}

run_chi_analysis <- function(data, measure_type, output_base_dir) {
  cat(paste0("\n--- Ki-kare Analizleri Başlıyor (Ölçüt: ", measure_type, ") ---\n"))
  
  output_dir <- file.path(output_base_dir, "senaryo_bazli", "ki_kare_analiz", measure_type)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  winner_subset <- data %>% dplyr::filter(measure == measure_type)
  if (nrow(winner_subset) == 0 || length(unique(winner_subset$winner)) <= 1) {
    cat(paste0("Uyarı: '", measure_type, "' için Ki-kare testi kurulamıyor. Yeterli veri yok.\n"))
    return(NULL) 
  }
  
  winner_subset$ds <- as.factor(winner_subset$ds)
  winner_subset$np <- as.factor(winner_subset$np)
  
  all_chi_results <- data.frame() 
  
  scenario_vars_for_chi <- c("ds", "np", "kk_tip")
  
  for (var in scenario_vars_for_chi) {
    contingency_table <- table(winner_subset$winner, winner_subset[[var]])
    
    x_label_map <- c(
      ds = "Number of Variables",
      np = "Sample-to-Variable Ratio",
      kk_tip = "Correlation Type"
    )
    
    measure_labels <- c(
      "max" = "Maximum loading",
      "min" = "Minimum loading",
      "range_max" = "Maximum loading range"
    )
    
    chi_test_result <- tryCatch({
      chisq.test(contingency_table)
    }, error = function(e) {
      cat(paste0("Hata: '", measure_type, "' için '", var, "' değişkeniyle Ki-kare testi yapılamadı. Hata: ", e$message, "\n"))
      return(NULL)
    })
    
    if (!is.null(chi_test_result)) {
      if (all(dim(contingency_table) > 1) && all(rowSums(contingency_table) > 0) && all(colSums(contingency_table) > 0)) {
        cramers_v <- vcd::assocstats(contingency_table)$cramer
      } else {
        cramers_v <- NA 
        cat(paste0("Uyarı: '", measure_type, "' için '", var, "' değişkeniyle Cramér's V hesaplanamadı (sıfır marjinal toplamlar).\n"))
      }
      
      chi_results_df <- data.frame(
        Measure = measure_type, 
        Test = paste0("Winner vs ", var),
        ChiSq_Stat = ifelse(!is.null(chi_test_result$statistic), chi_test_result$statistic, NA),
        Df = ifelse(!is.null(chi_test_result$parameter), chi_test_result$parameter, NA),
        p_value = ifelse(!is.null(chi_test_result$p.value), chi_test_result$p.value, NA),
        Cramers_V = cramers_v
      )
      all_chi_results <- rbind(all_chi_results, chi_results_df)
    }
    
    plot_data <- winner_subset %>%
      dplyr::group_by(!!sym(var), winner) %>% 
      dplyr::summarise(count = dplyr::n(), .groups = "drop_last") %>% 
      dplyr::mutate(percentage = count / sum(count) * 100) %>% 
      dplyr::ungroup() 
    
    plot_academic(
      data = plot_data,
      x_var = var,
      y_var = "percentage",
      group_var = "winner", 
      geom_type = "bar", 
      title = paste0("Distribution of Winning Methods (Criterion: ", measure_label(measure_type), ")"),
      x_label = x_label_map[[var]], 
      y_label = "Percentage",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_dir, paste0("ki_kare_gorsel_", measure_type, "_", var, ".tiff")),
      width = 15, height = 10, units = "in" 
    )
  }
  
  write.csv(all_chi_results, file.path(output_dir, paste0("ki_kare_sonuclari_", measure_type, ".csv")), row.names = FALSE)
  cat(paste0("'", measure_type, "' için Ki-kare testleri tamamlandı.\n"))
  
  plot_data_faceted <- winner_subset %>%
    dplyr::group_by(ds, np, kk_tip, winner) %>% 
    dplyr::summarise(count = dplyr::n(), .groups = "drop_last") %>% 
    dplyr::mutate(percentage = count / sum(count) * 100) %>% 
    dplyr::ungroup() 
  
  if (nrow(plot_data_faceted) == 0) {
    cat(paste0("Uyarı: '", measure_type, "' için facet'li kazanan dağılım grafiği oluşturulamıyor. Veri yok.\n"))
    return(NULL) 
  }
  
  plot_academic(
    data = plot_data_faceted, 
    x_var = "winner",
    y_var = "percentage",
    group_var = "winner", 
    facet_var1 = "kk_tip",
    facet_var2 = "ds", 
    geom_type = "bar", 
    title = paste0("Percentages of Winning Methods (Criterion: ", measure_label(measure_type), ")\nFactors: Correlation Structure and Number of Variables"),
    y_label = "Percentage",
    legend_title = "Rotation Method",
    legend_nrow = 2,
    save_path = file.path(output_dir, paste0("winner_distribution_faceted_kk_ds_", measure_type, ".tiff")),
    width = 18, 
    height = 12, 
    units = "in" 
  )
  cat(paste0("'", measure_type, "' için detaylı facetli KK Tipi vs DS grafiği kaydedildi.\n"))
  
  plot_academic(
    data = plot_data_faceted,
    x_var = "winner",
    y_var = "percentage",
    group_var = "winner", 
    facet_var1 = "kk_tip",
    facet_var2 = "np", 
    geom_type = "bar", 
    title = paste0("Percentages of Winning Methods (Criterion: ", measure_label(measure_type), ")\nFactors: Correlation Structure and Number of Observations"),
    y_label = "Percentage",
    legend_title = "Rotation Method",
    legend_nrow = 2,
    save_path = file.path(output_dir, paste0("winner_distribution_faceted_kk_np_", measure_type, ".tiff")),
    width = 20, 
    height = 12, 
    units = "in" 
  )
  cat(paste0("'", measure_type, "' için detaylı facetli KK Tipi vs NP grafiği kaydedildi.\n"))
  
  plot_academic(
    data = plot_data_faceted,
    x_var = "winner",
    y_var = "percentage",
    group_var = "winner", 
    facet_var1 = "ds", 
    facet_var2 = "np", 
    geom_type = "bar", 
    title = paste0("Percentages of Winning Methods (Criterion: ", measure_label(measure_type), ")\nFactors: Number of Variables and Sample Size"),
    y_label = "Percentage",
    legend_title = "Rotation Method",
    legend_nrow = 2,
    save_path = file.path(output_dir, paste0("winner_distribution_faceted_ds_np_", measure_type, ".tiff")),
    width = 22, 
    height = 15, 
    units = "in" 
  )
  cat(paste0("'", measure_type, "' için detaylı facetli DS vs NP grafiği kaydedildi.\n"))
  
}

analyze_by_group <- function(data, group_var, group_var_label, measure_type, output_base_dir_for_measure) {
  group_var_sym <- rlang::sym(group_var)
  
  summary_df <- data %>%
    dplyr::group_by(!!group_var_sym, winner) %>% 
    dplyr::summarise(count = n(), .groups = "drop_last") %>% 
    dplyr::mutate(percentage = count / sum(count) * 100) %>% 
    dplyr::ungroup() 
  
  if (nrow(summary_df) == 0) {
    cat(paste0("Uyarı: '", group_var, "' bazında kazanma oranları grafiği oluşturulamıyor. Veri yok.\n"))
    return(NULL) 
  }
  
  write.csv(summary_df, file.path(output_base_dir_for_measure, paste0("kazanma_oranlari_by_", group_var, ".csv")), row.names = FALSE)
  
  plot_academic(
    data = summary_df,
    x_var = "winner", 
    y_var = "percentage",
    group_var = "winner", 
    facet_var1 = group_var,
    geom_type = "bar", 
    title = paste0("Rotation Method Performance Under the Influence of ", group_var_label, "(Criterion: ", measure_label(measure_type), ")"),
    y_label = "Winning Percentage (%)",
    legend_title = "Rotation Method",
    legend_nrow = 2,
    save_path = file.path(output_base_dir_for_measure, paste0("kazanma_oranlari_gorsel_", group_var, ".tiff")),
    width = 15, height = 10, units = "in" 
  )
  cat(paste0(group_var, " değişkenine göre kazanma oranları analizi tamamlandı.\n"))
}

plot_most_frequent_winner_by_scenario <- function(data_filtered_by_measure, measure_type, output_base_dir_for_measure) {
  output_file_path <- file.path(output_base_dir_for_measure, paste0("en_cok_kazanan_yontem_", measure_type, ".tiff"))
  
  most_frequent_winner <- data_filtered_by_measure %>%
    dplyr::group_by(ds, np, kk_tip) %>% 
    dplyr::count(winner, name = "winner_count") %>%
    dplyr::slice_max(winner_count, n = 1, with_ties = FALSE) %>% 
    dplyr::ungroup()
  
  if (nrow(most_frequent_winner) == 0) {
    cat(paste0("Uyarı: Ölçüt '", measure_type, "' için en çok kazanan yöntem ısı haritası oluşturulamıyor. Veri yok.\n"))
    return(NULL)
  }
  
  p <- ggplot(most_frequent_winner, aes(x = factor(ds), y = factor(np))) +
    geom_tile(color = "white", linewidth = 1, aes(fill = winner)) +
    facet_wrap(~ kk_tip, ncol = 3, scales = "fixed") + 
    scale_fill_manual(values = c( 
      "bentlerT" = "#E41A1A",    
      "entropy" = "#4DAF4A",     
      "geominT" = "#377EB8",     
      "varimax" = "#FFD700",     
      "unrotated" = "#FF7F00",        
      "quartimax" = "#FFFF33",   
      "equamax" = "#A65628",     
      "cf_varimax" = "#F781BF",  
      "infomaxT" = "#999999",    
      "bifactorT" = "#E6AB02",   
      "varimin" = "#00CED1"      
    )) +
    labs(
      title = "Most Frequent Winning Rotation Method by Scenario",
      x = "Number of Variables", 
      y = "Sample-to-Variable Ratio (n/p)",
      fill = "Rotation Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10), 
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(face = "bold", size = 10), 
      panel.spacing = unit(0.5, "cm"),
      panel.background = element_rect(fill = NA, colour = NA), 
      plot.background = element_rect(fill = NA, colour = NA) 
    )+
    guides(
      fill = guide_legend(nrow = 2),
      color = guide_legend(nrow = 2),
      shape = guide_legend(nrow = 2),
      linetype = guide_legend(nrow = 2)
    )
  
  ggsave(output_file_path, p, width = 15, height = 10, units = "in", dpi = 600, compression = "lzw", bg="transparent") 
  cat(paste0("Ölçüt '", measure_type, "' için en çok kazanan yöntem grafiği kaydedildi: ", output_file_path, "\n"))
}

#### VARIABLE-BASED ANALYSIS FUNCTIONS ####

run_variable_based_analysis <- function(data, measure_type, base_output_dir) {
  cat(paste0("\n--- DEĞİŞKEN BAZLI ANALİZLER Başlıyor (Ölçüt: ", measure_type, ") ---\n"))
  
  output_dir <- file.path(base_output_dir, "degisken_bazli", measure_type)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  current_winner_subset <- data %>% dplyr::filter(measure == measure_type)
  
  if (nrow(current_winner_subset) == 0) {
    cat(paste0("Uyarı: '", measure_type, "' için değişken bazlı analiz yapılamıyor. Yeterli veri yok.\n"))
    return(NULL)
  }
  
  if (!"Variable" %in% colnames(current_winner_subset)) {
    cat(paste0("Hata: 'Variable' sütunu veri setinde bulunamadı. Değişken bazlı analiz atlanıyor. Lütfen simülasyon çıktısını ('senaryo_winner.rds') kontrol edin ve Simülasyon.txt dosyasındaki 'temp_dt_to_add' kısmını düzeltin.\n"))
    return(NULL)
  }
  
  current_winner_subset_for_plots <- current_winner_subset %>%
    dplyr::mutate(
      ds = as.factor(ds),
      np = as.factor(np)
    )
  
  variable_winner_summary <- current_winner_subset %>%
    dplyr::group_by(Variable, winner) %>% 
    dplyr::summarise(count = dplyr::n(), .groups = "drop_last") %>% 
    dplyr::mutate(percentage = count / sum(count) * 100) %>% 
    dplyr::ungroup() 
  
  if (nrow(variable_winner_summary) == 0) {
    cat(paste0("Uyarı: '", measure_type, "' için değişken bazlı genel kazanan dağılım grafiği oluşturulamıyor. Veri yok.\n"))
    return(NULL) 
  }
  
  write.csv(variable_winner_summary, file.path(output_dir, paste0("degisken_genel_kazanan_dagilimi_", measure_type, ".csv")), row.names = FALSE)
  
  plot_academic(
    data = variable_winner_summary,
    x_var = "winner",
    y_var = "percentage",
    group_var = "winner",
    facet_var1 = "Variable",
    geom_type = "bar",
    title = paste0("Distributions of Winning Methods by Variable (Criterion: ", measure_label(measure_type), ")"),
    y_label = "Percentage (%)",
    legend_title = "Rotation Method",
    legend_nrow = 2,
    save_path = file.path(output_dir, paste0("degisken_genel_kazanan_dagilimi_barplot_", measure_type, ".tiff")),
    width = 15, height = 10, units = "in"
  )
  cat(paste0("'", measure_type, "' için her değişkenin genel kazanan dağılımı analizi tamamlandı.\n"))
  
  most_frequent_winner_by_var_scenario <- current_winner_subset_for_plots %>%
    dplyr::group_by(Variable, ds, np, kk_tip) %>%
    dplyr::count(winner, name = "winner_count") %>%
    dplyr::slice_max(winner_count, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  unique_kk_tips_var <- unique(most_frequent_winner_by_var_scenario$kk_tip)
  for (k_tip_var in unique_kk_tips_var) {
    subset_var_data_for_tile <- most_frequent_winner_by_var_scenario %>%
      dplyr::filter(kk_tip == k_tip_var)
    
    if (nrow(subset_var_data_for_tile) == 0) {
      cat(paste0("Uyarı: Ölçüt '", measure_type, "', KK Tipi '", k_tip_var, "' için değişken bazlı ısı haritası oluşturulamıyor. Veri yok.\n"))
      next
    }
    
    plot_academic(
      data = subset_var_data_for_tile,
      x_var = "Variable",
      y_var = "ds",
      group_var = "winner",
      facet_var1 = "np",
      geom_type = "tile",
      title = paste0("Most Frequent Winning Method (Criterion: ", measure_label(measure_type), ", Correlation: ", k_tip_var, ")"),
      x_label = "Variable",
      y_label = "Number of variables",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_dir, paste0("degisken_ds_kazanan_tileplot_", measure_type, "_kk", gsub("\\.", "p", k_tip_var), ".tiff")),
      width = 30, height = 20, units = "in"
    )
    
    plot_academic(
      data = subset_var_data_for_tile,
      x_var = "Variable",
      y_var = "np",
      group_var = "winner",
      facet_var1 = "ds",
      geom_type = "tile",
      title = paste0("Most Frequent Winning Method (Criterion: ", measure_labels[m_type], ", Correlation: ", k_tip_var, ")"),
      x_label = "Variable",
      y_label = "Sample-to-Variable Ratio (n/p)",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_dir, paste0("degisken_np_kazanan_tileplot_", measure_type, "_kk", gsub("\\.", "p", k_tip_var), ".tiff")),
      width = 30, height = 20, units = "in"
    )
  }
  cat(paste0("\n--- Değişken Bazlı Analizler Tamamlandı (Ölçüt: ", measure_type, ") ---\n"))
}

#### ADDITIONAL METRICS ANALYSIS FUNCTIONS ####

run_ek_metrikler_analysis <- function(data, base_output_dir) {
  output_dir <- file.path(base_output_dir, "ek_metrikler", "yontem_bazli")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  conv_summary <- data %>%
    dplyr::filter(!is.na(converged_ratio), is.finite(converged_ratio)) %>%
    dplyr::group_by(method) %>% 
    dplyr::summarise(converged_mean = mean(converged_ratio, na.rm = TRUE), .groups = "drop")
  
  if (nrow(conv_summary) == 0) {
    cat(paste0("Uyarı: Konverjans oranı grafiği oluşturulamıyor. Veri yok.\n"))
  } else {
    write.csv(conv_summary, file.path(output_dir, "converged_ratio_summary.csv"), row.names = FALSE)
    plot_academic(
      data = conv_summary,
      y_is_percent = FALSE,
      x_var = "method",
      y_var = "converged_mean",
      group_var = "method", 
      geom_type = "bar", 
      title = "Average Convergence Rate by Method",
      y_label = "Average Convergence Rate",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_dir, "converged_ratio_barplot.tiff"),
      width = 15, height = 10, units = "in" 
    )
  }
  
  ssi_summary <- data %>%
    dplyr::filter(!is.na(mean_SSI), is.finite(mean_SSI)) %>%
    dplyr::group_by(method) %>% 
    dplyr::summarise(mean_SSI = mean(mean_SSI, na.rm = TRUE), .groups = "drop")
  
  if (nrow(ssi_summary) == 0) {
    cat(paste0("Uyarı: SSI grafiği oluşturulamıyor. Veri yok.\n"))
  } else {
    write.csv(ssi_summary, file.path(output_dir, "mean_SSI_summary.csv"), row.names = FALSE)
    plot_academic(
      data = ssi_summary,
      y_is_percent = FALSE,
      x_var = "method",
      y_var = "mean_SSI",
      group_var = "method", 
      geom_type = "bar", 
      title = "Average Simple Structure Index (SSI) by Method",
      y_label = "Average Simple Structure Index",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_dir, "mean_SSI_barplot.tiff"),
      width = 15, height = 10, units = "in" 
    )
  }
  
  prim_summary <- data %>%
    dplyr::filter(!is.na(mean_primary), is.finite(mean_primary)) %>%
    dplyr::group_by(method) %>% 
    dplyr::summarise(mean_primary = mean(mean_primary, na.rm = TRUE), .groups = "drop")
  
  if (nrow(prim_summary) == 0) {
    cat(paste0("Uyarı: Birincil yükler grafiği oluşturulamıyor. Veri yok.\n"))
  } else {
    write.csv(prim_summary, file.path(output_dir, "mean_primary_summary.csv"), row.names = FALSE)
    plot_academic(
      data = prim_summary,
      y_is_percent = FALSE,
      x_var = "method",
      y_var = "mean_primary",
      group_var = "method", 
      geom_type = "bar", 
      title = "Average Primary Factor Loading by Method",
      y_label = "Average Primary Loading Value",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_dir, "mean_primary_barplot.tiff"),
      width = 15, height = 10, units = "in" 
    )
  }
  
  if ("mean_cross" %in% names(data)) {
    cross_summary <- data %>%
      dplyr::filter(!is.na(mean_cross), is.finite(mean_cross)) %>%
      dplyr::group_by(method) %>% 
      dplyr::summarise(mean_cross = mean(mean_cross, na.rm = TRUE), .groups = "drop")
    
    if (nrow(cross_summary) == 0) {
      cat(paste0("Uyarı: Çapraz yükler grafiği oluşturulamıyor. Veri yok.\n"))
    } else {
      write.csv(cross_summary, file.path(output_dir, "mean_cross_summary.csv"), row.names = FALSE)
      plot_academic(
        data = cross_summary,
        y_is_percent = FALSE,
        x_var = "method",
        y_var = "mean_cross",
        group_var = "method", 
        geom_type = "bar", 
        title = "Average Cross-Loading by Method",
        y_label = "Average Cross-Loading Value",
        legend_title = "Rotation Method",
        legend_nrow = 2,
        save_path = file.path(output_dir, "mean_cross_barplot.tiff"),
        width = 15, height = 10, units = "in" 
      )
    }
  }
  
  library(dplyr)
  library(readr)
  
  conv_df <- read_csv(file.path(output_dir, "converged_ratio_summary.csv")) %>%
    rename(value = converged_mean) %>%
    mutate(metric = "Converged ratio")
  
  cross_df <- read_csv(file.path(output_dir, "mean_cross_summary.csv")) %>%
    rename(value = mean_cross) %>%
    mutate(metric = "Average cross-loading")
  
  prim_df <- read_csv(file.path(output_dir, "mean_primary_summary.csv")) %>%
    rename(value = mean_primary) %>%
    mutate(metric = "Average primary loading")
  
  ssi_df <- read_csv(file.path(output_dir, "mean_SSI_summary.csv")) %>%
    rename(value = mean_SSI) %>%
    mutate(metric = "Average simple structure index (SSI)")
  
  all_metrics_df <- bind_rows(conv_df, cross_df, prim_df, ssi_df)
  
  plot_academic(
    data = all_metrics_df,
    x_var = "method",
    y_var = "value",
    y_is_percent = FALSE,
    group_var = "method",
    facet_var1 = "metric",
    geom_type = "bar",
    title = "Comparison of Additional Metrics by Rotation Method",
    y_label = "Ratio",
    legend_title = "Rotation Method",
    legend_nrow = 2,
    save_path = file.path(output_dir, "ek_metrikler_yontem_bazli_facet.tiff"),
    width = 18, height = 12, units = "in"
  )
  cat("Ek metrikler facet'li grafik olarak kaydedildi.\n")
  
  
  cat("Ek metrikler (yöntem bazlı) analiz edildi, tablolar ve grafikler kaydedildi.\n")
}

run_ek_metrikler_analysis_by_rotation <- function(data, base_output_dir) {
  output_dir <- file.path(base_output_dir, "ek_metrikler", "donusum_bazli")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  metrics <- c("converged_ratio", "mean_SSI", "mean_primary", "mean_cross")
  metric_names_tr <- c(
    converged_ratio = "Convergence Rate",
    mean_SSI = "Average Simple Structure Index (SSI)",
    mean_primary = "Average Primary Loading",
    mean_cross = "Average Cross-Loading"
  )
  
  for (metric in metrics) {
    if (metric %in% names(data) && any(!is.na(data[[metric]]) & is.finite(data[[metric]]))) {
      cat(paste0("\n--- ", metric, " için analizler başlıyor ---\n"))
      
      summary_table <- data %>%
        dplyr::filter(!is.na(!!sym(metric)), is.finite(!!sym(metric))) %>%
        dplyr::group_by(rotation_method, ds, np, kk_tip) %>%
        dplyr::summarise(
          mean_value = mean(!!sym(metric), na.rm = TRUE),
          sd_value = sd(!!sym(metric), na.rm = TRUE),
          .groups = "drop"
        )
      
      if (nrow(summary_table) == 0) {
        cat(paste0("Uyarı: '", metric, "' için döndürme bazlı genel özet tablosu oluşturulamıyor. Veri yok.\n"))
        next 
      }
      
      write.csv(summary_table, file.path(output_dir, paste0(metric, "_summary_by_rotation.csv")), row.names = FALSE)
      
      clean_data <- data %>%
        dplyr::filter(!is.na(!!sym(metric)), is.finite(!!sym(metric))) %>%
        dplyr::mutate(
          rotation_method = forcats::fct_drop(rotation_method),
          ds = forcats::fct_drop(as.factor(ds)), 
          np = forcats::fct_drop(as.factor(np)), 
          kk_tip = forcats::fct_drop(kk_tip)
        )
      
      if (nrow(clean_data) < 2 || length(unique(clean_data$rotation_method)) < 2) {
        cat(paste0("Uyarı: '", metric, "' için ANOVA analizi yapılamıyor (yeterli veri veya rotation_method seviyesi yok). Analiz atlandı.\n"))
        next
      }
      
      aov_model <- NULL
      tryCatch({
        aov_formula_full <- reformulate(
          termlabels = c("rotation_method", "ds", "np", "kk_tip",
                         "rotation_method:ds", "rotation_method:np", "rotation_method:kk_tip"),
          response = metric
        )
        temp_aov_full <- aov(aov_formula_full, data = clean_data)
        if (is.null(temp_aov_full$df.residual) || temp_aov_full$df.residual <= 0) {
          stop("Tam modelin artık serbestlik derecesi sıfır veya negatif. Basit modele geçiliyor.")
        }
        aov_model <- temp_aov_full
      }, error = function(e) {
        cat(paste0("Hata (Tam Model): ", metric, " için ANOVA modeli kurulamadı veya artık SD < 1. Hata mesajı: ", e$message, "\n"))
        cat("   -> Sadece ana etkilerle daha basit bir model deneniyor.\n")
        
        tryCatch({
          aov_formula_main <- reformulate(
            termlabels = c("rotation_method", "ds", "np", "kk_tip"),
            response = metric
          )
          temp_aov_main <- aov(aov_formula_main, data = clean_data)
          if (is.null(temp_aov_main$df.residual) || temp_aov_main$df.residual <= 0) {
            stop("Basit modelin de artık serbestlik derecesi sıfır veya negatif.")
          }
          aov_model <- temp_aov_main
        }, error = function(e_main) {
          cat(paste0("Hata (Basit Model): ", metric, " için basit ANOVA modeli de kurulamadı. Hata mesajı: ", e_main$message, "\n"))
        })
      })
      
      if (is.null(aov_model)) {
        cat(paste0("Uyarı: '", metric, "' için hiçbir ANOVA modeli başarıyla kurulamadı. Post-hoc atlanıyor.\n"))
        next 
      }
      
      if (!is.null(aov_model)) {
        summary_aov <- summary(aov_model)
        capture.output(summary_aov, file = file.path(output_dir, paste0(metric, "_anova_results.txt")))
        cat(paste0(metric, " için ANOVA sonuçları kaydedildi.\n"))
        
        anova_table <- as.data.frame(summary_aov[[1]])
        
        rot_main_p <- NA
        if ("rotation_method" %in% rownames(anova_table)) {
          rot_main_p <- anova_table["rotation_method", "Pr(>F)"]
        }
        
        significant_interactions <- FALSE
        interaction_terms_in_anova <- rownames(anova_table)[grepl(":", rownames(anova_table))]
        
        for (int_term in interaction_terms_in_anova) {
          if (grepl("rotation_method", int_term) && !is.na(anova_table[int_term, "Pr(>F)"]) && anova_table[int_term, "Pr(>F)"] < 0.05) {
            significant_interactions <- TRUE
            break
          }
        }
        
        if (!is.na(rot_main_p) && rot_main_p < 0.05 && !significant_interactions) {
          cat(paste0("   -> rotation_method ana etkisi anlamlı (p = ", round(rot_main_p, 4), "). Tukey HSD ile post-hoc test yapılıyor.\n"))
          em_rot <- emmeans::emmeans(aov_model, ~ rotation_method)
          post_hoc_result <- pairs(em_rot, adjust = "tukey")
          capture.output(summary(post_hoc_result), file = file.path(output_dir, paste0(metric, "_posthoc_rotation_method_main.txt")))
          cat(paste0("   -> ", metric, " için rotation_method ana etkisi post-hoc sonuçları kaydedildi.\n"))
        } else if (significant_interactions) {
          cat(paste0("   -> rotation_method içeren anlamlı etkileşimler bulundu. Koşullu post-hoc testler yapılıyor.\n"))
          
          for (int_term in interaction_terms_in_anova) {
            if (grepl("rotation_method", int_term) && !is.na(anova_table[int_term, "Pr(>F)"]) && anova_table[int_term, "Pr(>F)"] < 0.05) {
              factors_in_interaction <- unlist(strsplit(int_term, ":"))
              conditioning_factors <- setdiff(factors_in_interaction, "rotation_method")
              
              if (length(conditioning_factors) > 0) {
                specs_str_for_emmeans <- paste0("~ rotation_method | ", paste(conditioning_factors, collapse = " * "))
                cat(paste0("     -> Etkileşim: ", int_term, " için koşullu post-hoc yapılıyor (Specs: ", specs_str_for_emmeans, ").\n"))
                em_conditional <- emmeans::emmeans(aov_model, as.formula(specs_str_for_emmeans))
                post_hoc_result_conditional <- pairs(em_conditional, adjust = "tukey")
                file_name_suffix <- gsub(":", "_", int_term)
                capture.output(summary(post_hoc_result_conditional), file = file.path(output_dir, paste0(metric, "_posthoc_", file_name_suffix, ".txt")))
                cat(paste0("     -> ", metric, " için ", int_term, " etkileşimi post-hoc sonuçları kaydedildi.\n"))
              } else {
                cat(paste0("     -> Uyarı: Etkileşim terimi '", int_term, "' için koşullama faktörü bulunamadı. Koşullu post-hoc atlanıyor.\n"))
              }
            }
          }
        } else {
          cat(paste0("   -> rotation_method ana etkisi veya etkileşimleri anlamlı değil (p = ", ifelse(is.na(rot_main_p), "NA", round(rot_main_p, 4)), "). Post-hoc test yapılmadı.\n"))
        }
      }
      
      plot_academic(
        data = data,
        y_is_percent = FALSE,
        x_var = "rotation_method",
        y_var = metric,
        group_var = "rotation_method",
        facet_var1 = "kk_tip",
        facet_var2 = c("ds", "np"),
        geom_type = "bar",
        title = paste0("Distribution of ", metric_names_tr[metric], " Values by Scenario"),
        y_label = "Convergence rate",
        legend_title = "Rotation Method",
        legend_nrow = 2,
        save_path = file.path(output_dir, paste0(metric, "_barplot_by_rotation_faceted.tiff")),
        width = 18, height = 12, units = "in"
      ) +
        scale_fill_manual(values = my_colors, drop = FALSE) +
        facet_grid(rows = vars(kk_tip), cols = vars(ds, np), drop = FALSE)
      
      cat(paste0(metric, " için tüm analizler tamamlandı.\n"))
    } else {
      cat(paste0("Uyarı: '", metric, "' değişkeni bulunamadı veya yeterli veri yok. Analiz atlandı.\\n"))
    }
  }
  
  cat("\nEk metrikler (döndürme yöntemlerine göre) analizi tamamlandı.\\n")
}


# ===============================================================
# --- RUN MAIN ANALYSIS ---
# ===============================================================

base_output_dir <- file.path(getwd(), "analiz_sonuclari")
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

yuk_ozellikleri <- c("max", "min", "range_max")
all_summary_winner <- list()  

for (m_type in yuk_ozellikleri) {
  
  cat(paste0("\n===================================================================\n"))
  cat(paste0("--- ÖLÇÜT: '", m_type, "' İÇİN ANALİZLER BAŞLIYOR ---\n"))
  cat(paste0("===================================================================\n"))
  
  measure_labels <- c(
    "max" = "Maximum loading",
    "min" = "Minimum loading",
    "range_max" = "Maximum loading range"
  )
  
  current_winner_subset <- winner_df %>% dplyr::filter(measure == m_type)
  
  cat(paste0("\n--- SENARYO BAZLI ANALİZLER (Ölçüt: ", m_type, ") ---\n"))
  
  output_summary_dir <- file.path(base_output_dir, "senaryo_bazli", "genel_ozet", m_type)
  dir.create(output_summary_dir, showWarnings = FALSE, recursive = TRUE)
  
  summary_winner <- current_winner_subset %>%
    dplyr::group_by(measure, winner) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(percentage = count / sum(count) * 100) 
  
  if (nrow(summary_winner) == 0) {
    cat(paste0("Uyarı: '", m_type, "' için genel kazanma oranları grafiği oluşturulamıyor. Veri yok.\n"))
  } else {
    write.csv(summary_winner, file.path(output_summary_dir, paste0("winner_summary_", m_type, ".csv")), row.names = FALSE)
    cat(paste0("Genel özet tablo oluşturuldu (winner_summary_", m_type, ".csv).\n"))
    
    summary_winner$measure_type <- m_type
    all_summary_winner[[m_type]] <- summary_winner
    
    plot_academic(
      data = summary_winner,
      x_var = "winner",
      y_var = "percentage",
      group_var = "winner",
      geom_type = "bar",
      title = paste0("Overall Success Percentages of Methods (Criterion: ", measure_label(m_type), ")"),
      y_label = "Winning Percentage (%)",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(output_summary_dir, paste0("genel_kazanma_oranlari_", m_type, ".tiff")),
      width = 15, height = 10, units = "in"
    )
    cat(paste0("'", m_type, "' için yöntemlerin genel kazanma oranları analizi tamamlandı.\n"))
  }
  
  if (length(all_summary_winner) > 0) {
    all_summary_winner_df <- dplyr::bind_rows(all_summary_winner)
    all_summary_winner_df$measure_type <- factor(
      all_summary_winner_df$measure_type,
      levels = c("min", "max", "range_max"),
      labels = c("Criterion: Minimum loading",
                 "Criterion: Maximum loading",
                 "Criterion: Maximum loading range")
    )
    
    genel_ozet_dir <- file.path(base_output_dir, "senaryo_bazli", "genel_ozet")
    dir.create(genel_ozet_dir, showWarnings = FALSE, recursive = TRUE)
    
    plot_academic(
      data = all_summary_winner_df,
      x_var = "winner",
      y_var = "percentage",
      group_var = "winner",
      facet_var1 = "measure_type",
      geom_type = "bar",
      title = "Overall Success Percentages of Methods Across All Criteria",
      y_label = "Winning Percentage (%)",
      legend_title = "Rotation Method",
      legend_nrow = 2,
      save_path = file.path(genel_ozet_dir, "genel_kazanma_oranlari_facet.tiff"),
      width = 15, height = 10, units = "in"
    ) 
    
    cat("Tüm ölçütler için facet'li genel kazanma oranı grafiği kaydedildi.\n")
  }
  
  run_mnl_analysis(freq_df, m_type, base_output_dir)
  
  run_chi_analysis(winner_df, m_type, base_output_dir) 
  
  cat(paste0("\n--- Yöntemlerin Kazanma Oranları Analizleri Başlıyor (Ölçüt: ", m_type, ") ---\n"))
  
  scenario_output_dir <- file.path(base_output_dir, "senaryo_bazli", m_type) 
  dir.create(scenario_output_dir, showWarnings = FALSE, recursive = TRUE)
  analyze_by_group(current_winner_subset, "ds", "Number of variable", m_type, scenario_output_dir)
  analyze_by_group(current_winner_subset, "np", "Sample-to-variable ratio(n/p)", m_type, scenario_output_dir)
  analyze_by_group(current_winner_subset, "kk_tip", "Correlation type", m_type, scenario_output_dir)
  
  cat(paste0("\n--- Yöntemlerin Kazanma Oranları Analizleri Tamamlandı (Ölçüt: ", m_type, ") ---\n"))
  
  plot_most_frequent_winner_by_scenario(current_winner_subset, m_type, scenario_output_dir)
  
  cat(paste0("\n--- SENARYO BAZLI ANALİZLER TAMAMLANDI (Ölçüt: ", m_type, ") ---\n"))
  
  cat(paste0("\n--- DEĞİŞKEN BAZLI ANALİZLER (Ölçüt: ", m_type, ") ---\n"))
  run_variable_based_analysis(winner_df, m_type, base_output_dir)
  cat(paste0("\n--- DEĞİŞKEN BAZLI ANALİZLER TAMAMLANDI (Ölçüt: ", m_type, ") ---\n"))
  
  cat(paste0("\n===================================================================\n"))
  cat(paste0("--- ÖLÇÜT: '", m_type, "' İÇİN TÜM ANALİZLER TAMAMLANDI ---\n"))
  cat(paste0("===================================================================\n"))
}

freq_df_rotated <- freq_df %>%
  dplyr::mutate(rotation_method = sapply(strsplit(as.character(method), "_"), `[`, 1))

if ("cf" %in% levels(as.factor(freq_df_rotated$rotation_method))) { 
  levels(freq_df_rotated$rotation_method)[levels(freq_df_rotated$rotation_method) == "cf"] <- "cf_varimax"
  cat("Döndürme yöntemi 'cf' -> 'cf_varimax' olarak düzeltildi.\n")
}
freq_df_rotated$rotation_method <- as.factor(freq_df_rotated$rotation_method) 

cat(paste0("\n--- EK METRİKLER ANALİZLERİ BAŞLIYOR ---\n"))

run_ek_metrikler_analysis(freq_df, base_output_dir)

run_ek_metrikler_analysis_by_rotation(freq_df_rotated, base_output_dir)
cat(paste0("\n--- EK METRİKLER ANALİZLERİ TAMAMLANDI ---\n"))

sim_runs <- 10000

scenario_totals <- freq_df %>%
  group_by(scenario = interaction(ds, np, kk_tip, measure, drop = TRUE)) %>%
  summarise(total = sum(freq), .groups = "drop") %>%
  mutate(reps = total / sim_runs,
         reps_int = as.integer(reps),
         is_integer = (reps_int == reps) & (reps > 0))

reps_map <- scenario_totals %>%
  dplyr::transmute(scenario = as.character(scenario), reps_int)

freq_df_adj <- freq_df %>%
  mutate(scenario = as.character(interaction(ds, np, kk_tip, measure, drop = TRUE))) %>%
  left_join(reps_map, by = "scenario") %>%
  mutate(
    reps_int = ifelse(is.na(reps_int), 1L, reps_int),
    freq_adj = freq / reps_int
  )

check <- freq_df_adj %>%
  group_by(scenario) %>%
  summarise(total_adj = sum(freq_adj), .groups = "drop") %>%
  mutate(matches = total_adj == sim_runs)

print(head(check, 10))

sim_runs <- 10000

interpretability <- freq_df_adj %>%
  group_by(scenario) %>%
  summarise(max_freq = max(freq_adj), .groups = "drop") %>%
  mutate(max_prop = max_freq / sim_runs) %>%
  summarise(mean_max_prop = mean(max_prop)) %>%
  pull(mean_max_prop)

stability <- winner_df %>%
  count(winner) %>%
  mutate(prop = n / sum(n)) %>%
  summarise(max_prop = max(prop)) %>%
  pull(max_prop)

ISI <- (interpretability + stability) / 2

cat(sprintf("\nInterpretability (0-1): %.4f\n", interpretability))
cat(sprintf("Stability (0-1): %.4f\n", stability))
cat(sprintf("ISI (0-1): %.4f\n", ISI))

normalized_entropy <- function(freqs) {
  freqs <- as.numeric(freqs)
  if (length(freqs) <= 1) return(0)
  p <- freqs / sum(freqs)
  p <- p[p > 0]
  H <- -sum(p * log2(p))
  H_norm <- H / log2(length(freqs))
  as.numeric(H_norm)
}

entropy_df <- freq_df_adj %>%
  group_by(scenario) %>%
  summarise(entropy = normalized_entropy(freq_adj), total = sum(freq_adj), .groups = "drop")

isi_df <- data.frame(
  Metric = c("Interpretability", "Stability", "ISI"),
  Value = c(interpretability, stability, ISI)
)

p_isi <- ggplot(isi_df, aes(x = Metric, y = Value, fill = Metric)) +
  geom_col(width = 0.6, show.legend = FALSE, na.rm = TRUE) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Interpretability-Stability Index (ISI)", x = "Metric", y = "Value (percent)") +
  theme_minimal(base_size = 14)

ggsave(file.path(base_output_dir, "ISI_barplot.tiff"), plot = p_isi,
       width = 9, height = 6, units = "in", dpi = 600, device = "tiff")

p_entropy <- ggplot(entropy_df, aes(x = entropy)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "#E74C3C", color = "white", alpha = 0.6) +
  geom_density(fill = "#E74C3C", alpha = 0.3, linewidth = 1) +
  labs(title = "Scenario Uncertainty Distribution", subtitle = "Normalized Entropy (0-1)",
       x = "Entropy", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10), hjust = 0.5)
  )

ggsave(file.path(base_output_dir, "Entropy_distribution.tiff"), plot = p_entropy,
       width = 15, height = 10, units = "in", dpi = 600, device = "tiff", compression = "lzw")


cat("\nTüm analizler tamamlandı. Lütfen 'analiz_sonuclari' klasörünü kontrol edin.\n")
