Sys.setlocale("LC_ALL", "tr_TR.UTF-8") # Türkçe karakter desteği için


temp_dir <- "C:/R_Temp" 
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(TMPDIR = temp_dir)

# ================================================================
# --- PACKAGES ---
# ================================================================

required_packages <- c("psych", "GPArotation", "MASS", "readxl", "furrr", "purrr", "openxlsx", "plyr", "corpcor", "data.table", "future")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

library(future)

# ================================================================
# --- PARALLEL PROCESSING SETUP ---
# ================================================================

# Reset previous parallel plan
plan(sequential) 

# Start parallel build (using workers equal to the number of CPU cores minus one)
plan(multisession, workers = parallel::detectCores() - 1)

# ================================================================
# --- FIXED PARAMETER DEFINITIONS ---
# ================================================================

num_iterations = 10000    

print_individual_loadings = FALSE

set.seed(20250821)

# ================================================================
# --- VARIABLE AND SAMPLE NUMBER SETTINGS ---
# ================================================================

ds_values <- c(5, 10, 15, 20, 30, 40, 50) 

os_multipliers <- c(2, 3, 5, 10, 15, 20, 30, 40, 50) 

kk_values <- c(0.2, 0.5, 0.8) 

scenario_types <- c("sabit") 

all_scenario_freq_dt_list <- list()
all_scenario_winner_dt_list <- list()

# ================================================================
### AUXILIARY FUNCTIONS ###
# ================================================================

pad_loadings <- function(loadings, expected_cols, expected_rows) {
  if (!is.null(loadings)) {
    loadings <- as.matrix(unclass(loadings))
  }
  if (is.null(loadings) || !is.matrix(loadings) || nrow(loadings) == 0) {
    return(matrix(NA_real_, nrow = expected_rows, ncol = expected_cols))
  }
  
  current_rows <- nrow(loadings)
  current_cols <- ncol(loadings)
  
  if (current_rows == 0 || current_cols == 0) {
    return(matrix(NA_real_, nrow = expected_rows, ncol = expected_cols))
  }
  
  if (current_cols < expected_cols || current_rows < expected_rows) { 
    padded_mat <- matrix(NA_real_, nrow = expected_rows, ncol = expected_cols)
    padded_mat[1:min(current_rows, expected_rows), 1:min(current_cols, expected_cols)] <-
      loadings[1:min(current_rows, expected_rows), 1:min(current_cols, expected_cols)]
    return(padded_mat)
  } else if (current_cols > expected_cols || current_rows > expected_rows) {
    return(loadings[1:expected_rows, 1:expected_cols])
  } else {
    return(loadings)
  }
}

# ================================================================

safe_which_max_colname <- function(matrix_row, col_names) {
  if (is.null(matrix_row) || all(is.na(matrix_row)) || length(matrix_row) == 0) {
    return(NA_character_)
  }
  max_val <- max(matrix_row, na.rm = TRUE)
  if (is.infinite(max_val) || is.na(max_val)) {
    return(NA_character_)
  }
  
  max_indices <- which(matrix_row == max_val)
  
  if (length(max_indices) == 0) {
    return(NA_character_)
  }
  
  tied_names <- unique(gsub("_[0-9]+$", "", col_names[max_indices]))
  if (length(tied_names) > 1) {
    tied_names <- sort(tied_names) 
  }
  return(paste(tied_names, collapse = "-"))
}

# ================================================================

safe_which_min_colname <- function(matrix_row, col_names) {
  if (is.null(matrix_row) || all(is.na(matrix_row)) || length(matrix_row) == 0) {
    return(NA_character_)
  }
  min_val <- min(matrix_row, na.rm = TRUE)
  if (is.infinite(min_val) || is.na(min_val)) {
    return(NA_character_)
  }
  
  min_indices <- which(matrix_row == min_val)
  
  if (length(min_indices) == 0) {
    return(NA_character_)
  }
  
  tied_names <- unique(gsub("_[0-9]+$", "", col_names[min_indices]))
  if (length(tied_names) > 1) {
    tied_names <- sort(tied_names) 
  }
  return(paste(tied_names, collapse = "-"))
}

# ================================================================

calculate_ssi <- function(loadings_matrix, threshold = 0.3) {
  if (is.null(loadings_matrix) || !is.matrix(loadings_matrix) || 
      nrow(loadings_matrix) == 0 || ncol(loadings_matrix) == 0 || 
      all(!is.finite(loadings_matrix))) {
    return(NA)
  }
  abs_loadings <- abs(loadings_matrix)
  num_high_loadings_per_var <- rowSums(abs_loadings >= threshold & is.finite(abs_loadings), na.rm = TRUE)
  simple_variables_count <- sum(num_high_loadings_per_var == 1, na.rm = TRUE)
  ssi <- simple_variables_count / nrow(loadings_matrix)
  return(ssi)
}

# ================================================================

calculate_interpretability_metrics <- function(loadings_matrix) {
  if (is.null(loadings_matrix) || !is.matrix(loadings_matrix) || 
      nrow(loadings_matrix) == 0 || ncol(loadings_matrix) == 0 || 
      all(!is.finite(loadings_matrix))) {
    return(list(avg_primary_loading = NA_real_, avg_cross_loading = NA_real_, primary_to_cross_ratio = NA_real_))
  }
  
  abs_loadings <- abs(loadings_matrix)
  
  primary_loadings <- numeric(nrow(abs_loadings))
  cross_loadings_values <- c() 
  
  for (i in 1:nrow(abs_loadings)) {
    row_loadings <- abs_loadings[i, ]
    finite_row_loadings <- row_loadings[is.finite(row_loadings)]
    
    if (length(finite_row_loadings) == 0) {
      primary_loadings[i] <- NA_real_
      next
    }
    
    max_val <- max(finite_row_loadings, na.rm = TRUE)
    if (is.infinite(max_val) || is.na(max_val)) {
      primary_loadings[i] <- NA_real_
      next
    }
    
    max_idx_orig <- which(row_loadings == max_val & is.finite(row_loadings))[1]
    
    if (length(max_idx_orig) == 0 || is.na(max_idx_orig)) {
      primary_loadings[i] <- NA_real_
      next
    }
    
    primary_loadings[i] <- row_loadings[max_idx_orig]
    
    cross_loadings_i <- row_loadings[-max_idx_orig]
    cross_loadings_i <- cross_loadings_i[is.finite(cross_loadings_i) & !is.na(cross_loadings_i)]
    if (length(cross_loadings_i) > 0) {
      cross_loadings_values <- c(cross_loadings_values, cross_loadings_i)
    }
  }
  
  avg_primary <- mean(primary_loadings, na.rm = TRUE)
  
  if (length(cross_loadings_values) > 0) {
    avg_cross <- mean(cross_loadings_values, na.rm = TRUE)
  } else {
    avg_cross <- 0
  }
  
  ratio <- if (!is.na(avg_cross) && avg_cross > 0) {
    avg_primary / avg_cross
  } else if (!is.na(avg_primary) && avg_primary > 0 && (is.na(avg_cross) || avg_cross == 0)) {
    Inf
  } else {
    NA_real_
  }
  
  return(list(avg_primary_loading = avg_primary, avg_cross_loading = avg_cross, primary_to_cross_ratio = ratio))
}

# ================================================================

safe_max <- function(x, na.rm = FALSE) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = na.rm)
}

# ================================================================

safe_min <- function(x, na.rm = FALSE) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = na.rm)
}

# ================================================================

format_frequencies_wide <- function(freq_table) { 
  filtered_freq_table <- freq_table[freq_table > 0 & !is.na(freq_table)]
  sorted_freq_table <- sort(filtered_freq_table, decreasing = TRUE)
  
  result_list <- list()
  current_pair_idx <- 1
  
  grouped_frequencies <- list()
  for(name in names(sorted_freq_table)) {
    freq <- as.numeric(sorted_freq_table[name])
    if (as.character(freq) %in% names(grouped_frequencies)) {
      grouped_frequencies[[as.character(freq)]] <- c(grouped_frequencies[[as.character(freq)]], gsub("_[0-9]+$", "", name))
    } else {
      grouped_frequencies[[as.character(freq)]] <- c(gsub("_[0-9]+$", "", name))
    }
  }
  
  sorted_grouped_frequencies <- grouped_frequencies[order(as.numeric(names(grouped_frequencies)), decreasing = TRUE)]
  
  for (freq_val_str in names(sorted_grouped_frequencies)) {
    methods <- sorted_grouped_frequencies[[freq_val_str]]
    methods <- sort(unique(methods)) 
    combined_method_name <- paste(methods, collapse = "-")
    
    name_suffix <- ifelse(current_pair_idx == 1, "", paste0("_", current_pair_idx))
    
    result_list[[paste0("Rotation", name_suffix)]] <- combined_method_name
    result_list[[paste0("Frequency", name_suffix)]] <- as.numeric(freq_val_str) 
    
    current_pair_idx <- current_pair_idx + 1
  }
  
  return(result_list)
}

# ================================================================

process_scenario_freq_table <- function(freq_table, measure_type, current_ds, current_np, current_kk, current_kk_tip, ssi_df, interpret_df, conv_df) {
  df <- data.frame(method = names(freq_table), freq = as.numeric(freq_table), stringsAsFactors = FALSE)
  df$ds <- current_ds
  df$np <- current_np
  df$kk <- current_kk
  df$kk_tip <- current_kk_tip
  df$measure <- measure_type
  
  conv_df$Rotation <- as.character(conv_df$Rotation)
  
  ssi_for_merge <- ssi_df[, c("Rotation", "Mean_SSI")]
  colnames(ssi_for_merge)[colnames(ssi_for_merge) == "Rotation"] <- "method"
  
  interpret_for_merge <- interpret_df[, c("Rotation", "Mean_Primary_Loading", "Mean_Cross_Loading")]
  colnames(interpret_for_merge)[colnames(interpret_for_merge) == "Rotation"] <- "method"
  
  conv_for_merge <- conv_df[, c("Rotation", "Convergence_Ratio")]
  colnames(conv_for_merge)[colnames(conv_for_merge) == "Rotation"] <- "method"
  
  df_merged <- merge(df, ssi_for_merge, by = "method", all.x = TRUE)
  df_merged <- merge(df_merged, interpret_for_merge, by = "method", all.x = TRUE)
  df_merged <- merge(df_merged, conv_for_merge, by = "method", all.x = TRUE)
  
  colnames(df_merged)[colnames(df_merged) == "Mean_SSI"] <- "mean_SSI"
  colnames(df_merged)[colnames(df_merged) == "Mean_Primary_Loading"] <- "mean_primary"
  colnames(df_merged)[colnames(df_merged) == "Mean_Cross_Loading"] <- "mean_cross"
  colnames(df_merged)[colnames(df_merged) == "Convergence_Ratio"] <- "converged_ratio"
  
  final_cols <- c("ds", "np", "kk", "kk_tip", "measure", "method", "freq", "converged_ratio", "mean_SSI", "mean_primary", "mean_cross")
  df_merged <- df_merged[, final_cols]
  
  return(df_merged)
}

# ================================================================

count_wins <- function(winner_vec, methods) {
  counts <- setNames(numeric(length(methods)), methods)
  for (w in winner_vec) {
    if (is.na(w)) next
    parts <- strsplit(w, "-", fixed = TRUE)[[1]]
    weight <- 1 / length(parts)
    for (p in parts) {
      if (p %in% methods) {
        counts[p] <- counts[p] + weight
      }
    }
  }
  counts
}

# ================================================================

run_one_iteration <- function(iter_idx, os, ds, n_factors, sample_meanvector, sample_covariance_matrix) {
  
  veri = mvrnorm(
    n = os,
    mu = sample_meanvector,
    Sigma = sample_covariance_matrix,
    empirical = FALSE
  )
  
  convergence_status <- list(
    varimax = c(converged = 1, non_converged = 0),
    quartimax = c(converged = 1, non_converged = 0),
    equamax = c(converged = 1, non_converged = 0),
    cf_varimax = c(converged = 1, non_converged = 0),
    geominT = c(converged = 1, non_converged = 0),
    entropy = c(converged = 1, non_converged = 0),
    infomaxT = c(converged = 1, non_converged = 0),
    bifactorT = c(converged = 1, non_converged = 0),
    bentlerT = c(converged = 1, non_converged = 0),
    varimin = c(converged = 1, non_converged = 0)
  )
  
# ================================================================
 
  run_rotation <- function(rotation_func, data, n_factors_val, rotation_type, unrotated_loadings_mat = NULL) {
    loadings <- NULL
    converged_status <- 0 
    non_converged_status <- 1
    
    tryCatch({
      if (rotation_type %in% c("varimax", "quartimax", "equamax", "varimin")) {
        res <- principal(data, nfactor = n_factors_val, rotate = rotation_type, covar = TRUE, scores = FALSE, cor = "cov")
        loadings <- res$loadings
      } else {
        if (is.null(unrotated_loadings_mat)) {
         stop("GPArotation fonksiyonları için döndürülmemiş yük matrisi sağlanmalıdır.")
        }
        
        if (rotation_type == "cfT") {
          res <- rotation_func(unrotated_loadings_mat, kappa = 1 / n_factors_val)
        } else {
          res <- rotation_func(unrotated_loadings_mat)
        }
        loadings <- res$loadings
      }
      
      if (!is.null(loadings)) {
        loadings <- as.matrix(unclass(loadings))
        if (all(is.finite(loadings))) {
          converged_status <- 1
          non_converged_status <- 0
        } else {
          loadings <- NULL 
        }
      } else {
        loadings <- NULL 
      }
      
    }, warning = function(w) {
      warn_msg <- conditionMessage(w)
      if (grepl("did not converge", warn_msg) || 
          grepl("Matrix was not positive definite", warn_msg) || 
          grepl("NaNs produced", warn_msg)) {
        loadings <- NULL 
        converged_status <- 0
        non_converged_status <- 1
        invokeRestart("muffleWarning") 
      } else {
        if (!is.null(loadings)) {
          loadings <- as.matrix(unclass(loadings))
          if (all(is.finite(loadings))) {
            converged_status <- 1
            non_converged_status <- 0
          } else {
            loadings <- NULL
            converged_status <- 0
            non_converged_status <- 1
          }
        } else {
          loadings <- NULL
          converged_status <- 0
          non_converged_status <- 1
        }
        invokeRestart("muffleWarning")
      }
    }, error = function(e) {
      loadings <- NULL 
      converged_status <- 0
      non_converged_status <- 1
    })
    
    return(list(loadings = loadings, converged = converged_status, non_converged = non_converged_status))
  }
  
  none_res_raw <- tryCatch(principal(veri, nfactor=n_factors, rotate="none", covar=TRUE, scores=FALSE, cor="cov"),
                           error = function(e) list(loadings = NULL))
  
  none_loadings_valid <- if(!is.null(none_res_raw$loadings)) {
    as.matrix(unclass(none_res_raw$loadings))
  } else {
    NULL
  }
  none_yuk_i <- abs(pad_loadings(none_loadings_valid, n_factors, ds))
  
 rotation_methods <- list(
    varimax = list(func = principal, type = "varimax"),
    quartimax = list(func = principal, type = "quartimax"),
    equamax = list(func = principal, type = "equamax"),
    cf_varimax = list(func = GPArotation::cfT, type = "cfT"),
    geominT = list(func = GPArotation::geominT, type = "geominT"),
    entropy = list(func = GPArotation::entropy, type = "entropy"),
    infomaxT = list(func = GPArotation::infomaxT, type = "infomaxT"),
    bifactorT = list(func = GPArotation::bifactorT, type = "bifactorT"),
    bentlerT = list(func = GPArotation::bentlerT, type = "bentlerT"),
    varimin = list(func = principal, type = "varimin")
  )
  
  rotated_results <- list()
  
  for (method_name in names(rotation_methods)) {
    method_info <- rotation_methods[[method_name]]
    
    unrotated_matrix_to_pass <- NULL
    if (!(method_info$type %in% c("varimax", "quartimax", "equamax", "varimin"))) {
      unrotated_matrix_to_pass <- none_loadings_valid 
    }
    
    res <- run_rotation(
      rotation_func = method_info$func, 
      data = veri, 
      n_factors_val = n_factors, 
      rotation_type = method_info$type, 
      unrotated_loadings_mat = unrotated_matrix_to_pass
    )
    
    rotated_results[[paste0(method_name, "_yuk")]] <- abs(pad_loadings(res$loadings, n_factors, ds))
    convergence_status[[method_name]] <- c(converged = res$converged, non_converged = res$non_converged)
  }
  
  all_ssi <- list(none_ssi = calculate_ssi(none_yuk_i))
  all_interpret_metrics <- list(none = calculate_interpretability_metrics(none_yuk_i))
  
  for (method_name in names(rotation_methods)) {
    loadings_matrix <- rotated_results[[paste0(method_name, "_yuk")]]
    all_ssi[[paste0(method_name, "_ssi")]] <- calculate_ssi(loadings_matrix)
    all_interpret_metrics[[method_name]] <- calculate_interpretability_metrics(loadings_matrix)
  }
  
  result_list <- list(
    none_yuk = none_yuk_i,
    none_ssi = all_ssi$none_ssi,
    none_avg_primary_loading = all_interpret_metrics$none$avg_primary_loading,
    none_avg_cross_loading = all_interpret_metrics$none$avg_cross_loading,
    none_primary_to_cross_ratio = all_interpret_metrics$none$primary_to_cross_ratio
  )
  
  for (method_name in names(rotation_methods)) {
    result_list[[paste0(method_name, "_yuk")]] <- rotated_results[[paste0(method_name, "_yuk")]]
    result_list[[paste0(method_name, "_ssi")]] <- all_ssi[[paste0(method_name, "_ssi")]]
    result_list[[paste0(method_name, "_avg_primary_loading")]] <- all_interpret_metrics[[method_name]]$avg_primary_loading
    result_list[[paste0(method_name, "_avg_cross_loading")]] <- all_interpret_metrics[[method_name]]$avg_cross_loading
    result_list[[paste0(method_name, "_primary_to_cross_ratio")]] <- all_interpret_metrics[[method_name]]$primary_to_cross_ratio
    result_list[[paste0(method_name, "_converged")]] <- convergence_status[[method_name]]["converged"]
    result_list[[paste0(method_name, "_non_converged")]] <- convergence_status[[method_name]]["non_converged"]
  }
  
  none_toplam_i <- none_yuk_i 
  dik_toplam_i <- do.call(cbind, lapply(names(rotation_methods), function(name) rotated_results[[paste0(name, "_yuk")]]))
  colnames(dik_toplam_i) <- unlist(lapply(names(rotation_methods), function(name) rep(name, n_factors))) 
  toplam_i <- cbind(none_yuk_i, dik_toplam_i)
  colnames(toplam_i) <- c(rep("none", n_factors), unlist(lapply(names(rotation_methods), function(name) rep(name, n_factors))))
  
  result_list$none_maks <- apply(none_toplam_i, 1, safe_which_max_colname, col_names = colnames(none_toplam_i))
  result_list$none_min <- apply(none_toplam_i, 1, safe_which_min_colname, col_names = colnames(none_toplam_i))
  
  result_list$dik_maks <- apply(dik_toplam_i, 1, safe_which_max_colname, col_names = colnames(dik_toplam_i))
  result_list$dik_min <- apply(dik_toplam_i, 1, safe_which_min_colname, col_names = colnames(dik_toplam_i))
  
  result_list$toplam_maks <- apply(toplam_i, 1, safe_which_max_colname, col_names = colnames(toplam_i))
  result_list$toplam_min <- apply(toplam_i, 1, safe_which_min_colname, col_names = colnames(toplam_i))
  
  result_list$none_maks_deger <- apply(none_yuk_i, 1, safe_max, na.rm = TRUE)
  result_list$none_min_deger <- apply(none_yuk_i, 1, safe_min, na.rm = TRUE)
  
  for (method_name in names(rotation_methods)) {
    yuk_matrix <- rotated_results[[paste0(method_name, "_yuk")]]
    result_list[[paste0(method_name, "_maks_deger")]] <- apply(yuk_matrix, 1, safe_max, na.rm = TRUE)
    result_list[[paste0(method_name, "_min_deger")]] <- apply(yuk_matrix, 1, safe_min, na.rm = TRUE)
  }
  
  fark_list <- list()
  fark_list$none <- result_list$none_maks_deger - result_list$none_min_deger
  for (method_name in names(rotation_methods)) {
    fark_list[[method_name]] <- result_list[[paste0(method_name, "_maks_deger")]] - result_list[[paste0(method_name, "_min_deger")]]
  }
  
  result_list$fark_none <- fark_list$none
  
  fark_toplam_i <- do.call(cbind, lapply(fark_list, function(x) matrix(x, ncol = 1))) 
  colnames(fark_toplam_i) <- c("none", names(rotation_methods)) 
  result_list$fark_toplam <- fark_toplam_i
  
  result_list$fark_maks <- apply(fark_toplam_i, 1, safe_which_max_colname, col_names = colnames(fark_toplam_i))
  result_list$fark_min <- apply(fark_toplam_i, 1, safe_which_min_colname, col_names = colnames(fark_toplam_i))
  
  return(result_list)
}

rotation_order <- c("none", "varimax", "quartimax", "equamax", "cf_varimax",
                    "geominT", "entropy", "infomaxT", "bifactorT", "bentlerT", "varimin")

# ================================================================
# --- MAIN SIMULATION CYCLE ---
# ================================================================

for (scenario_type in scenario_types) {
  
  if (scenario_type == "dengesiz_tip2") {
    kk_to_iterate <- kk_values[1] 
  } else {
    kk_to_iterate <- kk_values 
  }
  
  for (current_kk in kk_to_iterate) { 
    kk <- current_kk
    
    for (current_ds in ds_values) {
      ds <- current_ds
      
      
      if (scenario_type == "dengesiz_tip2") {
        kk_tip_label <- paste0(scenario_type, "_fixedpattern") 
      } else {
        kk_tip_label <- paste0(scenario_type, "_kk", gsub("\\.", "p", kk))
      }
      
      current_output_base_dir <- file.path(getwd(), "simülasyon_sonuçları", 
                                           paste0("scenario_", kk_tip_label, "_ds", ds))
      if (!dir.exists(current_output_base_dir)) {
        dir.create(current_output_base_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      sample_meanvector <- rep(0, ds) 
      original_sample_covariance_matrix <- diag(ds) 
      
      if (scenario_type == "sabit") {
        original_sample_covariance_matrix[!diag(ds)] <- kk
      } else if (scenario_type == "dengesiz_tip1") {
        for (i in 1:ds) {
          for (j in 1:ds) {
            if (i != j) {
              original_sample_covariance_matrix[i, j] <- kk * ((-1)^(abs(i-j) + 1))
            }
          }
        }
      } else if (scenario_type == "dengesiz_tip2") {
        num_elements <- ds * (ds - 1) / 2
        if (num_elements > 0) {
          repeating_pattern <- c(0.2, 0.4, 0.6, 0.8)
          original_sample_covariance_matrix[lower.tri(original_sample_covariance_matrix)] <- rep(repeating_pattern, length.out = num_elements)
          original_sample_covariance_matrix[upper.tri(original_sample_covariance_matrix)] <- t(original_sample_covariance_matrix)[upper.tri(original_sample_covariance_matrix)]
        }
      }
      
      sample_covariance_matrix_to_use <- original_sample_covariance_matrix 
      
      eigen_values <- eigen(original_sample_covariance_matrix, symmetric = TRUE)$values
      is_positive_definite <- !any(eigen_values <= .Machine$double.eps * max(abs(eigen_values)))
      
      matrix_output_dir <- file.path(current_output_base_dir, "kovaryans_matrisleri")
      if (!dir.exists(matrix_output_dir)) {
        dir.create(matrix_output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      wb_matrix <- createWorkbook()
      addWorksheet(wb_matrix, "Original_Covariance")
      writeData(wb_matrix, "Original_Covariance", original_sample_covariance_matrix, rowNames = TRUE, colNames = TRUE)
      
      if (!is_positive_definite) {
        warning(paste0("Oluşturulan kovaryans matrisi pozitif tanımlı değil veya neredeyse pozitif tanımlı değil. Düzeltme yapılıyor. (Senaryo=", scenario_type, ", ds=", ds, ", kk=", kk, ")"))
        sample_covariance_matrix_to_use <- corpcor::make.positive.definite(original_sample_covariance_matrix)
        addWorksheet(wb_matrix, "Corrected_Covariance")
        writeData(wb_matrix, "Corrected_Covariance", sample_covariance_matrix_to_use, rowNames = TRUE, colNames = TRUE)
      }
      
      saveWorkbook(wb_matrix, 
                   file.path(matrix_output_dir, paste0("ds", ds, "_kk", gsub("\\.", "p", kk), "_", scenario_type, "_covariance_matrices.xlsx")), 
                   overwrite = TRUE)
      
      for (multiplier in os_multipliers) {
        os <- current_ds * multiplier
        n_factors <- ds 
        
        np_ratio_label <- paste0(os / ds) 
        
        cat(paste0("\nSimülasyon Başlatılıyor: Senaryo = ", kk_tip_label, 
                   ", ds = ", ds, ", os = ", os, ", kk = ", kk, ", ts = ", num_iterations, "\n"))
        
        all_results_for_current_os <- future_map(seq_len(num_iterations), ~run_one_iteration(
          .x, os, ds, n_factors, sample_meanvector, sample_covariance_matrix_to_use
        ), .options = furrr_options(seed = TRUE))
        
        if (print_individual_loadings && num_iterations <= 10) {
          cat("\n--- İlk Birkaç İterasyon İçin Yük Matrisleri (Mutlak Değerler) ---\n")
          for (i in 1:min(num_iterations, 5)) {
            cat(paste0("\n--- İterasyon ", i, " ---\n"))
            
            cat("\nNone Yük Matrisi:\n")
            print(all_results_for_current_os[[i]]$none_yuk)
            
            rotation_methods_print <- c("varimax", "quartimax", "equamax", "cf_varimax",
                                        "geominT", "entropy", "infomaxT", "bifactorT", "bentlerT", "varimin")
            for (method_name in rotation_methods_print) {
              cat(paste0("\n", method_name, " Yük Matrisi:\n"))
              print(all_results_for_current_os[[i]][[paste0(method_name, "_yuk")]])
            }
          }
          cat("\n------------------------------------------------\n")
        }

# ================================================================
# --- RESULT SUMMARY AND OUTPUT ---
# ================================================================
        
        ssi_results_df <- data.frame(
          Rotation = factor(c("none", "varimax", "quartimax", "equamax", "cf_varimax",
                              "geominT", "entropy", "infomaxT", "bifactorT", "bentlerT", "varimin"), levels = rotation_order),
          np = np_ratio_label,
          Mean_SSI = c(
            round(mean(unlist(purrr::map(all_results_for_current_os, "none_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimax_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "quartimax_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "equamax_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "cf_varimax_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "geominT_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "entropy_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "infomaxT_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bifactorT_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bentlerT_ssi")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimin_ssi")), na.rm = TRUE), 4)
          ),
          SD_SSI = c(
            round(sd(unlist(purrr::map(all_results_for_current_os, "none_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimax_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "quartimax_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "equamax_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "cf_varimax_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "geominT_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "entropy_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "infomaxT_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bifactorT_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bentlerT_ssi")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimin_ssi")), na.rm = TRUE), 4)
          )
        )
        
        interpretability_metrics_summary_df <- data.frame(
          Rotation = factor(c("none", "varimax", "quartimax", "equamax", "cf_varimax",
                              "geominT", "entropy", "infomaxT", "bifactorT", "bentlerT", "varimin"), levels = rotation_order),
          np = np_ratio_label,
          Mean_Primary_Loading = c(
            round(mean(unlist(purrr::map(all_results_for_current_os, "none_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimax_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "quartimax_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "equamax_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "cf_varimax_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "geominT_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "entropy_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "infomaxT_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bifactorT_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bentlerT_avg_primary_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimin_avg_primary_loading")), na.rm = TRUE), 4)
          ),
          SD_Primary_Loading = c(
            round(sd(unlist(purrr::map(all_results_for_current_os, "none_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimax_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "quartimax_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "equamax_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "cf_varimax_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "geominT_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "entropy_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "infomaxT_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bifactorT_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bentlerT_avg_primary_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimin_avg_primary_loading")), na.rm = TRUE), 4)
          ),
          Mean_Cross_Loading = c(
            round(mean(unlist(purrr::map(all_results_for_current_os, "none_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimax_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "quartimax_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "equamax_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "cf_varimax_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "geominT_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "entropy_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "infomaxT_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bifactorT_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bentlerT_avg_cross_loading")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimin_avg_cross_loading")), na.rm = TRUE), 4)
          ),
          SD_Cross_Loading = c(
            round(sd(unlist(purrr::map(all_results_for_current_os, "none_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimax_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "quartimax_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "equamax_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "cf_varimax_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "geominT_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "entropy_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "infomaxT_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bifactorT_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bentlerT_avg_cross_loading")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimin_avg_cross_loading")), na.rm = TRUE), 4)
          ),
          Primary_to_Cross_Ratio = c(
            round(mean(unlist(purrr::map(all_results_for_current_os, "none_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "quartimax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "equamax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "cf_varimax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "geominT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "entropy_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "infomaxT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bifactorT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "bentlerT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(mean(unlist(purrr::map(all_results_for_current_os, "varimin_primary_to_cross_ratio")), na.rm = TRUE), 4)
          ),
          SD_Primary_to_Cross_Ratio = c(
            round(sd(unlist(purrr::map(all_results_for_current_os, "none_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "quartimax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "equamax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "cf_varimax_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "geominT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "entropy_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "infomaxT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bifactorT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "bentlerT_primary_to_cross_ratio")), na.rm = TRUE), 4),
            round(sd(unlist(purrr::map(all_results_for_current_os, "varimin_primary_to_cross_ratio")), na.rm = TRUE), 4)
          )
        )
        
        convergence_summary_df <- data.frame(
          Rotation = factor(c("varimax", "quartimax", "equamax", "cf_varimax",
                              "geominT", "entropy", "infomaxT", "bifactorT", "bentlerT", "varimin"), levels = rotation_order[-1]),
          np = np_ratio_label,
          Converged_Count = c(
            sum(unlist(purrr::map(all_results_for_current_os, "varimax_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "quartimax_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "equamax_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "cf_varimax_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "geominT_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "entropy_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "infomaxT_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "bifactorT_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "bentlerT_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "varimin_converged")), na.rm = TRUE)
          ),
          Non_Converged_Count = c(
            sum(unlist(purrr::map(all_results_for_current_os, "varimax_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "quartimax_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "equamax_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "cf_varimax_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "geominT_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "entropy_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "infomaxT_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "bifactorT_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "bentlerT_non_converged")), na.rm = TRUE),
            sum(unlist(purrr::map(all_results_for_current_os, "varimin_non_converged")), na.rm = TRUE)
          ),
          Convergence_Ratio = c(
            round(sum(unlist(purrr::map(all_results_for_current_os, "varimax_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "quartimax_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "equamax_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "cf_varimax_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "geominT_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "entropy_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "infomaxT_converged")), na.rm = TRUE) / num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "bifactorT_converged")), na.rm = TRUE)/ num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "bentlerT_converged")), na.rm = TRUE)/ num_iterations, 4),
            round(sum(unlist(purrr::map(all_results_for_current_os, "varimin_converged")), na.rm = TRUE) / num_iterations, 4)
          )
        )
        
        temp_summary_rows <- lapply(seq_len(ds), function(j) {
          var_name <- paste0("x", sprintf("%02d", j))
          
          temp_toplam_maks_vec <- unlist(purrr::map(all_results_for_current_os, ~ {
            val <- .x$toplam_maks[j]
            if (is.null(val) || is.na(val)) character(0) else as.character(val)
          }))
          
          temp_toplam_min_vec <- unlist(purrr::map(all_results_for_current_os, ~ {
            val <- .x$toplam_min[j]
            if (is.null(val) || is.na(val)) character(0) else as.character(val)
          }))
          
          temp_fark_maks_vec <- unlist(purrr::map(all_results_for_current_os, ~ {
            val <- .x$fark_maks[j]
            if (is.null(val) || is.na(val)) character(0) else as.character(val)
          }))
          
          max_freq_counts <- count_wins(temp_toplam_maks_vec, rotation_order)
          min_freq_counts <- count_wins(temp_toplam_min_vec, rotation_order)
          range_max_freq_counts <- count_wins(temp_fark_maks_vec, rotation_order)
          
          max_freq_table <- max_freq_counts[max_freq_counts > 0] 
          min_freq_table <- min_freq_counts[min_freq_counts > 0]
          range_max_freq_table <- range_max_freq_counts[range_max_freq_counts > 0]
          
          max_freq_row_data_list <- format_frequencies_wide(max_freq_table)
          min_freq_row_data_list <- format_frequencies_wide(min_freq_table)
          range_max_freq_row_data_list <- format_frequencies_wide(range_max_freq_table)
          
          list(
            data.frame(
              Variable = var_name,
              np = np_ratio_label,
              min_max = "min",
              min_freq_row_data_list,
              stringsAsFactors = FALSE,
              check.names = FALSE
            ),
            data.frame(
              Variable = var_name,
              np = np_ratio_label,
              min_max = "max",
              max_freq_row_data_list,
              stringsAsFactors = FALSE,
              check.names = FALSE
            ),
            data.frame(
              Variable = var_name,
              np = np_ratio_label,
              min_max = "range_max",
              range_max_freq_row_data_list,
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          )
        })
        
        final_summary_df <- do.call(plyr::rbind.fill, unlist(temp_summary_rows, recursive = FALSE))
        
        write.csv2(final_summary_df, 
                   file = file.path(current_output_base_dir, paste0("ds", ds, "_os", os, "_summary_rotations.csv")), 
                   row.names = FALSE)
        
        combined_metrics_df <- merge(ssi_results_df, interpretability_metrics_summary_df, by = c("Rotation", "np"), all = TRUE)
        combined_metrics_df <- merge(combined_metrics_df, convergence_summary_df, by = c("Rotation", "np"), all = TRUE)
        
        combined_metrics_df$Rotation <- factor(as.character(combined_metrics_df$Rotation), levels = rotation_order)
        combined_metrics_df <- combined_metrics_df[order(combined_metrics_df$Rotation, combined_metrics_df$np), ]
        
        all_metrics_wb <- createWorkbook()
        addWorksheet(all_metrics_wb, "Summary Metrics")
        writeData(all_metrics_wb, "Summary Metrics", combined_metrics_df)
        
        saveWorkbook(all_metrics_wb, 
                     file.path(current_output_base_dir, paste0("ds", ds, "_os", os, "_all_summary_metrics.xlsx")), 
                     overwrite = TRUE)
        
        cat(paste0("Simülasyon Tamamlandı: Senaryo = ", kk_tip_label, 
                   ", ds = ", ds, ", os = ", os, "\n"))
        
        all_max_winners_across_vars <- unlist(purrr::map(all_results_for_current_os, ~ .x$toplam_maks))
        all_min_winners_across_vars <- unlist(purrr::map(all_results_for_current_os, ~ .x$toplam_min))
        all_range_max_winners_across_vars <- unlist(purrr::map(all_results_for_current_os, ~ .x$fark_maks))
        
        total_max_freq_across_vars <- count_wins(all_max_winners_across_vars, rotation_order)
        total_min_freq_across_vars <- count_wins(all_min_winners_across_vars, rotation_order)
        total_range_max_freq_across_vars <- count_wins(all_range_max_winners_across_vars, rotation_order)
        
        df_max_freq <- process_scenario_freq_table(total_max_freq_across_vars, "max", ds, as.numeric(np_ratio_label), kk, kk_tip_label, ssi_results_df, interpretability_metrics_summary_df, convergence_summary_df)
        df_min_freq <- process_scenario_freq_table(total_min_freq_across_vars, "min", ds, as.numeric(np_ratio_label), kk, kk_tip_label, ssi_results_df, interpretability_metrics_summary_df, convergence_summary_df)
        df_range_max_freq <- process_scenario_freq_table(total_range_max_freq_across_vars, "range_max", ds, as.numeric(np_ratio_label), kk, kk_tip_label, ssi_results_df, interpretability_metrics_summary_df, convergence_summary_df)
        
        all_scenario_freq_dt_list[[length(all_scenario_freq_dt_list) + 1]] <- as.data.table(df_max_freq)
        all_scenario_freq_dt_list[[length(all_scenario_freq_dt_list) + 1]] <- as.data.table(df_min_freq)
        all_scenario_freq_dt_list[[length(all_scenario_freq_dt_list) + 1]] <- as.data.table(df_range_max_freq)
        
        temp_scenario_winner_dt_rows <- data.table(
          ds = character(0), np = character(0), kk = numeric(0), kk_tip = character(0), Variable = character(0), measure = character(0), winner = character(0)
        ) 
        
        for (iter_idx_inner in seq_len(num_iterations)) {
          for (j in seq_len(ds)) {
            var_name <- paste0("x", sprintf("%02d", j)) 
            temp_dt_to_add <- data.table(
              ds = ds,
              np = np_ratio_label,
              kk = kk,
              kk_tip = kk_tip_label,
              Variable = var_name, 
              measure = c("max", "min", "range_max"),
              winner = c(
                all_results_for_current_os[[iter_idx_inner]]$toplam_maks[j],
                all_results_for_current_os[[iter_idx_inner]]$toplam_min[j],
                all_results_for_current_os[[iter_idx_inner]]$fark_maks[j]
              )
            )
            temp_scenario_winner_dt_rows <- rbindlist(list(temp_scenario_winner_dt_rows, temp_dt_to_add), fill = TRUE)
          }
        }
        
        if (nrow(temp_scenario_winner_dt_rows) > 0) {
          all_scenario_winner_dt_list[[length(all_scenario_winner_dt_list) + 1]] <- temp_scenario_winner_dt_rows
        }
      }
    }
  }
}

# ================================================================
# --- SAVING ALL SIMULATION RESULTS ---
# ================================================================

if (length(all_scenario_freq_dt_list) > 0) {
  final_scenario_freq_dt <- rbindlist(all_scenario_freq_dt_list, fill = TRUE) 
  saveRDS(final_scenario_freq_dt, file = file.path(getwd(), "simülasyon_sonuçları", "senaryo_frekans.rds")) 
  cat("\n'senaryo_frekans.rds' dosyası oluşturuldu ve 'simülasyon_sonuçları' dizinine kaydedildi.\n") 
} else {
  cat("\n'senaryo_frekans.rds' oluşturulamadı çünkü toplanacak veri yok.\n") 
}

if (length(all_scenario_winner_dt_list) > 0) {
  final_scenario_winner_dt <- rbindlist(all_scenario_winner_dt_list, fill = TRUE) 
  saveRDS(final_scenario_winner_dt, file = file.path(getwd(), "simülasyon_sonuçları", "senaryo_winner.rds"))
  cat("\n'senaryo_winner.rds' dosyası oluşturuldu ve 'simülasyon_sonuçları' dizinine kaydedildi.\n")
} else {
  cat("\n'senaryo_winner.rds' oluşturulamadı çünkü toplanacak veri yok.\n")
}

cat("\nTüm parametre kombinasyonları ve senaryolar için simülasyonlar tamamlandı. Sonuçlar 'simülasyon_sonuçları' dizininde bulunmaktadır.\n")
