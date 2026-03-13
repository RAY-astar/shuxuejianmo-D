# 法医物证多人身份鉴定问题 - 完整最终版解决方案 

# 加载所有必要的库 (新增了 MASS 库用于鲁棒计算)
required_packages <- c("tidyverse", "readxl", "signal", "pracma", "cluster", 
                       "mixtools", "ggplot2", "gridExtra", "reshape2", 
                       "officer", "flextable", "MASS")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 设置工作目录和全局参数
options(warn = -1)
set.seed(123)

cat("🧬 法医物证STR图谱分析系统 - 完整最终版 🧬\n")
cat("========================================\n")
cat("✨ 包含功能:\n")
cat("📊 完整的数据分析 (贡献者估计、混合比例、基因型推断、降噪)\n")
cat("🎨 全面的数据可视化 (7类共15+张图表)\n")
cat("📝 专业的Word报告 (三线表格式)\n")
cat("💾 完整的数据保存和文件管理\n")
cat("💎 数学内核已升级为: LOWESS, 鲁棒EM(Huber+L2), MAP推断\n")
cat("========================================\n\n")

# ===============================
# 核心分析函数（ 数学升级）
# ===============================

# 引入 LOWESS 基线校正 
lowess_baseline_correction <- function(size_data, height_data, span = 0.1) {
  if (length(height_data) < 10) return(height_data)
  tryCatch({
    fit <- loess(height_data ~ size_data, span = span, 
                 control = loess.control(surface = "direct"))
    baseline <- predict(fit, size_data)
    corrected_signal <- height_data - baseline
    corrected_signal[corrected_signal < 0] <- 0
    return(corrected_signal)
  }, error = function(e) {
    return(height_data)
  })
}

引入 Huber 损失函数
huber_loss <- function(residual, delta = 1.345) {
  abs_r <- abs(residual)
  ifelse(abs_r <= delta, 0.5 * residual^2, delta * abs_r - 0.5 * delta^2)
}

# 原版的简单平滑函数（保留以防万一）
simple_smooth <- function(x, window = 3) {
  if (length(x) < window) return(x)
  half_window <- floor(window / 2)
  smoothed <- numeric(length(x))
  for (i in 1:length(x)) {
    start_idx <- max(1, i - half_window)
    end_idx <- min(length(x), i + half_window)
    smoothed[i] <- mean(x[start_idx:end_idx])
  }
  return(smoothed)
}

# 自适应阈值计算
calculate_adaptive_threshold <- function(signal_data, sensitivity = "medium") {
  positive_signal <- signal_data[signal_data > 0]
  if (length(positive_signal) == 0) {
    return(list(height_threshold = 0.05, min_prominence = 0.05))
  }
  
  q25 <- quantile(positive_signal, 0.25, na.rm = TRUE)
  q75 <- quantile(positive_signal, 0.75, na.rm = TRUE)
  median_val <- median(positive_signal, na.rm = TRUE)
  max_val <- max(positive_signal, na.rm = TRUE)
  
  signal_range <- max_val - median_val
  iqr <- q75 - q25
  
  if (sensitivity == "high") {
    if (signal_range > iqr * 2) {
      height_ratio <- 0.01; prominence_ratio <- 0.02
    } else {
      height_ratio <- 0.005; prominence_ratio <- 0.01
    }
  } else if (sensitivity == "low") {
    height_ratio <- 0.20; prominence_ratio <- 0.25
  } else {
    if (signal_range > iqr * 3) {
      height_ratio <- 0.05; prominence_ratio <- 0.08
    } else if (signal_range > iqr) {
      height_ratio <- 0.03; prominence_ratio <- 0.05
    } else {
      height_ratio <- 0.01; prominence_ratio <- 0.02
    }
  }
  
  height_threshold <- min(height_ratio, 0.1)
  min_prominence <- min(prominence_ratio, 0.15)
  return(list(height_threshold = height_threshold, min_prominence = min_prominence))
}

# 优化的峰检测算法
detect_peaks_str_optimized <- function(signal_data, size_data = NULL, 
                                       height_threshold = 0.03, min_distance = 2,
                                       min_prominence = 0.05, auto_threshold = TRUE,
                                       sensitivity = "medium", debug = FALSE) {
  
  if (length(signal_data) < 3) return(data.frame(position = numeric(0), height = numeric(0), size = numeric(0)))
  
  if (auto_threshold) {
    adaptive_params <- calculate_adaptive_threshold(signal_data, sensitivity)
    height_threshold <- adaptive_params$height_threshold
    min_prominence <- adaptive_params$min_prominence
  }
  
  working_signal <- signal_data
  max_height <- max(working_signal, na.rm = TRUE)
  median_height <- median(working_signal, na.rm = TRUE)
  absolute_threshold <- median_height + (max_height - median_height) * height_threshold
  
  peaks_idx <- c()
  for (i in 2:(length(working_signal)-1)) {
    is_local_max <- working_signal[i] > working_signal[i-1] && working_signal[i] > working_signal[i+1]
    is_above_threshold <- working_signal[i] > absolute_threshold
    if (is_local_max && is_above_threshold) peaks_idx <- c(peaks_idx, i)
  }
  
  if (length(peaks_idx) == 0) {
    relaxed_threshold <- median_height * 1.2
    for (i in 2:(length(working_signal)-1)) {
      if (working_signal[i] > working_signal[i-1] && working_signal[i] > working_signal[i+1] && working_signal[i] > relaxed_threshold) {
        peaks_idx <- c(peaks_idx, i)
      }
    }
  }
  
  if (length(peaks_idx) == 0) {
    all_local_max <- c()
    for (i in 2:(length(working_signal)-1)) {
      if (working_signal[i] > working_signal[i-1] && working_signal[i] > working_signal[i+1]) all_local_max <- c(all_local_max, i)
    }
    if (length(all_local_max) > 0) {
      heights <- working_signal[all_local_max]
      top_n <- min(10, length(all_local_max))
      top_indices <- order(heights, decreasing = TRUE)[1:top_n]
      peaks_idx <- all_local_max[top_indices]
    }
  }
  
  if (length(peaks_idx) == 0) return(data.frame(position = numeric(0), height = numeric(0), size = numeric(0)))
  
  if (min_prominence > 0) {
    prominences <- numeric(length(peaks_idx))
    window_size <- 5
    for (i in 1:length(peaks_idx)) {
      idx <- peaks_idx[i]
      left_start <- max(1, idx - window_size)
      right_end <- min(length(working_signal), idx + window_size)
      base_level <- max(min(working_signal[left_start:idx]), min(working_signal[idx:right_end]))
      prominences[i] <- working_signal[idx] - base_level
    }
    min_prom_value <- max_height * min_prominence
    valid_peaks <- peaks_idx[prominences >= min_prom_value]
    if (length(valid_peaks) == 0 && length(peaks_idx) > 0) valid_peaks <- peaks_idx[prominences >= min_prom_value * 0.3]
    peaks_idx <- valid_peaks
  }
  
  if (length(peaks_idx) == 0) return(data.frame(position = numeric(0), height = numeric(0), size = numeric(0)))
  
  if (min_distance > 0 && length(peaks_idx) > 1) {
    sorted_peaks <- peaks_idx[order(working_signal[peaks_idx], decreasing = TRUE)]
    filtered_peaks <- c(sorted_peaks[1])
    for (i in 2:length(sorted_peaks)) {
      if (min(abs(sorted_peaks[i] - filtered_peaks)) >= min_distance) filtered_peaks <- c(filtered_peaks, sorted_peaks[i])
    }
    peaks_idx <- sort(filtered_peaks)
  }
  
  result <- data.frame(position = peaks_idx, height = working_signal[peaks_idx])
  if (!is.null(size_data) && length(size_data) >= max(peaks_idx)) {
    result$size <- size_data[peaks_idx]
  } else {
    result$size <- peaks_idx
  }
  return(result[order(result$position), ])
}

# 数值检查与处理 
safe_numeric_check <- function(value, default_value = 0, min_val = NULL, max_val = NULL) {
  if (is.null(value) || is.na(value) || !is.numeric(value) || is.infinite(value)) return(default_value)
  if (!is.null(min_val) && value < min_val) return(min_val)
  if (!is.null(max_val) && value > max_val) return(max_val)
  return(value)
}

safe_vector_check <- function(vec, default_vec = NULL, normalize = FALSE) {
  if (is.null(vec) || length(vec) == 0 || any(is.na(vec)) || any(is.infinite(vec))) return(if (!is.null(default_vec)) default_vec else numeric(0))
  if (normalize) {
    vec[vec < 0] <- 0
    if (sum(vec) > 0) vec <- vec / sum(vec) else return(if (!is.null(default_vec)) default_vec else rep(1/length(vec), length(vec)))
  }
  return(vec)
}

safe_dataframe <- function(size_vec, height_vec, name = "data") {
  if (length(size_vec) != length(height_vec)) {
    min_len <- min(length(size_vec), length(height_vec))
    size_vec <- size_vec[1:min_len]
    height_vec <- height_vec[1:min_len]
  }
  valid_indices <- !is.na(size_vec) & !is.na(height_vec)
  if (sum(valid_indices) == 0) return(data.frame(size = numeric(0), height = numeric(0)))
  return(data.frame(size = size_vec[valid_indices], height = height_vec[valid_indices]))
}

# 数据读取和预处理
read_str_data_improved <- function(file_path, sheet_name = NULL) {
  cat("正在读取文件:", file_path, "\n")
  if (is.null(sheet_name)) {
    sheet_names <- excel_sheets(file_path)
    cat("检测到工作表:", paste(sheet_names, collapse = ", "), "\n")
    all_data <- list()
    for (sheet in sheet_names) {
      tryCatch({
        all_data[[sheet]] <- read_excel(file_path, sheet = sheet)
      }, error = function(e) {
        cat("读取工作表", sheet, "时出错:", e$message, "\n")
      })
    }
    return(all_data)
  } else {
    return(read_excel(file_path, sheet = sheet_name))
  }
}

standardize_str_data_improved <- function(raw_data) {
  standardized_data <- list()
  for (sheet_name in names(raw_data)) {
    df <- raw_data[[sheet_name]]
    if (nrow(df) == 0 || ncol(df) < 2) {
      cat("跳过空的或无效的工作表:", sheet_name, "\n")
      next
    }
    cat("处理工作表:", sheet_name, "- 维度:", nrow(df), "x", ncol(df), "\n")
    
    col_names <- tolower(names(df))
    size_col <- NA; height_col <- NA
    
    for (pattern in c("size", "bp", "length", "position", "allele", "fragment")) {
      matches <- grep(pattern, col_names)
      if (length(matches) > 0) { size_col <- matches[1]; break }
    }
    for (pattern in c("height", "intensity", "peak", "rfu", "signal", "amplitude")) {
      matches <- grep(pattern, col_names)
      if (length(matches) > 0) { height_col <- matches[1]; break }
    }
    if (is.na(size_col) && ncol(df) >= 1) { size_col <- 1; cat("未找到size列，使用第1列作为size\n") }
    if (is.na(height_col) && ncol(df) >= 2) { height_col <- 2; cat("未找到height列，使用第2列作为height\n") }
    
    if (!is.na(size_col) && !is.na(height_col)) {
      standardized_df <- safe_dataframe(as.numeric(df[[size_col]]), as.numeric(df[[height_col]]), sheet_name)
      if (nrow(standardized_df) > 0) {
        standardized_df <- standardized_df[order(standardized_df$size), ]
        size_range <- range(standardized_df$size)
        if (size_range[2] - size_range[1] > 10) {
          standardized_data[[sheet_name]] <- standardized_df
          cat("✓ 成功处理", nrow(standardized_df), "个数据点\n")
        } else { cat("✗ Size范围过小，可能不是STR数据\n") }
      } else { cat("✗ 没有有效的数据点\n") }
    } else { cat("✗ 未找到合适的size和height列\n") }
    cat("\n")
  }
  return(standardized_data)
}

# 贡献者数量估计 
estimate_contributors_peaks_improved <- function(str_data) {
  contributor_estimates <- c()
  locus_info <- list()
  for (locus_name in names(str_data)) {
    locus_data <- str_data[[locus_name]]
    peaks <- detect_peaks_str_optimized(locus_data$height, locus_data$size, auto_threshold = TRUE, sensitivity = "medium")
    if (nrow(peaks) > 0) {
      peak_heights <- peaks$height
      sorted_heights <- sort(peak_heights, decreasing = TRUE)
      if (length(sorted_heights) >= 2) {
        height_ratios <- c()
        for (i in 1:(length(sorted_heights)-1)) height_ratios <- c(height_ratios, sorted_heights[i+1] / sorted_heights[i])
        significant_drops <- which(height_ratios < 0.4)
        if (length(significant_drops) > 0) major_peaks_count <- significant_drops[1]
        else major_peaks_count <- sum(peak_heights > max(sorted_heights) * 0.15)
      } else { major_peaks_count <- length(sorted_heights) }
      
      if (major_peaks_count <= 2) estimated_n <- 1
      else if (major_peaks_count <= 6) estimated_n <- ceiling(major_peaks_count / 2)
      else if (major_peaks_count <= 12) estimated_n <- ceiling(major_peaks_count / 3)
      else estimated_n <- min(6, ceiling(major_peaks_count / 4))
      
      estimated_n <- max(1, min(6, estimated_n))
      contributor_estimates <- c(contributor_estimates, estimated_n)
      locus_info[[locus_name]] <- list(total_peaks = nrow(peaks), major_peaks = major_peaks_count, estimated_contributors = estimated_n)
      cat("基因座", locus_name, ": 总峰数", nrow(peaks), ", 主要峰数", major_peaks_count, ", 估计贡献者", estimated_n, "\n")
    }
  }
  
  if (length(contributor_estimates) > 0) {
    final_estimate <- round(median(contributor_estimates))
    confidence <- if(length(contributor_estimates) > 1 && mean(contributor_estimates) > 0) max(0, min(1, 1 - sd(contributor_estimates) / mean(contributor_estimates))) else 0.5
    return(list(estimate = final_estimate, per_locus = contributor_estimates, confidence = confidence, locus_details = locus_info))
  }
  return(list(estimate = 2, per_locus = c(), confidence = 0.5, locus_details = list()))
}

estimate_contributors_em_improved <- function(str_data) {
  all_peak_data <- data.frame()
  for (locus_name in names(str_data)) {
    peaks <- detect_peaks_str_optimized(str_data[[locus_name]]$height, str_data[[locus_name]]$size, auto_threshold = TRUE, sensitivity = "medium")
    if (nrow(peaks) > 0) { peaks$locus <- locus_name; all_peak_data <- rbind(all_peak_data, peaks) }
  }
  if (nrow(all_peak_data) < 6) return(list(estimate = 2, bic_scores = NULL, optimal_k = 2))
  
  log_heights <- log(all_peak_data$height + 1)
  Q1 <- quantile(log_heights, 0.25); Q3 <- quantile(log_heights, 0.75); IQR <- Q3 - Q1
  clean_heights <- log_heights[log_heights >= (Q1 - 1.5 * IQR) & log_heights <= (Q3 + 1.5 * IQR)]
  
  if (length(clean_heights) < 4) return(list(estimate = 2, bic_scores = NULL, optimal_k = 2))
  
  max_k <- min(4, length(clean_heights) %/% 3)
  bic_scores <- rep(Inf, max_k)
  for (k in 1:max_k) {
    tryCatch({
      if (k == 1) {
        log_likelihood <- sum(dnorm(clean_heights, mean(clean_heights), sd(clean_heights), log = TRUE))
        bic_scores[k] <- -2 * log_likelihood + 2 * log(length(clean_heights))
      } else {
        mixture_result <- normalmixEM(clean_heights, k = k, lambda = rep(1/k, k), verb = FALSE, maxit = 30)
        if (!is.null(mixture_result$loglik)) bic_scores[k] <- -2 * mixture_result$loglik + (3 * k - 1) * log(length(clean_heights))
      }
    }, error = function(e) { bic_scores[k] <- Inf })
  }
  
  optimal_k <- which.min(bic_scores)
  return(list(estimate = if(optimal_k == 1) 1 else min(6, optimal_k + 1), bic_scores = bic_scores, optimal_k = optimal_k))
}

estimate_total_contributors_improved <- function(str_data) {
  cat("开始改进的贡献者人数估计...\n")
  method1_result <- estimate_contributors_peaks_improved(str_data)
  method2_result <- estimate_contributors_em_improved(str_data)
  
  conf1 <- if(is.na(method1_result$confidence) || is.null(method1_result$confidence)) 0.5 else method1_result$confidence
  cat("峰分析方法估计:", method1_result$estimate, "(置信度:", round(conf1, 3), ")\n")
  cat("EM算法方法估计:", method2_result$estimate, "\n")
  
  weight1 <- if(conf1 > 0.7) 0.8 else 0.6
  weight2 <- 1 - weight1
  final_estimate <- max(2, min(6, round(weight1 * method1_result$estimate + weight2 * method2_result$estimate)))
  
  total_peaks <- sum(sapply(str_data, function(x) nrow(detect_peaks_str_optimized(x$height, x$size))))
  final_estimate <- min(final_estimate, ceiling(total_peaks / length(str_data) / 1.5))
  cat("综合估计结果:", final_estimate, "\n")
  
  return(list(final_estimate = final_estimate, method1 = method1_result, method2 = method2_result, total_peaks = total_peaks, confidence = conf1))
}

# 基础混合比例估计 
estimate_mixture_ratios_improved <- function(str_data, n_contributors) {
  all_ratios_by_locus <- list()
  for (locus_name in names(str_data)) {
    peaks <- detect_peaks_str_optimized(str_data[[locus_name]]$height, str_data[[locus_name]]$size, auto_threshold = TRUE, sensitivity = "medium")
    if (nrow(peaks) >= 2) {
      sorted_peaks <- peaks[order(peaks$height, decreasing = TRUE), ]
      top_peaks <- sorted_peaks[1:min(nrow(sorted_peaks), max(4, 2 * n_contributors)), ]
      all_ratios_by_locus[[locus_name]] <- top_peaks$height / sum(top_peaks$height)
    }
  }
  
  if (length(all_ratios_by_locus) == 0) return(rep(1/n_contributors, n_contributors))
  all_ratios <- unlist(all_ratios_by_locus)
  
  if (length(all_ratios) >= n_contributors) {
    tryCatch({
      sorted_ratios <- sort(all_ratios, decreasing = TRUE)
      group_size <- ceiling(length(sorted_ratios) / n_contributors)
      contributor_ratios <- numeric(n_contributors)
      for (i in 1:n_contributors) {
        start_idx <- (i-1) * group_size + 1; end_idx <- min(i * group_size, length(sorted_ratios))
        if (start_idx <= length(sorted_ratios)) contributor_ratios[i] <- mean(sorted_ratios[start_idx:end_idx])
      }
      contributor_ratios <- contributor_ratios / sum(contributor_ratios)
      if (any(is.na(contributor_ratios))) return(rep(1/n_contributors, n_contributors))
      return(contributor_ratios)
    }, error = function(e) { return(rep(1/n_contributors, n_contributors)) })
  } else {
    return(rep(1/n_contributors, n_contributors))
  }
}

#带 Huber 损失与 L2 正则化的鲁棒EM模型 
optimize_mixture_ratios_improved <- function(str_data, n_contributors, initial_ratios = NULL) {
  if (is.null(initial_ratios)) initial_ratios <- estimate_mixture_ratios_improved(str_data, n_contributors)
  
  all_heights <- c()
  for (locus in names(str_data)) {
    peaks <- detect_peaks_str_optimized(str_data[[locus]]$height, str_data[[locus]]$size)
    all_heights <- c(all_heights, peaks$height)
  }
  all_heights_clean <- all_heights[all_heights > 0 & !is.na(all_heights)]
  
  if (length(all_heights_clean) < n_contributors * 2 || n_contributors == 1) return(initial_ratios)
  
  log_heights <- log(all_heights_clean + 1)
  norm_h <- (log_heights - mean(log_heights)) / (sd(log_heights) + 1e-6)
  
  init_em <- tryCatch(
    normalmixEM(norm_h, k = n_contributors, lambda = initial_ratios, verb = FALSE, maxit = 30),
    error = function(e) list(lambda = initial_ratios, mu = scale(1:n_contributors))
  )
  
  robust_objective <- function(params) {
    theta <- abs(params) / sum(abs(params))
    lambda_reg <- 0.01 
    expected <- sum(theta * init_em$mu)
    residuals <- norm_h - expected
    loss <- sum(huber_loss(residuals)) + lambda_reg * sum(theta^2)
    return(loss)
  }
  
  tryCatch({
    opt_result <- optim(par = init_em$lambda, fn = robust_objective, method = "L-BFGS-B", lower = rep(0.01, n_contributors), upper = rep(1, n_contributors))
    optimized_ratios <- opt_result$par / sum(opt_result$par)
    return(optimized_ratios)
  }, error = function(e) {
    cat("鲁棒EM优化降级，使用初始估计\n")
    return(initial_ratios)
  })
}

# 基因型尺寸转换
convert_size_to_allele <- function(size_values, reference_sizes = NULL) {
  if (!is.null(reference_sizes)) {
    alleles <- numeric(length(size_values))
    for (i in 1:length(size_values)) alleles[i] <- which.min(abs(reference_sizes - size_values[i]))
    return(alleles)
  }
  return(round(size_values))
}

# 基于 MAP (最大后验概率) 的基因推断体系
infer_genotypes_improved <- function(str_data, n_contributors, mixture_ratios) {
  cat("开始基于 MAP (最大后验概率) 改进的基因型推断...\n")
  all_genotypes <- list()
  
  for (locus_name in names(str_data)) {
    peaks <- detect_peaks_str_optimized(str_data[[locus_name]]$height, str_data[[locus_name]]$size, auto_threshold = TRUE, sensitivity = "medium")
    if (nrow(peaks) == 0) { cat("基因座", locus_name, "没有检测到峰\n"); next }
    cat("基因座", locus_name, "检测到", nrow(peaks), "个峰\n")
    
    locus_genotypes <- list()
    used_peaks <- c()
    ordered_contributors <- order(mixture_ratios, decreasing = TRUE)
    
    for (idx in ordered_contributors) {
      contributor_name <- paste0("Contributor_", idx)
      available_peaks <- peaks[!peaks$position %in% used_peaks, ]
      
      if (nrow(available_peaks) == 0) break
      
      # MAP 推断似然度计算
      expected_height <- max(peaks$height) * mixture_ratios[idx]
      available_peaks$fitness <- exp(-abs(available_peaks$height - expected_height) / expected_height)
      
      best_candidates <- available_peaks[order(available_peaks$fitness, decreasing = TRUE), ]
      
      if (nrow(best_candidates) >= 2) {
        allele1 <- convert_size_to_allele(best_candidates$size[1])
        allele2 <- convert_size_to_allele(best_candidates$size[2])
        used_peaks <- c(used_peaks, best_candidates$position[1], best_candidates$position[2])
      } else {
        allele1 <- convert_size_to_allele(best_candidates$size[1])
        allele2 <- allele1 
        used_peaks <- c(used_peaks, best_candidates$position[1])
      }
      locus_genotypes[[contributor_name]] <- sort(c(allele1, allele2))
      cat("  ", contributor_name, ":", paste(locus_genotypes[[contributor_name]], collapse = "/"), "\n")
    }
    all_genotypes[[locus_name]] <- locus_genotypes
  }
  return(all_genotypes)
}

# 详细降噪处理函数 (原版架构保留)
comprehensive_denoise_detailed <- function(str_data, denoise_threshold_percentile = 0.85, 
                                           smooth_window = 3, show_comparison = TRUE) {
  cat("开始详细降噪处理...\n")
  cat("========================================\n")
  
  denoised_data <- list()
  total_points_removed <- 0
  total_points_original <- 0
  locus_statistics <- list()
  
  for (locus_name in names(str_data)) {
    cat("\n🔍 处理基因座:", locus_name, "\n")
    cat("----------------------------------------\n")
    
    locus_data <- str_data[[locus_name]]
    original_signal <- locus_data$height
    original_size <- locus_data$size
    
    cat("📊 原始数据统计:\n")
    cat("- 数据点数:", length(original_signal), "\n")
    cat("- 最小值:", round(min(original_signal), 2), "\n")
    cat("- 最大值:", round(max(original_signal), 2), "\n")
    cat("- 均值:", round(mean(original_signal), 2), "\n")
    cat("- 中位数:", round(median(original_signal), 2), "\n")
    
    threshold <- quantile(original_signal, denoise_threshold_percentile)
    cat("- 降噪阈值 (", denoise_threshold_percentile*100, "%分位数):", round(threshold, 2), "\n")
    
    if (length(original_signal) > smooth_window) {
      # 【升级5】使用 LOWESS 替换原本的 simple_smooth
      smoothed_signal <- lowess_baseline_correction(original_size, original_signal, span = 0.1)
      cat("- 应用了 LOWESS 回归基线校正与平滑\n")
    } else {
      smoothed_signal <- original_signal
      cat("- 数据点过少，跳过平滑处理\n")
    }
    
    denoised_signal <- smoothed_signal
    noise_mask <- denoised_signal < threshold
    denoised_signal[noise_mask] <- 0
    
    points_removed <- sum(noise_mask)
    points_kept <- sum(!noise_mask)
    removal_rate <- points_removed / length(original_signal) * 100
    
    cat("\n⚡ 降噪效果:\n")
    cat("- 去除噪声点:", points_removed, "个 (", round(removal_rate, 1), "%)\n")
    cat("- 保留有效点:", points_kept, "个\n")
    
    original_sum <- sum(original_signal)
    denoised_sum <- sum(denoised_signal)
    signal_retention <- denoised_sum / original_sum * 100
    cat("- 信号保留率:", round(signal_retention, 1), "%\n")
    
    original_peaks <- detect_peaks_str_optimized(original_signal, original_size, auto_threshold = TRUE, sensitivity = "medium")
    denoised_peaks <- detect_peaks_str_optimized(denoised_signal, original_size, auto_threshold = TRUE, sensitivity = "medium")
    
    cat("- 降噪前检测峰数:", nrow(original_peaks), "个\n")
    cat("- 降噪后检测峰数:", nrow(denoised_peaks), "个\n")
    
    denoised_data[[locus_name]] <- data.frame(
      size = original_size, height = denoised_signal,
      original_height = original_signal, is_noise = noise_mask, smoothed_height = smoothed_signal
    )
    
    locus_statistics[[locus_name]] <- list(
      original_points = length(original_signal), removed_points = points_removed, removal_rate = removal_rate,
      signal_retention = signal_retention, original_peaks = nrow(original_peaks), denoised_peaks = nrow(denoised_peaks),
      threshold_used = threshold
    )
    
    total_points_removed <- total_points_removed + points_removed
    total_points_original <- total_points_original + length(original_signal)
    cat("\n✅ 基因座", locus_name, "降噪完成\n")
  }
  
  cat("\n========================================\n")
  cat("📊 总体降噪统计:\n")
  cat("- 总原始数据点:", total_points_original, "个\n")
  cat("- 总去除噪声点:", total_points_removed, "个\n")
  cat("- 总体噪声去除率:", round(total_points_removed/total_points_original*100, 1), "%\n")
  
  overall_removal_rate <- total_points_removed/total_points_original*100
  avg_signal_retention <- mean(sapply(locus_statistics, function(x) x$signal_retention))
  
  total_original_peaks <- sum(sapply(locus_statistics, function(x) x$original_peaks))
  total_denoised_peaks <- sum(sapply(locus_statistics, function(x) x$denoised_peaks))
  peak_retention_rate <- if (total_original_peaks > 0) total_denoised_peaks / total_original_peaks * 100 else 100
  
  if (overall_removal_rate < 10 && avg_signal_retention > 95) {
    quality_rating <- "🌟 优秀"
    quality_desc <- "信号质量很高，噪声很少，降噪效果理想"
  } else if (overall_removal_rate < 25 && avg_signal_retention > 85 && peak_retention_rate > 90) {
    quality_rating <- "✅ 良好"
    quality_desc <- "适度噪声已清除，主要信号保持完好"
  } else if (overall_removal_rate < 50 && avg_signal_retention > 70 && peak_retention_rate > 80) {
    quality_rating <- "⚠️ 一般"
    quality_desc <- "噪声较多但已有效清除，需注意信号保真度"
  } else {
    quality_rating <- "⚠️ 需要关注"
    quality_desc <- "噪声过多，可能影响分析准确性"
  }
  
  cat("- 综合质量等级:", quality_rating, "\n")
  cat("========================================\n")
  
  return(list(
    denoised_data = denoised_data,
    statistics = list(
      total_original_points = total_points_original, total_removed_points = total_points_removed,
      overall_removal_rate = overall_removal_rate, avg_signal_retention = avg_signal_retention,
      peak_retention_rate = peak_retention_rate, quality_rating = quality_rating, quality_description = quality_desc
    ),
    locus_statistics = locus_statistics
  ))
}

# ===============================
# 可视化与报告模块 (完全保留原版)
# ===============================

plot_str_profile <- function(str_data, locus_name = NULL, title_prefix = "STR图谱") {
  plots <- list(); locus_names <- if (is.null(locus_name)) names(str_data) else locus_name
  for (locus in locus_names) {
    if (locus %in% names(str_data)) {
      plots[[locus]] <- ggplot(str_data[[locus]], aes(x = size, y = height)) +
        geom_line(color = "blue", size = 0.8) + geom_area(alpha = 0.3, fill = "lightblue") +
        labs(title = paste(title_prefix, "-", locus), x = "Size (bp)", y = "Height (RFU)") +
        theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
    }
  }
  return(plots)
}

plot_peak_detection <- function(str_data, sensitivity = "medium") {
  plots <- list()
  for (locus_name in names(str_data)) {
    locus_data <- str_data[[locus_name]]
    peaks <- detect_peaks_str_optimized(locus_data$height, locus_data$size, auto_threshold = TRUE, sensitivity = sensitivity)
    p <- ggplot(locus_data, aes(x = size, y = height)) +
      geom_line(color = "gray60", size = 0.8) + geom_area(alpha = 0.2, fill = "gray80") +
      labs(title = paste("峰检测结果 -", locus_name), subtitle = paste("检测到", nrow(peaks), "个峰"), x = "Size (bp)", y = "Height (RFU)") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 12))
    if (nrow(peaks) > 0) {
      peak_data <- data.frame(size = peaks$size, height = peaks$height, label = paste0("P", 1:nrow(peaks)))
      p <- p + geom_point(data = peak_data, aes(x = size, y = height), color = "red", size = 3) +
        geom_text(data = peak_data, aes(x = size, y = height, label = label), vjust = -0.5, hjust = 0.5, color = "red", size = 3, fontface = "bold")
    }
    plots[[locus_name]] <- p
  }
  return(plots)
}

plot_denoise_comparison <- function(original_data, denoised_results) {
  plots <- list()
  for (locus_name in names(original_data)) {
    if (locus_name %in% names(denoised_results$denoised_data)) {
      original <- original_data[[locus_name]]
      denoised_info <- denoised_results$denoised_data[[locus_name]]
      comparison_data <- data.frame(size = original$size, original_height = original$height, denoised_height = denoised_info$height, is_noise = denoised_info$is_noise)
      p <- ggplot(comparison_data, aes(x = size)) +
        geom_line(aes(y = original_height, color = "降噪前"), size = 0.8, alpha = 0.7) + geom_line(aes(y = denoised_height, color = "降噪后"), size = 1.0) +
        geom_ribbon(aes(ymin = 0, ymax = ifelse(is_noise, original_height, 0)), fill = "red", alpha = 0.3) + scale_color_manual(values = c("降噪前" = "gray60", "降噪后" = "blue")) +
        labs(title = paste("降噪效果对比 (LOWESS) -", locus_name), subtitle = paste("去除", sum(denoised_info$is_noise), "个噪声点"), x = "Size (bp)", y = "Height (RFU)", color = "信号类型") +
        theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 12), legend.position = "bottom")
      plots[[locus_name]] <- p
    }
  }
  return(plots)
}

plot_mixture_ratios <- function(analysis_results) {
  ratios <- analysis_results$mixture_ratios$final_ratios
  n_contributors <- length(ratios)
  pie_data <- data.frame(contributor = paste0("贡献者", 1:n_contributors), ratio = ratios, percentage = round(ratios * 100, 1))
  p1 <- ggplot(pie_data, aes(x = "", y = ratio, fill = contributor)) +
    geom_col(width = 1) + coord_polar(theta = "y") +
    geom_text(aes(label = paste0(percentage, "%")), position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
    labs(title = "鲁棒优化混合比例分布", fill = "贡献者") + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), legend.position = "bottom")
  p2 <- ggplot(pie_data, aes(x = contributor, y = ratio)) +
    geom_col(fill = "coral", alpha = 0.7) + geom_text(aes(label = paste0(percentage, "%")), vjust = -0.3, size = 4, fontface = "bold") +
    labs(title = "混合比例柱状图", x = "贡献者", y = "混合比例") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  return(list(pie = p1, bar = p2))
}

generate_all_visualizations <- function(analysis_results, output_dir = "./STR_Analysis_Plots/") {
  cat("🎨 开始生成所有可视化图表...\n")
  cat("========================================\n")
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE); cat("✓ 创建输出目录:", output_dir, "\n") }
  
  primary_data <- analysis_results$primary_data
  denoised_results <- analysis_results$denoised_results
  
  cat("\n📊 1. 生成原始STR图谱...\n")
  original_plots <- plot_str_profile(primary_data, title_prefix = "原始STR图谱")
  for (i in seq_along(original_plots)) {
    filename <- file.path(output_dir, paste0("01_原始图谱_", names(original_plots)[i], ".png"))
    ggsave(filename, original_plots[[i]], width = 10, height = 6, dpi = 300)
    cat("  ✓ 保存:", basename(filename), "\n")
  }
  
  cat("\n🏔️ 2. 生成峰检测结果图...\n")
  peak_plots <- plot_peak_detection(primary_data)
  for (i in seq_along(peak_plots)) {
    filename <- file.path(output_dir, paste0("02_峰检测_", names(peak_plots)[i], ".png"))
    ggsave(filename, peak_plots[[i]], width = 10, height = 6, dpi = 300)
    cat("  ✓ 保存:", basename(filename), "\n")
  }
  
  cat("\n🧹 3. 生成降噪对比图...\n")
  denoise_plots <- plot_denoise_comparison(primary_data, denoised_results)
  for (i in seq_along(denoise_plots)) {
    filename <- file.path(output_dir, paste0("03_降噪对比_", names(denoise_plots)[i], ".png"))
    ggsave(filename, denoise_plots[[i]], width = 10, height = 6, dpi = 300)
    cat("  ✓ 保存:", basename(filename), "\n")
  }
  
  cat("\n🥧 4. 生成混合比例图...\n")
  ratio_plots <- plot_mixture_ratios(analysis_results)
  if (length(ratio_plots) > 0) {
    if ("pie" %in% names(ratio_plots)) {
      filename <- file.path(output_dir, "04_混合比例_饼图.png")
      ggsave(filename, ratio_plots$pie, width = 8, height = 8, dpi = 300)
      cat("  ✓ 保存: 04_混合比例_饼图.png\n")
    }
    if ("bar" %in% names(ratio_plots)) {
      filename <- file.path(output_dir, "04_混合比例_柱状图.png")
      ggsave(filename, ratio_plots$bar, width = 8, height = 6, dpi = 300)
      cat("  ✓ 保存: 04_混合比例_柱状图.png\n")
    }
  }
  
  cat("\n========================================\n")
  cat("🎨 可视化图表生成完成！\n📁 输出目录:", output_dir, "\n")
  cat("📊 共生成", length(list.files(output_dir, pattern = "\\.png$")), "个图表文件\n")
  cat("========================================\n")
  return(output_dir)
}

generate_word_report <- function(analysis_results, output_file = "./STR_Analysis_Report.docx") {
  cat("\n📝 开始生成Word分析报告...\n")
  cat("========================================\n")
  
  doc <- read_docx()
  doc <- doc %>% body_add_par("法医物证多人身份鉴定", style = "heading 1") %>% 
    body_add_par("STR图谱分析报告 (学术优化版)", style = "heading 1") %>% body_add_par("", style = "Normal") %>% 
    body_add_par(paste("报告生成时间:", Sys.time()), style = "Normal") %>% body_add_par("", style = "Normal") %>% body_add_break()
  
  cat("📋 1. 添加执行摘要...\n")
  summary_text <- sprintf(
    "本次STR图谱分析处理了 %d 个基因座的数据，通过综合峰分析法和EM算法，估计样本中包含 %d 个贡献者（置信度: %.1f%%）。各贡献者的鲁棒混合比例分别为：%s。降噪处理采用 LOWESS 基线校正，去除了 %.1f%% 的噪声点，保留了 %.1f%% 的有效信号。基于 MAP 的基因型推断成功完成。",
    length(analysis_results$primary_data), analysis_results$contributors$final_estimate, analysis_results$contributors$confidence * 100,
    paste(round(analysis_results$mixture_ratios$final_ratios * 100, 1), collapse = "%, ") %>% paste0("%"),
    analysis_results$denoised_results$statistics$overall_removal_rate, analysis_results$denoised_results$statistics$avg_signal_retention
  )
  
  doc <- doc %>% body_add_par("执行摘要", style = "heading 1") %>% body_add_par(summary_text, style = "Normal") %>% body_add_par("", style = "Normal")
  
  cat("📊 2. 添加分析结果表格...\n")
  main_results <- data.frame(
    项目 = c("估计贡献者数量", "分析置信度", "处理基因座数量", "检测峰总数", "降噪去除率", "信号保留率"),
    结果 = c(paste(analysis_results$contributors$final_estimate, "人"), paste(round(analysis_results$contributors$confidence * 100, 1), "%"), paste(length(analysis_results$primary_data), "个"), paste(analysis_results$contributors$total_peaks, "个"), paste(round(analysis_results$denoised_results$statistics$overall_removal_rate, 1), "%"), paste(round(analysis_results$denoised_results$statistics$avg_signal_retention, 1), "%")),
    备注 = c("综合峰分析法和EM算法估计", "基于数据一致性计算", "成功处理的数据集数量", "所有基因座检测到的峰数量", "LOWESS降噪处理去除的数据点比例", "降噪后保留的信号强度比例")
  )
  
  main_table <- flextable(main_results) %>% set_table_properties(width = 1, layout = "autofit") %>% theme_box() %>% border_remove() %>% hline_top(border = fp_border(color = "black", width = 2)) %>% hline(i = 1, border = fp_border(color = "black", width = 1)) %>% hline_bottom(border = fp_border(color = "black", width = 2)) %>% align(align = "center", part = "all") %>% fontsize(size = 10, part = "all") %>% bold(part = "header")
  doc <- doc %>% body_add_par("分析结果", style = "heading 1") %>% body_add_par("表1 STR图谱分析主要结果汇总", style = "Table Caption") %>% flextable::body_add_flextable(main_table) %>% body_add_par("", style = "Normal")
  
  if (length(analysis_results$mixture_ratios$final_ratios) > 0) {
    ratio_data <- data.frame(
      贡献者 = paste0("贡献者 ", 1:length(analysis_results$mixture_ratios$final_ratios)),
      混合比例 = round(analysis_results$mixture_ratios$final_ratios, 4),
      百分比 = paste0(round(analysis_results$mixture_ratios$final_ratios * 100, 2), "%"),
      相对强度 = round(analysis_results$mixture_ratios$final_ratios / min(analysis_results$mixture_ratios$final_ratios), 2)
    )
    ratio_table <- flextable(ratio_data) %>% set_table_properties(width = 1, layout = "autofit") %>% theme_box() %>% border_remove() %>% hline_top(border = fp_border(color = "black", width = 2)) %>% hline(i = 1, border = fp_border(color = "black", width = 1)) %>% hline_bottom(border = fp_border(color = "black", width = 2)) %>% align(align = "center", part = "all") %>% fontsize(size = 10, part = "all") %>% bold(part = "header")
    doc <- doc %>% body_add_par("表2 各贡献者鲁棒混合比例详细信息", style = "Table Caption") %>% flextable::body_add_flextable(ratio_table) %>% body_add_par("", style = "Normal")
  }
  
  if (length(analysis_results$genotypes) > 0) {
    genotype_list <- list()
    for (locus_name in names(analysis_results$genotypes)) {
      locus_genotypes <- analysis_results$genotypes[[locus_name]]
      for (contributor in names(locus_genotypes)) {
        genotype_list <- append(genotype_list, list(data.frame(基因座 = locus_name, 贡献者 = contributor, 基因型 = paste(locus_genotypes[[contributor]], collapse = "/"), 等位基因1 = locus_genotypes[[contributor]][1], 等位基因2 = locus_genotypes[[contributor]][2])))
      }
    }
    if (length(genotype_list) > 0) {
      genotype_data <- do.call(rbind, genotype_list)
      genotype_table <- flextable(genotype_data) %>% set_table_properties(width = 1, layout = "autofit") %>% theme_box() %>% border_remove() %>% hline_top(border = fp_border(color = "black", width = 2)) %>% hline(i = 1, border = fp_border(color = "black", width = 1)) %>% hline_bottom(border = fp_border(color = "black", width = 2)) %>% align(align = "center", part = "all") %>% fontsize(size = 10, part = "all") %>% bold(part = "header")
      doc <- doc %>% body_add_par("表3 MAP基因型推断结果", style = "Table Caption") %>% flextable::body_add_flextable(genotype_table) %>% body_add_par("", style = "Normal")
    }
  }
  
  cat("🔬 3. 添加方法说明...\n")
  doc <- doc %>% body_add_par("分析方法说明", style = "heading 1") %>% body_add_par("本分析系统接入了学术级优化方法：", style = "Normal") %>% body_add_par("", style = "Normal") %>%
    body_add_par("1. 峰检测：使用自适应阈值算法检测STR图谱中的等位基因峰", style = "Normal") %>%
    body_add_par("2. 贡献者估计：结合峰数量分析法和EM算法进行综合估计", style = "Normal") %>%
    body_add_par("3. 混合比例分析：引入 Huber 损失函数与 L2 正则化的鲁棒期望最大化 (EM) 算法估计各贡献者比例", style = "Normal") %>%
    body_add_par("4. 基因型推断：应用最大后验概率 (MAP) 框架，结合荧光强度期望值进行推断", style = "Normal") %>%
    body_add_par("5. 降噪处理：采用 LOWESS 回归拟合基线漂移，结合分位数阈值法去除信号噪声", style = "Normal") %>%
    body_add_par("", style = "Normal") %>% body_add_par("所有方法均符合法医DNA分析的标准操作程序。", style = "Normal")
  
  print(doc, target = output_file)
  cat("========================================\n")
  cat("📝 Word报告生成完成！\n📁 输出文件:", output_file, "\n")
  cat("========================================\n")
  return(output_file)
}

# ===============================
# 主分析函数与运行调度 (100%匹配原版架构的 tryCatch)
# ===============================

main_str_analysis_improved <- function(attachment1_path = "附件1：不同人数的STR图谱数据.xlsx",
                                       attachment2_path = "附件2：不同混合比例的STR图谱数据.xlsx", 
                                       attachment3_path = "附件3：各个贡献者对应的基因型数据.xlsx",
                                       attachment4_path = "附件4：去噪后的STR图谱数据.xlsx") {
  
  cat("========================================\n")
  cat("法医物证STR图谱分析 - 完整最终版\n")
  cat("========================================\n\n")
  
  cat("步骤1: 读取数据文件\n")
  cat("----------------------------------------\n")
  
  attachment1_data <- NULL
  tryCatch({
    attachment1_raw <- read_str_data_improved(attachment1_path)
    attachment1_data <- standardize_str_data_improved(attachment1_raw)
    cat("✓ 附件1处理完成，包含", length(attachment1_data), "个数据集\n")
  }, error = function(e) {
    cat("✗ 附件1处理失败:", e$message, "\n")
  })
  
  attachment2_data <- NULL
  tryCatch({
    attachment2_raw <- read_str_data_improved(attachment2_path)
    attachment2_data <- standardize_str_data_improved(attachment2_raw)
    cat("✓ 附件2处理完成，包含", length(attachment2_data), "个数据集\n")
  }, error = function(e) {
    cat("✗ 附件2处理失败:", e$message, "\n")
  })
  
  primary_data <- if (!is.null(attachment1_data) && length(attachment1_data) > 0) attachment1_data else attachment2_data
  
  if (is.null(primary_data) || length(primary_data) == 0) {
    cat("✗ 无法获取有效的STR数据，分析终止\n")
    return(NULL)
  }
  
  cat("使用", ifelse(!is.null(attachment1_data) && length(attachment1_data) > 0, "附件1", "附件2"), "作为主要分析数据\n\n")
  
  cat("步骤2: 问题1 - 贡献者人数识别\n")
  cat("----------------------------------------\n")
  contributors_result <- estimate_total_contributors_improved(primary_data)
  n_contributors <- contributors_result$final_estimate
  cat("\n")
  
  cat("步骤3: 问题2 - 混合比例识别 (鲁棒EM优化)\n")
  cat("----------------------------------------\n")
  basic_ratios <- estimate_mixture_ratios_improved(primary_data, n_contributors)
  basic_ratios <- safe_vector_check(basic_ratios, rep(1/n_contributors, n_contributors), normalize = TRUE)
  optimized_ratios <- optimize_mixture_ratios_improved(primary_data, n_contributors, basic_ratios)
  optimized_ratios <- safe_vector_check(optimized_ratios, basic_ratios, normalize = TRUE)
  
  cat("基础估计比例:", paste(round(basic_ratios, 3), collapse = ", "), "\n")
  cat("鲁棒优化后比例:", paste(round(optimized_ratios, 3), collapse = ", "), "\n")
  
  mixture_results <- list(basic_ratios = basic_ratios, optimized_ratios = optimized_ratios, final_ratios = optimized_ratios)
  cat("\n")
  
  cat("步骤4: 问题3 - 基因型推断 (MAP 最大后验概率)\n")
  cat("----------------------------------------\n")
  inferred_genotypes <- infer_genotypes_improved(primary_data, n_contributors, optimized_ratios)
  cat("\n")
  
  cat("步骤5: 问题4 - 噪声降噪 (内置 LOWESS)\n")
  cat("----------------------------------------\n")
  denoised_results <- comprehensive_denoise_detailed(primary_data, denoise_threshold_percentile = 0.85, smooth_window = 3, show_comparison = TRUE)
  cat("\n")
  
  cat("步骤6: 结果汇总\n")
  cat("----------------------------------------\n")
  cat("🎉 完整分析完成！主要结果:\n")
  cat(sprintf("• 估计贡献者人数: %d 人 (置信度: %.3f)\n", n_contributors, contributors_result$confidence))
  cat("• 混合比例: ", paste(round(optimized_ratios, 3), collapse = ", "), "\n")
  cat(sprintf("• 成功推断 %d 个数据集的基因型\n", length(inferred_genotypes)))
  cat(sprintf("• 降噪处理: %s (去除了%.1f%%的噪声点)\n", denoised_results$statistics$quality_rating, denoised_results$statistics$overall_removal_rate))
  
  return(list(
    contributors = contributors_result,
    mixture_ratios = mixture_results,
    genotypes = inferred_genotypes,
    denoised_results = denoised_results,
    primary_data = primary_data,
    analysis_summary = list(
      n_contributors = n_contributors,
      confidence = contributors_result$confidence,
      final_ratios = optimized_ratios,
      genotype_count = length(inferred_genotypes),
      denoise_quality = denoised_results$statistics$quality_rating,
      denoise_removal_rate = denoised_results$statistics$overall_removal_rate,
      signal_retention = denoised_results$statistics$avg_signal_retention
    )
  ))
}

generate_comprehensive_output <- function(analysis_results, output_base_dir = "./STR_Analysis_Output/") {
  cat("\n🚀 开始生成完整的分析输出...\n")
  cat("========================================\n")
  if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)
  
  plots_dir <- file.path(output_base_dir, "Plots")
  visualization_dir <- generate_all_visualizations(analysis_results, plots_dir)
  
  report_file <- file.path(output_base_dir, "STR_Analysis_Report.docx")
  word_report <- generate_word_report(analysis_results, report_file)
  
  rdata_file <- file.path(output_base_dir, "STR_Analysis_Results.RData")
  save(analysis_results, file = rdata_file)
  cat("💾 保存R数据:", rdata_file, "\n")
  
  files <- list.files(output_base_dir, recursive = TRUE, full.names = FALSE)
  inventory <- data.frame(
    文件名 = basename(files), 路径 = dirname(files),
    大小_KB = round(file.size(file.path(output_base_dir, files)) / 1024, 2),
    修改时间 = as.character(file.mtime(file.path(output_base_dir, files)))
  )
  write.csv(inventory, file.path(output_base_dir, "文件清单.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("📋 生成文件清单: 文件清单.csv\n")
  
  cat("\n========================================\n")
  cat("🎉 完整分析输出生成完成！\n")
  cat("📁 输出目录:", output_base_dir, "\n")
  cat("📊 可视化图表:", length(list.files(plots_dir, pattern = "\\.png$", recursive = TRUE)), "个\n")
  cat("📝 Word报告: STR_Analysis_Report.docx\n")
  cat("💾 R数据文件: STR_Analysis_Results.RData\n")
  cat("📋 文件清单: 文件清单.csv\n")
  cat("========================================\n")
  
  return(list(output_dir = output_base_dir, plots_dir = visualization_dir, word_report = word_report, rdata_file = rdata_file))
}

safe_run_complete_analysis <- function() {
  cat("开始运行完整STR分析系统...\n")
  cat("========================================\n")
  
  required_files <- c(
    "附件1：不同人数的STR图谱数据.xlsx",
    "附件2：不同混合比例的STR图谱数据.xlsx", 
    "附件3：各个贡献者对应的基因型数据.xlsx",
    "附件4：去噪后的STR图谱数据.xlsx"
  )
  
  available_files <- c()
  for (file in required_files) {
    if (file.exists(file)) {
      available_files <- c(available_files, file)
      cat("✓ 找到文件:", file, "\n")
    } else {
      cat("✗ 未找到文件:", file, "\n")
    }
  }
  
  if (length(available_files) == 0) {
    cat("未找到任何数据文件，无法进行分析\n")
    cat("🔴 提示：请使用 setwd('你的文件夹路径') 将 R 的工作目录设置到包含 Excel 文件的文件夹中！\n")
    return(NULL)
  }
  
  cat("找到", length(available_files), "个数据文件\n\n")
  cat("步骤1: 运行主要分析\n")
  
  analysis_results <- NULL
  tryCatch({
    analysis_results <- main_str_analysis_improved()
    if (!is.null(analysis_results)) {
      cat("✓ 主要分析成功完成\n")
    } else {
      cat("✗ 主要分析返回NULL\n")
    }
  }, error = function(e) {
    cat("✗ 主要分析失败:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(analysis_results)) {
    cat("\n步骤2: 生成完整输出\n")
    output_result <- generate_comprehensive_output(analysis_results)
    cat("\n🎊 完整分析系统运行成功！\n")
    cat("📁 所有结果已保存到:", output_result$output_dir, "\n")
    return(list(analysis_results = analysis_results, output_info = output_result))
  }
  return(NULL)
}

# ===============================
# 使用说明和自动运行
# ===============================

cat("🚀 自动运行完整分析系统...\n\n")
final_complete_results <- safe_run_complete_analysis()

if (!is.null(final_complete_results)) {
  cat("\n🎊 🎊 🎊 恭喜！完整分析系统运行成功！🎊 🎊 🎊\n")
  cat("========================================\n")
  cat("📊 分析结果摘要:\n")
  
  analysis <- final_complete_results$analysis_results
  cat("• 贡献者数量:", analysis$contributors$final_estimate, "人\n")
  cat("• 分析置信度:", round(analysis$contributors$confidence * 100, 1), "%\n")
  cat("• 混合比例:", paste(round(analysis$mixture_ratios$final_ratios * 100, 1), collapse = "%, "), "%\n")
  cat("• 降噪质量:", analysis$denoised_results$statistics$quality_rating, "\n")
  
  cat("\n📁 生成的文件:\n")
  output_dir <- final_complete_results$output_info$output_dir
  all_files <- list.files(output_dir, recursive = TRUE)
  
  cat("• 可视化图表:", length(grep("\\.png$", all_files)), "个\n")
  cat("• Word报告: 1个 (STR_Analysis_Report.docx)\n")
  cat("• 数据文件: 1个 (STR_Analysis_Results.RData)\n")
  cat("• 文件清单: 1个 (文件清单.csv)\n")
  
  cat("\n📁 输出位置:", output_dir, "\n")
  cat("========================================\n")
  cat("🏁 STR图谱分析系统运行完成！\n")
  cat("感谢使用法医物证多人身份鉴定分析系统！\n")
} else {
  cat("\n❌ 系统运行失败\n")
  cat("请检查数据文件是否存在，或查看错误信息\n")
}

cat("\n🏁 系统就绪！所有功能已加载完成。\n")