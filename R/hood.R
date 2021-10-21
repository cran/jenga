#' support functions for jenga
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @importFrom narray split
#' @importFrom modeest mlv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd
#' @importFrom scales number
#' @importFrom lubridate is.Date
#' @import purrr


#' @param ts A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param k Positive integer. Number of neighbors to consider when applying kernel average. Min number is 3.
#' @param method Positive integer. Distance method for calculating neighbors. Possibile options are: "euclidean", "manhattan", "chebyshev", "sorensen", "gower", "soergel", "kulczynski_d", "canberra", "lorentzian", "intersection", "non-intersection", "wavehedges", "czekanowski", "motyka", "kulczynski_s", "tanimoto", "ruzicka", "inner_product", "harmonic_mean", "cosine", "hassebrook", "jaccard", "dice",  "fidelity",  "bhattacharyya", "squared_chord", "squared_euclidean", "pearson", "neyman", "squared_chi", "prob_symm", "divergence", "clark", "additive_symm", "taneja", "kumar-johnson", "avg".
#' @param kernel String. DIstribution used to calculate kernel densities. Possible options are: "norm", "cauchy", "logis", "unif", "t", "exp", "lnorm".
#' @param ci Confidence interval.
#' @param deriv Integer vector. Number of differentiation operations to perform for each original time feature. 0 = no change; 1: one diff; 2: two diff.
#' @param n_windows Positive integer. Number of validation tests to measure/sample error.
#' @param mode String. Sequencing method: deterministic ("segmented"), or non-deterministic ("sampled").
#' @param dates Date. Vector with dates for time features.

  hood <- function(ts, seq_len, k, method, kernel, ci, deriv, n_windows, mode, dates)
  {
    ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                         forcat_band = "darkorange", forcat_line = "darkorange", hist_line = "gray43",
                         label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
    {
      all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = c(y_hist, y_forcat))
      forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)

      if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- lower; forcat_data$upper <- upper}

      plot <- ggplot()+geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
      if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes(x = x_forcat, ymin = lower, ymax = upper), alpha = 0.3, fill = forcat_band)}
      plot <- plot + geom_line(data = forcat_data, aes(x = x_forcat, y = y_forcat), color = forcat_line, size = line_size)
      if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
      if(is.null(dbreak)){plot <- plot + xlab(label_x)}
      plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)
      plot <- plot + ylab(label_y)  + theme_bw()
      plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

      return(plot)
    }

    ###
    sequential_kld <- function(m)
    {
      matrix <- t(as.matrix(m))
      n <- nrow(matrix)
      if(n == 1){return(NA)}
      dens <- apply(matrix, 1, function(x) tryCatch(density(x[is.finite(x)], from = min(matrix[is.finite(matrix)]), to = max(matrix[is.finite(matrix)])), error = function(e) NA))
      backward <- dens[-n]
      forward <- dens[-1]

      finite_index <- map2_lgl(forward, backward, ~ is.finite(sum(.x$y * log(.x$y/.y$y))))
      seq_kld <- pmap_dbl(list(forward[finite_index], backward[finite_index]), ~ sum(..1$y * log(..1$y/..2$y)[..3]))
      avg_seq_kld <- tryCatch(round(mean(seq_kld, na.rm = TRUE), 3), error = function(e) NA)

      ratios <- log(dens[[n]]$y/dens[[1]]$y)
      finite_index <- is.finite(ratios)

      end_to_end_kld <- dens[[n]]$y * log(dens[[n]]$y/dens[[1]]$y)
      end_to_end_kld <- tryCatch(round(sum(end_to_end_kld[finite_index]), 3), error = function(e) NA)
      kld_stats <- c(avg_seq_kld, end_to_end_kld)

      return(kld_stats)
    }

    ###
    upside_probability <- function(m)
    {
      matrix <- t(as.matrix(m))
      n <- nrow(matrix)
      if(n == 1){return(NA)}
      growths <- matrix[-1,]/matrix[-n,] - 1
      dens <- apply(growths, 1, function(x) tryCatch(density(x[is.finite(x)], from = min(x[is.finite(x)]), to = max(x[is.finite(x)])), error = function(e) NA))
      not_na <- !is.na(dens)
      avg_upp <- tryCatch(round(mean(map_dbl(dens[not_na], ~ sum(.x$y[.x$x>0])/sum(.x$y))), 3), error = function(e) NA)
      end_growth <- matrix[n,]/matrix[1,] - 1
      end_to_end_dens <- tryCatch(density(end_growth[is.finite(end_growth)], from = min(end_growth[is.finite(end_growth)]), to = max(end_growth[is.finite(end_growth)])), error = function(e) NA)
      if(class(end_to_end_dens) == "density"){last_to_first_upp <- round(sum(end_to_end_dens$y[end_to_end_dens$x>0])/sum(end_to_end_dens$y), 3)} else {last_to_first_upp <- NA}
      upp_stats <- c(avg_upp, last_to_first_upp)
      return(upp_stats)
    }

    ts <- as.data.frame(ts)
    n_ts <- nrow(ts)
    feat_names <- colnames(ts)

    all_positive_check <- map_lgl(ts, ~ ifelse(all(.x >= 0), TRUE, FALSE))
    all_negative_check <- map_lgl(ts, ~ ifelse(all(.x < 0), TRUE, FALSE))

    if(seq_len * k > n_ts){stop("vector length too short for testing with seq_len and k")}
    test_index <- unique(round(seq(seq_len * k, n_ts, length.out = n_windows)))
    if(length(test_index) < n_windows){message("testing on ", length(test_index), " windows")}

    results <- map(test_index, ~ tryCatch(engine(ts[1:.x,, drop = FALSE], seq_len, k, method, kernel, error_measurement = TRUE, deriv, mode), error = function(e) NA))
    not_na <- !is.na(results)
    results <- results[not_na]
    raw_errors <- map(transpose(map(results, ~ narray::split(.x$raw_error, along = 2))), ~ as.data.frame(.x))

    last_prediction <- engine(ts, seq_len, k, method, kernel, error_measurement = FALSE, deriv, mode)$prediction

    window_metrics <- map(transpose(map(results, ~ .x$test_metrics)), ~ Reduce(rbind, .x))
    test_metrics <- map(window_metrics, ~ apply(.x, 2, function(x) mean(x[is.finite(x)], na.rm = TRUE)))
    names(test_metrics) <- feat_names

    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

    sample_pred <- map2(last_prediction, raw_errors, ~ map2_dfc(.x, narray::split(.y, along = 1), ~ .x + sample(unlist(.y), size = 1000, replace = TRUE)))

    sample_pred <- map_if(sample_pred, all_positive_check, ~ {.x[.x < 0] <- 0; return(.x)})
    sample_pred <- map_if(sample_pred, all_negative_check, ~ {.x[.x > 0] <- 0; return(.x)})

    p_stats <- function(x){stats <- round(c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = tryCatch(suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA), skewness = tryCatch(suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA)), 3); return(stats)}

    prediction <- map(sample_pred, ~ as.data.frame(t(apply(.x, 2, p_stats))))
    prediction <- map(prediction, ~ {rownames(.x) <- NULL; return(.x)})
    names(prediction) <- feat_names


    avg_iqr_to_range <- round(map_dbl(prediction, ~ mean((.x[,"75%"] - .x[,"25%"])/(.x[,"max"] - .x[,"min"]))), 3)
    last_to_first_iqr <- round(map_dbl(prediction, ~ (.x[seq_len,"75%"] - .x[seq_len,"25%"])/(.x[1,"75%"] - .x[1,"25%"])), 3)

    pred_stats <- as.data.frame(rbind(avg_iqr_to_range, last_to_first_iqr))
    rownames(pred_stats) <- c("avg_iqr_to_range", "terminal_iqr_ratio")

    kld_stats <- map(sample_pred, ~ tryCatch(sequential_kld(.x), error = function(e) c(NA, NA)))
    kld_stats <- as.data.frame(map(kld_stats, ~.x))
    rownames(kld_stats) <- c("avg_kl_divergence", "terminal_kl_divergence")

    upp_stats <- map(sample_pred, ~ tryCatch(upside_probability(.x), error = function(e) c(NA, NA)))
    upp_stats <- as.data.frame(map(upp_stats, ~.x))
    rownames(upp_stats) <- c("avg_upside_prob", "terminal_upside_prob")

    pred_stats <- rbind(pred_stats, kld_stats, upp_stats)

    x_hist <- 1:n_ts
    x_forcat <- (n_ts + 1):(n_ts + seq_len)

    if(is.Date(dates) & length(dates) == n_ts)
    {
      pred_dates <- tail(dates, seq_len) + seq_len
      prediction <- map(prediction, ~ {rownames(.x) <- pred_dates; return(.x)})
      x_hist <- dates
      x_forcat <- pred_dates
    }

    plot <- pmap(list(prediction, ts, feat_names), ~ ts_graph(x_hist = x_hist, y_hist = ..2, x_forcat = x_forcat, y_forcat = ..1[, 3], lower = ..1[, 1], upper = ..1[, 5], label_y = ..3))

    outcome <- list(prediction = prediction, test_metrics = test_metrics, pred_stats = pred_stats, plot = plot)

    return(outcome)
  }
