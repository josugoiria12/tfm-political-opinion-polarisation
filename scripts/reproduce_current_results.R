
#   1. read the raw POL-AXES CSV;
#   2. estimate the retained ordinal CFA model separately by country;
#   3. save model fit, standardized loadings, factor scores, and validation;
#   4. compute KDE-based DER-style polarisation indices and sensitivity checks;
#   5. compute party-choice validation summaries;
#   6. create compact thesis-ready tables;
#   7. compare key outputs with the previous step-by-step result files.

required_packages <- c("lavaan")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Please install the missing package(s) first: ",
    paste(missing_packages, collapse = ", "),
    "\nExample: install.packages('lavaan')",
    call. = FALSE
  )
}

get_script_dir <- function() {
  file_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- file_args[startsWith(file_args, "--file=")][1]

  if (!is.na(file_arg)) {
    return(dirname(normalizePath(sub("--file=", "", file_arg), mustWork = TRUE)))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) {
      return(frame$ofile)
    }
    NA_character_
  }, character(1))
  frame_file <- frame_files[!is.na(frame_files)][1]
  if (!is.na(frame_file)) {
    return(dirname(normalizePath(frame_file, mustWork = TRUE)))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
data_path <- file.path(script_dir, "POL-AXES data.csv")

if (!file.exists(data_path)) {
  project_root <- normalizePath(getwd(), mustWork = FALSE)
  data_path <- file.path(project_root, "Scripts", "POL-AXES data.csv")
}

if (!file.exists(data_path)) {
  stop(
    "Could not find TFM/Scripts/POL-AXES data.csv. ",
    "Place the POL-AXES CSV in the Scripts folder and run this script from the TFM project folder.",
    call. = FALSE
  )
}

output_dir <- file.path(project_root, "thesis", "results", "reproduce_current_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)

model_items <- c(paste0("opn_1_", 1:6), paste0("opn_2_", 1:5))
required_columns <- c(
  "id", "country", "ideol_2", "weight",
  "ideol_1_spa", "ideol_1_uk", "ideol_1_ger",
  model_items
)

missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0) {
  stop(
    "The data file is missing expected column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

analysis_data <- raw_data[, required_columns]
analysis_data$country_label <- factor(
  analysis_data$country,
  levels = c(0, 1, 2),
  labels = c("Spain", "United Kingdom", "Germany")
)

for (item in model_items) {
  values <- analysis_data[[item]]
  unexpected_values <- sort(unique(values[!is.na(values) & !values %in% 1:5]))

  if (length(unexpected_values) > 0) {
    stop(
      item,
      " has values outside the expected 1-5 ordered response scale: ",
      paste(unexpected_values, collapse = ", "),
      call. = FALSE
    )
  }

  analysis_data[[item]] <- ordered(values, levels = 1:5)
}

parse_numeric_prefix <- function(x) {
  suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*", "\\1", as.character(x))))
}

analysis_data$ideol_2_numeric <- parse_numeric_prefix(analysis_data$ideol_2)

candidate_model <- "
  ideas =~ opn_1_1 + opn_1_2 + opn_1_3 + opn_1_4 + opn_1_5 + opn_1_6
  beliefs =~ opn_1_6 + opn_2_1 + opn_2_2 + opn_2_3 + opn_2_4 + opn_2_5
  ideas ~~ beliefs
"

country_fits <- list()
country_scores <- list()
fit_rows <- list()
loading_rows <- list()

for (country_name in levels(analysis_data$country_label)) {
  country_data <- analysis_data[analysis_data$country_label == country_name, , drop = FALSE]

  fit <- lavaan::cfa(
    model = candidate_model,
    data = country_data,
    ordered = model_items,
    estimator = "WLSMV",
    parameterization = "theta",
    std.lv = TRUE
  )

  country_fits[[country_name]] <- fit

  fit_measures <- lavaan::fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    sample = country_name,
    cfi = as.numeric(fit_measures["cfi"]),
    tli = as.numeric(fit_measures["tli"]),
    rmsea = as.numeric(fit_measures["rmsea"]),
    srmr = as.numeric(fit_measures["srmr"]),
    stringsAsFactors = FALSE
  )

  standardized <- lavaan::standardizedSolution(fit)
  loadings <- standardized[standardized$op == "=~", c("lhs", "rhs", "est.std")]
  loadings$sample <- country_name
  loading_rows[[length(loading_rows) + 1]] <- loadings[, c("sample", "lhs", "rhs", "est.std")]

  scores <- as.data.frame(lavaan::lavPredict(fit, type = "lv"))
  scores$id <- country_data$id
  scores$country <- as.character(country_data$country_label)
  scores$ideol_2 <- country_data$ideol_2_numeric
  scores$weight <- country_data$weight
  country_scores[[country_name]] <- scores[, c("id", "country", "ideol_2", "weight", "ideas", "beliefs")]
}

model_fit <- do.call(rbind, fit_rows)
standardized_loadings <- do.call(rbind, loading_rows)
scores_all <- do.call(rbind, country_scores)

write.csv(model_fit, file.path(output_dir, "candidate_model_fit.csv"), row.names = FALSE)
write.csv(standardized_loadings, file.path(output_dir, "candidate_standardized_loadings.csv"), row.names = FALSE)
write.csv(scores_all, file.path(output_dir, "baseline_factor_scores.csv"), row.names = FALSE)

summarize_factor <- function(data, country_name, factor_name) {
  values <- data[[factor_name]]
  values <- values[!is.na(values)]
  data.frame(
    country = country_name,
    factor = factor_name,
    n = length(values),
    mean = mean(values),
    sd = sd(values),
    variance = var(values),
    iqr = IQR(values),
    mad = mad(values),
    p05 = as.numeric(quantile(values, 0.05)),
    p10 = as.numeric(quantile(values, 0.10)),
    p25 = as.numeric(quantile(values, 0.25)),
    median = median(values),
    p75 = as.numeric(quantile(values, 0.75)),
    p90 = as.numeric(quantile(values, 0.90)),
    p95 = as.numeric(quantile(values, 0.95)),
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
validation_rows <- list()

for (country_name in unique(scores_all$country)) {
  country_data <- scores_all[scores_all$country == country_name, , drop = FALSE]

  for (factor_name in c("ideas", "beliefs")) {
    summary_rows[[length(summary_rows) + 1]] <- summarize_factor(country_data, country_name, factor_name)

    complete_rows <- !is.na(country_data[[factor_name]]) & !is.na(country_data$ideol_2)
    validation_rows[[length(validation_rows) + 1]] <- data.frame(
      country = country_name,
      factor = factor_name,
      validation_variable = "ideol_2",
      n_complete = sum(complete_rows),
      pearson = cor(country_data[[factor_name]][complete_rows], country_data$ideol_2[complete_rows], method = "pearson"),
      spearman = cor(country_data[[factor_name]][complete_rows], country_data$ideol_2[complete_rows], method = "spearman"),
      stringsAsFactors = FALSE
    )
  }
}

score_summaries <- do.call(rbind, summary_rows)
score_validation <- do.call(rbind, validation_rows)

write.csv(score_summaries, file.path(output_dir, "baseline_score_summaries.csv"), row.names = FALSE)
write.csv(score_validation, file.path(output_dir, "baseline_score_validation_ideol_2.csv"), row.names = FALSE)

weighted_quantile <- function(x, w, probs) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  approx(cw, x, xout = probs, ties = "ordered", rule = 2)$y
}

weighted_mean <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_var <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  mu <- sum(x * w) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
}

duclos_bw <- function(x, alpha, w = NULL) {
  keep <- !is.na(x)
  if (!is.null(w)) {
    keep <- keep & !is.na(w) & w > 0
  }
  x <- x[keep]
  n <- length(x)
  if (n < 2 || alpha <= 0) {
    return(NA_real_)
  }
  sigma <- if (is.null(w)) sd(x) else sqrt(weighted_var(x, w[keep]))
  4.7 * n^(-1 / 5) * sigma * alpha^(1 / 5)
}

weighted_kde <- function(x, w = NULL, alpha = 0.50, bw_multiplier = 1, n_grid = 512) {
  keep_x <- !is.na(x)
  x <- x[keep_x]

  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- w[keep_x]
  }

  keep_w <- !is.na(w) & w > 0
  x <- x[keep_w]
  w <- w[keep_w] / sum(w[keep_w])

  bw <- duclos_bw(x, alpha, w) * bw_multiplier
  if (is.na(bw) || bw <= 0) {
    stop("The Duclos bandwidth could not be computed.", call. = FALSE)
  }
  x_range <- range(x)
  padding <- 3 * bw
  grid <- seq(x_range[1] - padding, x_range[2] + padding, length.out = n_grid)
  dx <- grid[2] - grid[1]

  z <- outer(grid, x, "-") / bw
  density_values <- as.vector((dnorm(z) %*% w) / bw)
  density_values <- density_values / sum(density_values * dx)

  list(grid = grid, density = density_values, dx = dx, bw = bw)
}

der_index <- function(kde, alpha) {
  x <- kde$grid
  f <- kde$density
  dx <- kde$dx
  distance_matrix <- abs(outer(x, x, "-"))
  sum((f^(1 + alpha)) * as.vector(distance_matrix %*% f) * dx * dx)
}

mode_count <- function(kde) {
  f <- kde$density
  if (length(f) < 3) {
    return(NA_integer_)
  }
  sum(f[2:(length(f) - 1)] > f[1:(length(f) - 2)] & f[2:(length(f) - 1)] > f[3:length(f)])
}

alphas <- c(0.25, 0.50, 0.75)
factors <- c("ideas", "beliefs")
polarization_rows <- list()
kde_rows <- list()

for (country_name in unique(scores_all$country)) {
  country_data <- scores_all[scores_all$country == country_name, , drop = FALSE]

  for (factor_name in factors) {
    x <- country_data[[factor_name]]
    w <- country_data$weight
    # Alpha 0.50 is the reference density saved for plots and mode counts.
    kde_unweighted <- weighted_kde(x, alpha = 0.50)
    kde_weighted <- weighted_kde(x, w, alpha = 0.50)

    base_row <- data.frame(
      country = country_name,
      factor = factor_name,
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      variance = var(x, na.rm = TRUE),
      iqr = IQR(x, na.rm = TRUE),
      p10 = as.numeric(quantile(x, 0.10, na.rm = TRUE)),
      median = median(x, na.rm = TRUE),
      p90 = as.numeric(quantile(x, 0.90, na.rm = TRUE)),
      weighted_mean = weighted_mean(x, w),
      weighted_sd = sqrt(weighted_var(x, w)),
      weighted_p10 = weighted_quantile(x, w, 0.10),
      weighted_median = weighted_quantile(x, w, 0.50),
      weighted_p90 = weighted_quantile(x, w, 0.90),
      kde_bw = kde_unweighted$bw,
      weighted_kde_bw = kde_weighted$bw,
      kde_modes = mode_count(kde_unweighted),
      weighted_kde_modes = mode_count(kde_weighted),
      stringsAsFactors = FALSE
    )

    for (alpha in alphas) {
      alpha_key <- gsub("\\.", "_", as.character(alpha))
      alpha_kde_unweighted <- weighted_kde(x, alpha = alpha)
      alpha_kde_weighted <- weighted_kde(x, w, alpha = alpha)
      base_row[[paste0("der_alpha_", alpha_key)]] <- der_index(alpha_kde_unweighted, alpha)
      base_row[[paste0("weighted_der_alpha_", alpha_key)]] <- der_index(alpha_kde_weighted, alpha)
    }

    polarization_rows[[length(polarization_rows) + 1]] <- base_row
    kde_rows[[length(kde_rows) + 1]] <- data.frame(
      country = country_name,
      factor = factor_name,
      weighted = FALSE,
      x = kde_unweighted$grid,
      density = kde_unweighted$density
    )
    kde_rows[[length(kde_rows) + 1]] <- data.frame(
      country = country_name,
      factor = factor_name,
      weighted = TRUE,
      x = kde_weighted$grid,
      density = kde_weighted$density
    )
  }
}

polarization_table <- do.call(rbind, polarization_rows)
kde_table <- do.call(rbind, kde_rows)

write.csv(polarization_table, file.path(output_dir, "baseline_polarization_indices.csv"), row.names = FALSE)
write.csv(kde_table, file.path(output_dir, "baseline_kde_values.csv"), row.names = FALSE)

rank_rows <- list()
rank_cols <- c("der_alpha_0_25", "der_alpha_0_5", "der_alpha_0_75")
for (factor_name in factors) {
  factor_table <- polarization_table[polarization_table$factor == factor_name, ]
  for (col in rank_cols) {
    ord <- order(factor_table[[col]], decreasing = TRUE)
    rank_rows[[length(rank_rows) + 1]] <- data.frame(
      factor = factor_name,
      index = col,
      rank = seq_along(ord),
      country = factor_table$country[ord],
      value = factor_table[[col]][ord],
      stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, rank_rows), file.path(output_dir, "baseline_polarization_rankings.csv"), row.names = FALSE)

bandwidth_multipliers <- c(0.75, 1.00, 1.25, 1.50)
sensitivity_rows <- list()
for (country_name in unique(scores_all$country)) {
  country_data <- scores_all[scores_all$country == country_name, , drop = FALSE]

  for (factor_name in factors) {
    x <- country_data[[factor_name]]
    w <- country_data$weight

    for (bw_multiplier in bandwidth_multipliers) {
      for (weighted in c(FALSE, TRUE)) {
        for (alpha in alphas) {
          kde <- if (weighted) {
            weighted_kde(x, w, alpha = alpha, bw_multiplier = bw_multiplier)
          } else {
            weighted_kde(x, alpha = alpha, bw_multiplier = bw_multiplier)
          }

          sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.frame(
            country = country_name,
            factor = factor_name,
            weighted = weighted,
            bw_multiplier = bw_multiplier,
            bandwidth = kde$bw,
            alpha = alpha,
            der_index = der_index(kde, alpha),
            kde_modes = mode_count(kde),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

sensitivity_table <- do.call(rbind, sensitivity_rows)
write.csv(sensitivity_table, file.path(output_dir, "polarization_bandwidth_sensitivity.csv"), row.names = FALSE)

sensitivity_rank_rows <- list()
for (factor_name in factors) {
  for (weighted in c(FALSE, TRUE)) {
    for (bw_multiplier in bandwidth_multipliers) {
      for (alpha in alphas) {
        subset_table <- sensitivity_table[
          sensitivity_table$factor == factor_name &
            sensitivity_table$weighted == weighted &
            sensitivity_table$bw_multiplier == bw_multiplier &
            sensitivity_table$alpha == alpha,
        ]
        ord <- order(subset_table$der_index, decreasing = TRUE)

        sensitivity_rank_rows[[length(sensitivity_rank_rows) + 1]] <- data.frame(
          factor = factor_name,
          weighted = weighted,
          bw_multiplier = bw_multiplier,
          alpha = alpha,
          rank = seq_along(ord),
          country = subset_table$country[ord],
          value = subset_table$der_index[ord],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
write.csv(
  do.call(rbind, sensitivity_rank_rows),
  file.path(output_dir, "polarization_bandwidth_rankings.csv"),
  row.names = FALSE
)

spain_labels <- c(
  "1" = "PP", "2" = "PSOE", "3" = "Vox", "4" = "Sumar", "5" = "ERC",
  "6" = "JxCAT", "7" = "EH Bildu", "8" = "EAJ-PNV", "9" = "Podemos",
  "10" = "BNG", "11" = "Coalicion Canaria", "12" = "UPN",
  "13" = "Se Acabo la Fiesta", "14" = "Other", "15" = "Blank vote",
  "16" = "Null vote", "17" = "Would not vote", "18" = "Don't know",
  "19" = "Prefer not to answer"
)

uk_labels <- c(
  "1" = "Labour", "2" = "Conservative", "3" = "Liberal Democrats",
  "4" = "SNP", "5" = "Sinn Fein", "6" = "DUP", "7" = "Reform UK",
  "8" = "Green", "9" = "Plaid Cymru", "10" = "SDLP", "11" = "Alliance",
  "12" = "TUV", "13" = "UUP", "14" = "Other", "15" = "Blank vote",
  "17" = "Would not vote", "18" = "Don't know", "19" = "Prefer not to answer"
)

germany_labels <- c(
  "1" = "CDU/CSU", "2" = "SPD", "3" = "AfD", "4" = "Greens",
  "5" = "Die Linke", "6" = "FDP", "7" = "BSW", "8" = "SSW",
  "14" = "Other", "17" = "Would not vote", "18" = "Don't know",
  "19" = "Prefer not to answer"
)

substantive_codes <- list(
  "Spain" = as.character(1:13),
  "United Kingdom" = as.character(1:13),
  "Germany" = as.character(c(1:8, 14))
)

get_party <- function(country, spa, uk, ger) {
  if (country == "Spain") return(as.character(spa))
  if (country == "United Kingdom") return(as.character(uk))
  if (country == "Germany") return(as.character(ger))
  NA_character_
}

get_label <- function(country, code) {
  if (is.na(code) || code == "") return(NA_character_)
  if (country == "Spain") return(ifelse(code %in% names(spain_labels), spain_labels[[code]], paste0("Code ", code)))
  if (country == "United Kingdom") return(ifelse(code %in% names(uk_labels), uk_labels[[code]], paste0("Code ", code)))
  if (country == "Germany") return(ifelse(code %in% names(germany_labels), germany_labels[[code]], paste0("Code ", code)))
  NA_character_
}

party_data <- analysis_data[, c("id", "country", "country_label", "ideol_1_spa", "ideol_1_uk", "ideol_1_ger", "ideol_2")]
party_data$party_code <- mapply(
  get_party,
  as.character(party_data$country_label),
  party_data$ideol_1_spa,
  party_data$ideol_1_uk,
  party_data$ideol_1_ger,
  USE.NAMES = FALSE
)
party_data$party_label <- mapply(
  get_label,
  as.character(party_data$country_label),
  party_data$party_code,
  USE.NAMES = FALSE
)
party_data$is_substantive_party <- mapply(
  function(country, code) !is.na(code) && code %in% substantive_codes[[country]],
  as.character(party_data$country_label),
  party_data$party_code,
  USE.NAMES = FALSE
)
party_data$country <- as.character(party_data$country_label)
party_data$ideol_2_numeric <- parse_numeric_prefix(party_data$ideol_2)

merged <- merge(
  scores_all,
  party_data[, c("id", "country", "party_code", "party_label", "is_substantive_party", "ideol_2_numeric")],
  by = c("id", "country"),
  all.x = TRUE
)
write.csv(merged, file.path(output_dir, "baseline_scores_with_party.csv"), row.names = FALSE)

summarize_party <- function(data, country_name, party_code, party_label, substantive) {
  subset_data <- data[data$country == country_name & data$party_code == party_code, , drop = FALSE]
  data.frame(
    country = country_name,
    party_code = party_code,
    party_label = party_label,
    is_substantive_party = substantive,
    n = nrow(subset_data),
    mean_ideas = mean(subset_data$ideas, na.rm = TRUE),
    median_ideas = median(subset_data$ideas, na.rm = TRUE),
    sd_ideas = sd(subset_data$ideas, na.rm = TRUE),
    mean_beliefs = mean(subset_data$beliefs, na.rm = TRUE),
    median_beliefs = median(subset_data$beliefs, na.rm = TRUE),
    sd_beliefs = sd(subset_data$beliefs, na.rm = TRUE),
    mean_ideol_2 = mean(subset_data$ideol_2_numeric, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

party_keys <- unique(merged[, c("country", "party_code", "party_label", "is_substantive_party")])
party_keys <- party_keys[!is.na(party_keys$party_code), ]
party_summary_rows <- list()

for (i in seq_len(nrow(party_keys))) {
  party_summary_rows[[i]] <- summarize_party(
    merged,
    party_keys$country[i],
    party_keys$party_code[i],
    party_keys$party_label[i],
    party_keys$is_substantive_party[i]
  )
}

party_summary <- do.call(rbind, party_summary_rows)
party_summary <- party_summary[order(party_summary$country, party_summary$mean_beliefs), ]
party_summary_main <- party_summary[party_summary$is_substantive_party & party_summary$n >= 30, ]

write.csv(party_summary, file.path(output_dir, "party_factor_score_summary.csv"), row.names = FALSE)
write.csv(party_summary_main, file.path(output_dir, "party_factor_score_summary_substantive_n30.csv"), row.names = FALSE)

round_cols <- function(data, cols, digits = 3) {
  for (col in cols) {
    data[[col]] <- round(data[[col]], digits)
  }
  data
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x
}

write_latex_table <- function(data, file, caption, label, align = NULL) {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", ncol(data) - 1), collapse = ""))
  }

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(latex_escape(names(data)), collapse = " & "),
    "\\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    row_values <- vapply(row, function(value) {
      if (is.numeric(value)) {
        format(value, nsmall = 3, trim = TRUE)
      } else {
        latex_escape(as.character(value))
      }
    }, character(1))
    lines <- c(lines, paste(row_values, collapse = " & "), "\\\\")
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file)
}

fit_table <- model_fit[, c("sample", "cfi", "tli", "rmsea", "srmr")]
names(fit_table) <- c("Country", "CFI", "TLI", "RMSEA", "SRMR")
fit_table <- round_cols(fit_table, c("CFI", "TLI", "RMSEA", "SRMR"))
write.csv(fit_table, file.path(output_dir, "table_model_fit.csv"), row.names = FALSE)
write_latex_table(fit_table, file.path(output_dir, "table_model_fit.tex"), "Fit of the retained baseline measurement model", "tab:baseline-fit")

validation_wide <- reshape(
  score_validation[, c("country", "factor", "pearson")],
  idvar = "country",
  timevar = "factor",
  direction = "wide"
)
names(validation_wide) <- c("Country", "Ideas", "Beliefs")
validation_wide <- round_cols(validation_wide, c("Ideas", "Beliefs"))
write.csv(validation_wide, file.path(output_dir, "table_ideol2_validation.csv"), row.names = FALSE)
write_latex_table(validation_wide, file.path(output_dir, "table_ideol2_validation.tex"), "Correlation between factor scores and left-right self-placement", "tab:ideol2-validation")

pol_table <- polarization_table[, c("country", "factor", "sd", "iqr", "der_alpha_0_5", "weighted_der_alpha_0_5")]
names(pol_table) <- c("Country", "Factor", "SD", "IQR", "DER alpha 0.50", "Weighted DER alpha 0.50")
pol_table <- round_cols(pol_table, c("SD", "IQR", "DER alpha 0.50", "Weighted DER alpha 0.50"))
write.csv(pol_table, file.path(output_dir, "table_polarization.csv"), row.names = FALSE)
write_latex_table(
  pol_table,
  file.path(output_dir, "table_polarization.tex"),
  "Baseline dispersion and polarisation indices",
  "tab:baseline-polarization",
  align = "llrrrr"
)

party_extremes <- list()
for (country_name in unique(party_summary_main$country)) {
  country_party <- party_summary_main[party_summary_main$country == country_name, ]
  low <- country_party[which.min(country_party$mean_beliefs), ]
  high <- country_party[which.max(country_party$mean_beliefs), ]
  party_extremes[[length(party_extremes) + 1]] <- data.frame(
    Country = country_name,
    `Lowest party` = low$party_label,
    `Lowest beliefs` = low$mean_beliefs,
    `Highest party` = high$party_label,
    `Highest beliefs` = high$mean_beliefs,
    check.names = FALSE
  )
}
party_extremes <- do.call(rbind, party_extremes)
party_extremes <- round_cols(party_extremes, c("Lowest beliefs", "Highest beliefs"))
write.csv(party_extremes, file.path(output_dir, "table_party_extremes.csv"), row.names = FALSE)
write_latex_table(
  party_extremes,
  file.path(output_dir, "table_party_extremes.tex"),
  "Party groups at the extremes of the beliefs factor",
  "tab:party-extremes",
  align = "llrlr"
)

write.csv(
  standardized_loadings[, c("sample", "lhs", "rhs", "est.std")],
  file.path(output_dir, "table_standardized_loadings.csv"),
  row.names = FALSE
)

compare_csv <- function(label, new_file, old_file, keys = NULL, tolerance = 1e-9) {
  if (!file.exists(old_file)) {
    return(data.frame(
      output = label,
      status = "old file missing",
      max_abs_diff = NA_real_,
      differing_numeric_cells = NA_integer_,
      old_file = old_file,
      new_file = new_file,
      stringsAsFactors = FALSE
    ))
  }

  new_data <- read.csv(new_file, stringsAsFactors = FALSE, check.names = FALSE)
  old_data <- read.csv(old_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (!is.null(keys)) {
    new_data <- new_data[do.call(order, new_data[keys]), , drop = FALSE]
    old_data <- old_data[do.call(order, old_data[keys]), , drop = FALSE]
  }

  common_cols <- intersect(names(new_data), names(old_data))
  numeric_cols <- common_cols[vapply(new_data[common_cols], is.numeric, logical(1)) &
                                vapply(old_data[common_cols], is.numeric, logical(1))]

  if (nrow(new_data) != nrow(old_data)) {
    status <- "different row count"
  } else if (!setequal(names(new_data), names(old_data))) {
    status <- "different columns"
  } else {
    status <- "same shape"
  }

  max_abs_diff <- 0
  differing_numeric_cells <- 0
  if (length(numeric_cols) > 0 && nrow(new_data) == nrow(old_data)) {
    diffs <- unlist(lapply(numeric_cols, function(col) abs(new_data[[col]] - old_data[[col]])))
    max_abs_diff <- max(diffs, na.rm = TRUE)
    differing_numeric_cells <- sum(diffs > tolerance, na.rm = TRUE)
  }

  if (status == "same shape" && differing_numeric_cells == 0) {
    status <- "match"
  } else if (status == "different columns" && differing_numeric_cells == 0) {
    status <- "numeric match; different columns"
  }

  data.frame(
    output = label,
    status = status,
    max_abs_diff = max_abs_diff,
    differing_numeric_cells = differing_numeric_cells,
    old_file = old_file,
    new_file = new_file,
    stringsAsFactors = FALSE
  )
}

old_results <- file.path(project_root, "thesis", "results")
comparison_rows <- list(
  compare_csv(
    "model fit",
    file.path(output_dir, "candidate_model_fit.csv"),
    file.path(old_results, "fourth_lavaan_candidate_validation", "candidate_model_fit.csv"),
    keys = c("sample")
  ),
  compare_csv(
    "standardized loadings",
    file.path(output_dir, "candidate_standardized_loadings.csv"),
    file.path(old_results, "fourth_lavaan_candidate_validation", "candidate_standardized_loadings.csv"),
    keys = c("sample", "lhs", "rhs")
  ),
  compare_csv(
    "score summaries",
    file.path(output_dir, "baseline_score_summaries.csv"),
    file.path(old_results, "seventh_baseline_score_descriptives", "baseline_score_summaries.csv"),
    keys = c("country", "factor")
  ),
  compare_csv(
    "ideol_2 validation",
    file.path(output_dir, "baseline_score_validation_ideol_2.csv"),
    file.path(old_results, "seventh_baseline_score_descriptives", "baseline_score_validation_ideol_2.csv"),
    keys = c("country", "factor")
  ),
  compare_csv(
    "polarisation indices",
    file.path(output_dir, "baseline_polarization_indices.csv"),
    file.path(old_results, "eighth_baseline_polarization_indices", "baseline_polarization_indices.csv"),
    keys = c("country", "factor")
  ),
  compare_csv(
    "party validation",
    file.path(output_dir, "party_factor_score_summary_substantive_n30.csv"),
    file.path(old_results, "tenth_party_validation", "party_factor_score_summary_substantive_n30.csv"),
    keys = c("country", "party_code")
  ),
  compare_csv(
    "thesis table model fit",
    file.path(output_dir, "table_model_fit.csv"),
    file.path(old_results, "eleventh_results_tables_figures", "table_model_fit.csv"),
    keys = c("Country"),
    tolerance = 1e-12
  ),
  compare_csv(
    "thesis table ideol_2 validation",
    file.path(output_dir, "table_ideol2_validation.csv"),
    file.path(old_results, "eleventh_results_tables_figures", "table_ideol2_validation.csv"),
    keys = c("Country"),
    tolerance = 1e-12
  ),
  compare_csv(
    "thesis table polarisation",
    file.path(output_dir, "table_polarization.csv"),
    file.path(old_results, "eleventh_results_tables_figures", "table_polarization.csv"),
    keys = c("Country", "Factor"),
    tolerance = 1e-12
  )
)

comparison_table <- do.call(rbind, comparison_rows)
write.csv(comparison_table, file.path(output_dir, "comparison_with_previous_outputs.csv"), row.names = FALSE)

cat("\nReproducible current-results pipeline complete.\n")
cat("Input data: ", data_path, "\n", sep = "")
cat("Outputs: ", output_dir, "\n", sep = "")
cat("\nKey files:\n")
cat("- candidate_model_fit.csv\n")
cat("- baseline_score_validation_ideol_2.csv\n")
cat("- baseline_polarization_indices.csv\n")
cat("- polarization_bandwidth_sensitivity.csv\n")
cat("- party_factor_score_summary_substantive_n30.csv\n")
cat("- table_model_fit.csv / table_ideol2_validation.csv / table_polarization.csv\n")
cat("- comparison_with_previous_outputs.csv\n")
