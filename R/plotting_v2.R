#' Plot number of valid precursor values per sample
#'
#' @description
#' Creates a bar plot showing, for each sample/column in the selected matrix
#' of a `pSILAC` object (`RIA`, `hol`, `NLI`, `light`, or `heavy`), how many
#' non-missing precursor values are available. This is useful for visually
#' inspecting data completeness and identifying samples with unusually low
#' coverage.
#'
#' @param o
#'   A `pSILAC` object. The function expects one of `o$RIA`, `o$hol`, `o$NLI`,
#'   `o$light`, or `o$heavy` (depending on `method`) to be a matrix or data frame
#'   with precursors in rows and samples in columns.
#' @param c
#'   Character string specifying the fill color of the bars. Defaults to
#'   `"lightblue"`.
#' @param method
#'   Character string specifying which data slot to use. Must be one of
#'   `"RIA"`, `"hol"`, `"NLI"`, `"light"`, or `"heavy"`. Defaults to `"RIA"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts `o$RIA`, `o$hol`, `o$NLI`, `o$light`, or `o$heavy` depending on `method`.
#'   \item Reshapes the selected matrix into long format and removes `NA` values.
#'   \item Computes, for each sample (column), the number of non-missing values (`count`).
#'   \item For `method` in `"NLI"`, `"light"`, or `"heavy"`, log2-transforms values
#'         to compute a per-sample median (kept in the internal summary table but
#'         not visualized).
#' }
#'
#' Only the counts are visualized (one bar per sample). If all values are
#' missing, the function returns an empty plot with labels and issues a warning.
#'
#' The function returns a `ggplot2` object so it can be further customized by
#' the user (e.g., adding themes or modifying labels).
#'
#' @return
#' A `ggplot` object showing the number of valid values per sample (bar chart)
#' for the selected `method`.
#'
#' @examples
#' \dontrun{
#'   plotValidValuesPeptide(ps)                     # default: RIA
#'   plotValidValuesPeptide(ps, method = "hol")
#'   plotValidValuesPeptide(ps, method = "NLI")
#'   plotValidValuesPeptide(ps, method = "light")
#'   plotValidValuesPeptide(ps, method = "heavy")
#'
#'   # Custom bar color:
#'   plotValidValuesPeptide(ps, c = "steelblue")
#' }
#'
#' @export
plotValidValuesPeptide <- function(o, c = "lightblue", method = "RIA") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  
  # 1. Extract the correct data slot based on method
  dat <- if (method == "RIA") {
    o$RIA
  } else if (method == "hol") {
    o$hol
  } else if (method == "NLI") {
    o$NLI
  } else if (method == "light") {
    o$light %>% as.data.frame()
  } else if (method == "heavy") {
    o$heavy %>% as.data.frame()
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', 'NLI', 'light', or 'heavy'.")
  }
  
  # 2. Check if data exists
  if (is.null(dat)) stop(paste0("'o$", method, "' is NULL; no data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, "' has zero rows or columns."))
  
  # 3. Long format + remove NAs
  df_long <- dat %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  if (nrow(df_long) == 0L) {
    warning(paste("All", method, "values are NA; returning an empty plot."))
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("Number of precursors (", method, " values)"),
          x     = "Sample",
          y     = "# of valid values"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # 4. Summarise: count always on non-NA; median optionally on log2 scale
  summary_data <- df_long %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      count = dplyr::n(),
      .groups = "drop"
    )
  
  # Median diagnostic (not plotted): log2 only for NLI/light/heavy
  if (method %in% c("NLI", "light", "heavy")) {
    df_for_median <- df_long
    n_nonpos <- sum(df_for_median$value <= 0, na.rm = TRUE)
    if (n_nonpos > 0L) {
      warning(paste0(
        "Found ", n_nonpos, " non-positive ", method, " values (<= 0). ",
        "These cannot be log2-transformed and will be excluded from the median calculation."
      ))
      df_for_median <- dplyr::filter(df_for_median, value > 0)
    }
    
    median_tbl <- df_for_median %>%
      dplyr::mutate(value = log2(value)) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(median = stats::median(value), .groups = "drop")
    
    summary_data <- dplyr::left_join(summary_data, median_tbl, by = "name")
  } else {
    # For RIA/hol, you can keep median on original scale if you want it available:
    median_tbl <- df_long %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(median = stats::median(value), .groups = "drop")
    summary_data <- dplyr::left_join(summary_data, median_tbl, by = "name")
  }
  
  # 5. Plot
  ymax <- max(summary_data$count) * 1.05
  
  ggplot2::ggplot(summary_data, ggplot2::aes(x = name, y = count, label = count)) +
    ggplot2::geom_col(color = "black", fill = c) +
    ggplot2::geom_text(angle = 90, hjust = -0.1) +
    ggplot2::ylim(c(0, ymax)) +
    ggplot2::labs(
      title = paste0("Number of precursors (", method, " values)"),
      x     = "Sample",
      y     = "# of valid values"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title  = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot number of valid protein values per sample
#'
#' @description
#' Aggregates precursor-level values (RIA, HoL, or NLI) to the protein level by
#' taking the median across peptides per protein and sample. It then creates a
#' bar plot showing, for each sample, how many proteins have at least one
#' quantified (non-missing) value.
#'
#' @param o
#'   A `pSILAC` object. The function expects:
#'   \itemize{
#'     \item A data slot defined by `method` (`o$RIA`, `o$hol`, or `o$NLI`)
#'           as a matrix or data frame with peptides/precursors in rows and
#'           samples in columns.
#'     \item `o$peptides` as a data frame containing at least the columns
#'           `peptide` and `protein` to map rows to protein groups.
#'   }
#' @param c
#'   Character string specifying the fill color of the bars. Defaults to
#'   `"lightyellow"`.
#' @param method
#'   Character string specifying which data slot to use. Must be one of
#'   `"RIA"`, `"hol"`, or `"NLI"`. Defaults to `"RIA"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Joins the selected matrix with `o$peptides` to map peptides to proteins.
#'   \item Aggregates to protein level by computing the median across peptides
#'         for each protein and sample (`na.rm = TRUE`).
#'   \item Reshapes the protein-by-sample matrix into long format and removes `NA` values.
#'   \item Computes, for each sample, the number of proteins with a valid value
#'         (`count`) and the median of log2-transformed values (`median`) for
#'         summary purposes.
#' }
#'
#' Only `count` is visualized (one bar per sample). If all protein-level values
#' are missing, the function returns an empty plot with labels and issues a warning.
#'
#' @return
#' A `ggplot` object showing the number of proteins with valid values per sample
#' (bar chart) for the selected `method`.
#'
#' @export
plotValidValuesProtein <- function(o, c = "lightyellow", method = "RIA") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  
  # 1. Extract the correct data slot
  dat <- if (method == "RIA") {
    o$RIA
  } else if (method == "hol") {
    o$hol
  } else if (method == "NLI") {
    o$NLI
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', or 'NLI'.")
  }
  
  # 2. Safety Checks
  if (is.null(dat)) stop(paste0("'o$", method, "' is NULL; no data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, "' has zero rows/cols."))
  if (is.null(o$peptides)) stop("'o$peptides' is NULL; mapping is required.")
  if (!all(c("peptide", "protein") %in% colnames(o$peptides))) {
    stop("'o$peptides' must contain 'peptide' and 'protein' columns.")
  }
  
  # 3. Aggregate to Protein Level
  summary_data <- dat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("peptide") %>%
    dplyr::left_join(
      y  = o$peptides %>% dplyr::select(peptide, protein),
      by = "peptide"
    ) %>%
    dplyr::select(-peptide) %>%
    dplyr::filter(!is.na(protein)) %>% 
    dplyr::group_by(protein) %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), ~ stats::median(., na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::select(-protein) %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "name",
                        values_to = "value") %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  # 4. Handle empty data
  if (nrow(summary_data) == 0L) {
    warning(paste("All protein-level", method, "values are NA; returning empty plot."))
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste0("Number of proteins (", method, ")"),
                           x = "Sample", y = "# of valid values") +
             ggplot2::theme_classic())
  }
  
  # 5. Final Summarization
  summary_data <- summary_data %>%
    dplyr::mutate(value = log2(value)) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      count  = dplyr::n(),
      median = stats::median(value),
      .groups = "drop"
    )
  
  ymax <- max(summary_data$count) * 1.05
  
  # 6. Plotting
  ggplot2::ggplot(summary_data, ggplot2::aes(x = name, y = count, label = count)) +
    ggplot2::geom_col(color = "black", fill = c) +
    ggplot2::geom_text(angle = 90, hjust = -0.1) +
    ggplot2::ylim(c(0, ymax)) +
    ggplot2::labs(
      title = paste0("Number of proteins (", method, " values)"),
      x     = "Sample",
      y     = "# of valid values"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title    = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot precursor-level value distribution per sample
#'
#' @description
#' Creates boxplots of the precursor-level values for each sample/column in the
#' selected matrix (`RIA`, `hol`, `NLI`, `light`, or `heavy`) of a `pSILAC` object.
#' This helps visually assess the distribution, spread, and potential outliers
#' across samples.
#'
#' @param o
#'   A `pSILAC` object. The function expects one of `o$RIA`, `o$hol`, `o$NLI`,
#'   `o$light`, or `o$heavy` (depending on `method`) to be a matrix or data frame
#'   with precursors (peptides) in rows and samples in columns.
#' @param c
#'   Character string specifying the fill color of the boxplots. Defaults to
#'   `"lightblue"`.
#' @param method
#'   Character string specifying which data slot to plot. Must be one of
#'   `"RIA"`, `"hol"`, `"NLI"`, `"light"`, or `"heavy"`. Defaults to `"RIA"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts `o$RIA`, `o$hol`, `o$NLI`, `o$light`, or `o$heavy` depending on `method`.
#'   \item Reshapes the selected matrix into long format and removes `NA` values.
#'   \item If `method` is `"NLI"`, `"light"`, or `"heavy"`, log2-transforms values prior to plotting.
#'   \item Creates a boxplot for each sample (one box per column); outliers are shown as open circles.
#' }
#'
#' RIA and HoL values are plotted on their original scale. NLI, light, and heavy
#' intensities are plotted on the log2 scale. If all values are missing, the
#' function returns an empty plot with labels and issues a warning.
#'
#' @return
#' A `ggplot` object showing the distribution of precursor-level values per sample
#' (boxplots) for the selected `method`.
#'
#' @examples
#' \dontrun{
#'   plotDistributionPeptide(ps)                     # default: RIA
#'   plotDistributionPeptide(ps, method = "hol")
#'   plotDistributionPeptide(ps, method = "NLI")
#'   plotDistributionPeptide(ps, method = "light")
#'   plotDistributionPeptide(ps, method = "heavy")
#'
#'   # Custom color:
#'   plotDistributionPeptide(ps, c = "skyblue3", method = "RIA")
#' }
#'
#' @export
plotDistributionPeptide <- function(o, c = "lightblue", method = "RIA") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  
  # 1. Extract the correct data slot based on method
  dat <- if (method == "RIA") {
    o$RIA
  } else if (method == "hol") {
    o$hol
  } else if (method == "NLI") {
    o$NLI
  } else if (method == "light") {
    o$light %>% as.data.frame()
  } else if (method == "heavy") {
    o$heavy %>% as.data.frame()
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', 'NLI', 'light', or 'heavy'.")
  }
  
  # 2. Safety checks
  if (is.null(dat)) stop(paste0("'o$", method, "' is NULL; no data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, "' has zero rows or columns."))
  
  # 3. Long format + remove NAs
  df_long <- dat %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  # 4. Handle empty data (all NA)
  if (nrow(df_long) == 0L) {
    warning(paste("All", method, "values are NA; returning an empty plot."))
    
    ylab_empty <- if (method == "RIA") {
      "RIA = L / (L + H)"
    } else if (method == "hol") {
      "HoL = (H / L) + 1"
    } else if (method == "NLI") {
      "log2(NLI)"
    } else if (method == "light") {
      "log2(Light intensity)"
    } else { # heavy
      "log2(Heavy intensity)"
    }
    
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("Data distribution (peptide ", method, ")"),
          x     = "Sample",
          y     = ylab_empty
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # 5. Transform prior to plotting where required
  if (method %in% c("NLI", "light", "heavy")) {
    n_nonpos <- sum(df_long$value <= 0, na.rm = TRUE)
    if (n_nonpos > 0L) {
      warning(paste0(
        "Found ", n_nonpos, " non-positive ", method, " values (<= 0). ",
        "These cannot be log2-transformed and will be removed."
      ))
      df_long <- dplyr::filter(df_long, value > 0)
      if (nrow(df_long) == 0L) {
        warning(paste0("After removing non-positive ", method, " values, no data remain; returning an empty plot."))
        return(
          ggplot2::ggplot() +
            ggplot2::labs(
              title = paste0("Data distribution (peptide ", method, ")"),
              x     = "Sample",
              y     = paste0("log2(", method, ")")
            ) +
            ggplot2::theme_classic()
        )
      }
    }
    df_long <- dplyr::mutate(df_long, value = log2(value))
  }
  
  # 6. Axis labels
  ylab <- if (method == "RIA") {
    "RIA = L / (L + H)"
  } else if (method == "hol") {
    "HoL = (H / L) + 1"
  } else if (method == "NLI") {
    "log2(NLI)"
  } else if (method == "light") {
    "log2(Light intensity)"
  } else { # heavy
    "log2(Heavy intensity)"
  }
  
  # 7. Plot
  ggplot2::ggplot(df_long, ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_boxplot(
      color         = "black",
      fill          = c,
      outlier.shape = 1,
      outlier.alpha = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Data distribution (peptide ", method, ")"),
      x     = "Sample",
      y     = ylab
    ) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot protein-level value distribution per sample
#'
#' @description
#' Aggregates precursor-level values (RIA, HoL, or NLI) to the protein level by
#' taking the median across peptides per protein and sample, and creates boxplots
#' of the resulting protein-level values for each sample. This helps assess the
#' distribution and spread of protein-level values across samples.
#'
#' @param o
#'   A `pSILAC` object. The function expects:
#'   \itemize{
#'     \item A data slot defined by `method` (`o$RIA`, `o$hol`, or `o$NLI`) as a
#'           matrix or data frame with precursors (peptides) in rows and samples
#'           in columns.
#'     \item `o$peptides` as a data frame containing at least the columns
#'           `peptide` and `protein` for mapping precursors to proteins.
#'   }
#' @param c
#'   Character string specifying the fill color of the boxplots. Defaults to
#'   `"lightyellow"`.
#' @param method
#'   Character string specifying which data slot to plot. Must be one of
#'   `"RIA"`, `"hol"`, or `"NLI"`. Defaults to `"RIA"`.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts `o$RIA`, `o$hol`, or `o$NLI` depending on `method`.
#'   \item Joins the selected matrix (rows = peptides) with `o$peptides` to map
#'         each row to a protein.
#'   \item For each protein and sample, computes the median value across all
#'         peptides assigned to that protein (ignoring `NA`s).
#'   \item Reshapes the resulting protein-by-sample matrix into long format,
#'         removes `NA`s, and creates one boxplot per sample/column.
#'   \item If `method = "NLI"`, log2-transforms values prior to plotting.
#' }
#'
#' RIA and HoL values are plotted on their original scale. NLI values are plotted
#' on the log2 scale. Outliers are shown as open circles. If all aggregated
#' protein-level values are missing, the function returns an empty plot with
#' labels and issues a warning.
#'
#' @return
#' A `ggplot` object showing the distribution of protein-level values per sample
#' (boxplots) for the selected `method`.
#'
#' @examples
#' \dontrun{
#'   plotDistributionProtein(ps)                 # default: RIA
#'   plotDistributionProtein(ps, method = "hol")
#'   plotDistributionProtein(ps, method = "NLI")
#'
#'   # Custom color:
#'   plotDistributionProtein(ps, c = "gold", method = "RIA")
#' }
#'
#' @export
plotDistributionProtein <- function(o, c = "lightyellow", method = "RIA") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  
  # 1. Extract the correct data slot based on method
  dat <- if (method == "RIA") {
    o$RIA
  } else if (method == "hol") {
    o$hol
  } else if (method == "NLI") {
    o$NLI
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', or 'NLI'.")
  }
  
  # 2. Safety checks
  if (is.null(dat)) stop(paste0("'o$", method, "' is NULL; no data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, "' has zero rows or columns."))
  if (is.null(o$peptides)) stop("'o$peptides' is NULL; peptide-to-protein mapping is required.")
  if (!all(c("peptide", "protein") %in% colnames(o$peptides))) {
    stop("'o$peptides' must contain 'peptide' and 'protein' columns.")
  }
  
  # 3. Aggregate to protein level
  df_long <- dat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("peptide") %>%
    dplyr::left_join(
      y  = o$peptides %>% dplyr::select(peptide, protein),
      by = "peptide"
    ) %>%
    dplyr::filter(!is.na(protein)) %>%          # avoid protein = NA group
    dplyr::select(-peptide) %>%
    dplyr::group_by(protein) %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), ~ stats::median(., na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::select(-protein) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  # 4. Handle empty data
  if (nrow(df_long) == 0L) {
    warning(paste("All protein-level", method, "values are NA; returning an empty plot."))
    
    ylab_empty <- if (method == "RIA") {
      "RIA = L / (L + H)"
    } else if (method == "hol") {
      "HoL = (H / L) + 1"
    } else {
      "log2(NLI)"
    }
    
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("Data distribution (protein ", method, ")"),
          x     = "Sample",
          y     = ylab_empty
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # 5. Transform NLI prior to plotting
  if (method == "NLI") {
    n_nonpos <- sum(df_long$value <= 0, na.rm = TRUE)
    if (n_nonpos > 0L) {
      warning(paste0(
        "Found ", n_nonpos, " non-positive NLI values (<= 0). ",
        "These cannot be log2-transformed and will be removed."
      ))
      df_long <- dplyr::filter(df_long, value > 0)
      if (nrow(df_long) == 0L) {
        warning("After removing non-positive NLI values, no data remain; returning an empty plot.")
        return(
          ggplot2::ggplot() +
            ggplot2::labs(
              title = "Data distribution (protein NLI)",
              x     = "Sample",
              y     = "log2(NLI)"
            ) +
            ggplot2::theme_classic()
        )
      }
    }
    df_long <- dplyr::mutate(df_long, value = log2(value))
  }
  
  # 6. Axis labels
  ylab <- if (method == "RIA") {
    "RIA = L / (L + H)"
  } else if (method == "hol") {
    "HoL = (H / L) + 1"
  } else {
    "log2(NLI)"
  }
  
  # 7. Plot
  ggplot2::ggplot(df_long, ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_boxplot(
      color         = "black",
      fill          = c,
      outlier.shape = 1,
      outlier.alpha = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Data distribution (protein ", method, ")"),
      x     = "Sample",
      y     = ylab
    ) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot distribution of peptide-level k_loss values per sample
#'
#' @description
#' Creates boxplots of peptide-level k_loss values for each sample/column in the
#' selected peptide-level k_loss table of a `pSILAC` object (`RIA.kloss`,
#' `hol.kloss`, or `NLI.kloss`). This is useful for inspecting the distribution,
#' spread, and potential outliers of k_loss at the peptide level across samples.
#'
#' @param o
#'   A `pSILAC` object. The function expects one of `o$RIA.kloss`, `o$hol.kloss`,
#'   or `o$NLI.kloss` (depending on `method`) to be a data frame containing
#'   peptide-level k_loss estimates in columns ending with `.kloss`.
#' @param c
#'   Character string specifying the fill color of the boxplots. Defaults to
#'   `"lightblue"`.
#' @param method
#'   Character string specifying which peptide-level k_loss table to use. Must be
#'   one of `"RIA"`, `"hol"`, or `"NLI"`. Defaults to `"RIA"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts `o$RIA.kloss`, `o$hol.kloss`, or `o$NLI.kloss` depending on `method`.
#'   \item Selects only columns ending with `.kloss`.
#'   \item Reshapes the selected k_loss columns into long format.
#'   \item Removes `NA` values.
#'   \item Applies a log2 transformation to k_loss values (plots `log2(kloss)`).
#'   \item Produces one boxplot per `.kloss` column, showing the distribution of
#'         log2-transformed k_loss values across peptides.
#' }
#'
#' It is assumed that k_loss values are positive. Non-positive values (<= 0)
#' are incompatible with the log2 transformation and will be removed with a warning.
#'
#' The function returns a `ggplot2` object, which can be further customized.
#'
#' @return
#' A `ggplot` object showing the distribution of peptide-level log2(k_loss)
#' values per sample (boxplots).
#'
#' @examples
#' \dontrun{
#'   # Assuming `ps` is a pSILAC object with peptide-level k_loss tables:
#'   plotDistributionKlossPeptide(ps)                 # default: RIA
#'   plotDistributionKlossPeptide(ps, method = "hol")
#'   plotDistributionKlossPeptide(ps, method = "NLI")
#'
#'   # Custom color:
#'   plotDistributionKlossPeptide(ps, c = "skyblue3")
#' }
#'
#' @export
plotDistributionKlossPeptide <- function(o, c = "lightblue", method = "RIA") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  
  # 1. Extract the correct peptide-level kloss table based on method
  dat <- if (method == "RIA") {
    o$RIA.kloss
  } else if (method == "hol") {
    o$hol.kloss
  } else if (method == "NLI") {
    o$NLI.kloss
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', or 'NLI'.")
  }
  
  # 2. Safety checks
  if (is.null(dat)) stop(paste0("'o$", method, ".kloss' is NULL; no peptide kloss data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, ".kloss' has zero rows or columns."))
  
  # 3. Keep only kloss columns
  kloss_cols <- grep("\\.kloss$", colnames(dat), value = TRUE)
  if (length(kloss_cols) == 0L) {
    stop(paste0("'o$", method, ".kloss' contains no columns ending with '.kloss'."))
  }
  
  df_long <- dat %>%
    dplyr::select(dplyr::all_of(kloss_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  # 4. Handle empty data
  if (nrow(df_long) == 0L) {
    warning(paste0("All peptide ", method, " kloss values are NA; returning an empty plot."))
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("Data distribution (peptide kloss, ", method, ")"),
          x     = "Sample",
          y     = "log2(kloss)"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # 5. Log2 transform (remove non-positive)
  n_nonpos <- sum(df_long$value <= 0, na.rm = TRUE)
  if (n_nonpos > 0L) {
    warning(paste0(
      "Found ", n_nonpos, " non-positive kloss values (<= 0). ",
      "These cannot be log2-transformed and will be removed."
    ))
    df_long <- dplyr::filter(df_long, value > 0)
    if (nrow(df_long) == 0L) {
      warning("After removing non-positive kloss values, no data remain; returning an empty plot.")
      return(
        ggplot2::ggplot() +
          ggplot2::labs(
            title = paste0("Data distribution (peptide kloss, ", method, ")"),
            x     = "Sample",
            y     = "log2(kloss)"
          ) +
          ggplot2::theme_classic()
      )
    }
  }
  
  df_long <- dplyr::mutate(df_long, value = log2(value))
  
  # 6. Plot
  ggplot2::ggplot(df_long, ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_boxplot(
      color         = "black",
      fill          = c,
      outlier.shape = 1,
      outlier.alpha = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Data distribution (peptide kloss, ", method, ")"),
      x     = "Sample",
      y     = "log2(kloss)"
    ) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot peptide-level k_loss QC metrics per sample
#'
#' @description
#' Creates diagnostic QC plots for peptide-level k_loss fitting. For continuous
#' metrics (`kloss.stderr`, `kloss.SSR`, `R2`), the function produces violin plots
#' with a median line per sample. Standard errors and SSR are log2-transformed
#' prior to plotting. R2 is plotted on its original scale (bounded).
#'
#' @param o
#'   A `pSILAC` object. The function expects one of `o$RIA.kloss`, `o$hol.kloss`,
#'   or `o$NLI.kloss` (depending on `method`) to be a data frame containing
#'   peptide-level metrics in columns ending with the relevant suffix.
#' @param c
#'   Character string specifying the fill color of the violins. Defaults to
#'   `"lightblue"`.
#' @param method
#'   Character string specifying which peptide-level k_loss table to use. Must be
#'   one of `"RIA"`, `"hol"`, or `"NLI"`. Defaults to `"RIA"`.
#' @param metric
#'   Character string specifying which QC metric to plot. Must be one of:
#'   `"kloss.stderr"`, `"kloss.SSR"`, or `"R2"`.
#'
#' @details
#' - `kloss.stderr` is extracted from columns ending with `.kloss.stderr` and plotted as `log2(kloss.stderr)`.
#' - `kloss.SSR` is extracted from columns ending with `.kloss.SSR` and plotted as `log2(kloss.SSR)`.
#' - `R2` is extracted from columns ending with `.R_squared` (HoL only) and plotted on the original scale.
#' - Metric suffixes are stripped from column names so the x-axis shows clean sample names.
#'
#' @return
#' A `ggplot` object showing the distribution of the selected QC metric per sample.
#'
#' @export
plotQCkLossPeptide <- function(o,
                               c = "lightblue",
                               method = "RIA",
                               metric = c("kloss.stderr", "kloss.SSR", "R2")) {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  metric <- match.arg(metric)
  
  dat <- if (method == "RIA") {
    o$RIA.kloss
  } else if (method == "hol") {
    o$hol.kloss
  } else if (method == "NLI") {
    o$NLI.kloss
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', or 'NLI'.")
  }
  
  if (is.null(dat)) stop(paste0("'o$", method, ".kloss' is NULL; no peptide QC data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, ".kloss' has zero rows or columns."))
  
  # Map metric -> suffix pattern + transform behavior + labels
  suffix_pattern <- switch(
    metric,
    "kloss.stderr" = "\\.kloss\\.stderr$",
    "kloss.SSR"    = "\\.kloss\\.SSR$",
    "R2"           = "\\.R_squared$"
  )
  
  if (metric == "R2" && method != "hol") {
    stop("metric = 'R2' is only available for method = 'hol'.")
  }
  
  metric_cols <- grep(suffix_pattern, colnames(dat), value = TRUE)
  if (length(metric_cols) == 0L) {
    stop(paste0(
      "No columns found for metric '", metric, "' in 'o$", method,
      ".kloss'. Expected columns matching: ", suffix_pattern
    ))
  }
  
  df_long <- dat %>%
    dplyr::select(dplyr::all_of(metric_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value))
  
  if (nrow(df_long) == 0L) {
    warning(paste0("All peptide ", method, " values for metric '", metric, "' are NA; returning an empty plot."))
    ylab_empty <- if (metric == "R2") "R-squared" else paste0("log2(", metric, ")")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("QC distribution (peptide ", method, " ", metric, ")"),
          x     = "Sample",
          y     = ylab_empty
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # Strip suffix -> clean sample names
  df_long <- df_long %>%
    dplyr::mutate(
      sample = sub(suffix_pattern, "", name),
      sample = factor(sample, levels = unique(sample))
    )
  
  # Transform as requested
  ylab <- metric
  if (metric %in% c("kloss.stderr", "kloss.SSR")) {
    n_nonpos <- sum(df_long$value <= 0, na.rm = TRUE)
    if (n_nonpos > 0L) {
      warning(paste0(
        "Found ", n_nonpos, " non-positive values (<= 0) for metric '", metric, "'. ",
        "These cannot be log2-transformed and will be removed."
      ))
      df_long <- dplyr::filter(df_long, value > 0)
      if (nrow(df_long) == 0L) {
        warning(paste0("After removing non-positive values for metric '", metric, "', no data remain; returning an empty plot."))
        return(
          ggplot2::ggplot() +
            ggplot2::labs(
              title = paste0("QC distribution (peptide ", method, " ", metric, ")"),
              x     = "Sample",
              y     = paste0("log2(", metric, ")")
            ) +
            ggplot2::theme_classic()
        )
      }
    }
    df_long <- dplyr::mutate(df_long, value = log2(value))
    ylab <- paste0("log2(", metric, ")")
  } else if (metric == "R2") {
    ylab <- "R-squared"
  }
  
  # Violin + median line
  ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_violin(fill = c, color = "black", trim = TRUE, scale = "width") +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = 0.35, color = "black") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("QC distribution (peptide ", method, " ", metric, ")"),
      x     = "Sample",
      y     = ylab
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title  = ggplot2::element_text(hjust = 0.5)
    )
}


#' Plot peptide-level nbpoints QC as 100 percent stacked bars per sample
#'
#' @description
#' Creates a 100 percent stacked bar plot showing, for each sample, the proportion of
#' peptides fitted with each discrete nbpoints value (e.g., 3, 4, 5, ...).
#' This scales well to thousands of peptides and enables clear sample-to-sample
#' comparison. The legend is placed at the top.
#'
#' @param o
#'   A `pSILAC` object. The function expects one of `o$RIA.kloss`, `o$hol.kloss`,
#'   or `o$NLI.kloss` (depending on `method`) to be a data frame containing
#'   peptide-level nbpoints in columns ending with `.nbpoints` or named `nbpoints`.
#' @param method
#'   Character string specifying which peptide-level k_loss table to use. Must be
#'   one of `"RIA"`, `"hol"`, or `"NLI"`. Defaults to `"RIA"`.
#'
#' @return
#' A `ggplot` object showing per-sample nbpoints proportions (100 percent stacked bars).
#'
#' @export
plotQCkLossPeptideNbpoints <- function(o, method = "RIA") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  
  dat <- if (method == "RIA") {
    o$RIA.kloss
  } else if (method == "hol") {
    o$hol.kloss
  } else if (method == "NLI") {
    o$NLI.kloss
  } else {
    stop("Invalid method. Choose 'RIA', 'hol', or 'NLI'.")
  }
  
  if (is.null(dat)) stop(paste0("'o$", method, ".kloss' is NULL; no peptide QC data to plot."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", method, ".kloss' has zero rows or columns."))
  
  suffix_pattern <- "(^nbpoints$|\\.nbpoints$)"
  metric_cols <- grep(suffix_pattern, colnames(dat), value = TRUE)
  if (length(metric_cols) == 0L) {
    stop(paste0(
      "'o$", method, ".kloss' contains no nbpoints columns ",
      "(expected suffix '.nbpoints' or name 'nbpoints')."
    ))
  }
  
  df_long <- dat %>%
    dplyr::select(dplyr::all_of(metric_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "nbpoints"
    ) %>%
    dplyr::filter(!is.na(nbpoints))
  
  if (nrow(df_long) == 0L) {
    warning(paste0("All peptide ", method, " nbpoints values are NA; returning an empty plot."))
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = paste0("QC nbpoints distribution (peptide ", method, ")"),
          x     = "Sample",
          y     = "Proportion"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # Clean sample names: strip '.nbpoints' if present, otherwise keep the name.
  df_long <- df_long %>%
    dplyr::mutate(
      sample = sub("\\.nbpoints$", "", name),
      sample = factor(sample, levels = unique(sample)),
      nbpoints = as.integer(nbpoints)
    )
  
  # Summarise to proportions
  df_prop <- df_long %>%
    dplyr::count(sample, nbpoints, name = "n") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(nbpoints = factor(nbpoints, levels = sort(unique(nbpoints))))
  
  ggplot2::ggplot(df_prop, ggplot2::aes(x = sample, y = prop, fill = nbpoints)) +
    ggplot2::geom_col(color = "black", width = 0.9) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("QC nbpoints distribution (peptide ", method, ")"),
      x     = "Sample",
      y     = "Proportion",
      fill  = "nbpoints"
    ) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x     = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title      = ggplot2::element_text(hjust = 0.5)
    )
}


#' Plot binned peptide-level R2 as 100 percent stacked bars per sample (HoL only)
#'
#' @description
#' Creates a 100 percent stacked bar plot showing, for each sample, the proportion of
#' peptides falling into R2 quality bins. This is useful when R2 values are
#' concentrated near 1 (e.g., median ~0.95) and differences are mainly in the tails.
#'
#' @param o
#'   A `pSILAC` object. The function expects `o$hol.kloss` to contain peptide-level
#'   R2 values in columns ending with `.R_squared`.
#'
#' @return
#' A `ggplot` object showing per-sample proportions of binned R2 (100 percent stacked bars).
#'
#' @export
plotQCkLossPeptideR2Binned <- function(o) {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  if (is.null(o$hol.kloss)) stop("'o$hol.kloss' is NULL; no HoL QC data to plot.")
  if (ncol(o$hol.kloss) == 0L || nrow(o$hol.kloss) == 0L) stop("'o$hol.kloss' has zero rows or columns.")
  
  dat <- o$hol.kloss
  
  suffix_pattern <- "\\.R_squared$"
  r2_cols <- grep(suffix_pattern, colnames(dat), value = TRUE)
  if (length(r2_cols) == 0L) stop("'o$hol.kloss' contains no columns ending with '.R_squared'.")
  
  df_long <- dat %>%
    dplyr::select(dplyr::all_of(r2_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "R2"
    ) %>%
    dplyr::filter(!is.na(R2))
  
  if (nrow(df_long) == 0L) {
    warning("All peptide HoL R2 values are NA; returning an empty plot.")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = "QC R2 binned distribution (peptide hol)",
          x     = "Sample",
          y     = "Proportion"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # Clean sample names + bin R2
  # NOTE: Use plain hyphen '-' (avoid Unicode en dash) for Rd/build robustness.
  bin_levels <- c(">=0.90", "0.80-0.90", "0.60-0.80", "<0.60")
  
  df_long <- df_long %>%
    dplyr::mutate(
      sample = sub(suffix_pattern, "", name),
      R2_bin = dplyr::case_when(
        R2 >= 0.90 ~ ">=0.90",
        R2 >= 0.80 ~ "0.80-0.90",
        R2 >= 0.60 ~ "0.60-0.80",
        TRUE       ~ "<0.60"
      ),
      sample = factor(sample, levels = unique(sample)),
      R2_bin = factor(R2_bin, levels = bin_levels)
    )
  
  # Summarise to proportions
  df_prop <- df_long %>%
    dplyr::count(sample, R2_bin, name = "n") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(df_prop, ggplot2::aes(x = sample, y = prop, fill = R2_bin)) +
    ggplot2::geom_col(color = "black", width = 0.9) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "QC R-squared binned distribution (peptide hol)",
      x     = "Sample",
      y     = "Proportion",
      fill  = "R2 bin"
    ) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x     = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title      = ggplot2::element_text(hjust = 0.5)
    )
}


#' Plot distribution of protein-level k_loss values per sample
#'
#' @description
#' Creates boxplots of protein-level k_loss values for each sample/column in
#' the `protein.kloss` matrix of a `pSILAC` object. This is useful for
#' inspecting the distribution, spread, and potential outliers of k_loss at
#' the protein level across samples.
#'
#' @param o
#'   A `pSILAC` object. The function expects:
#'   \itemize{
#'     \item `o$protein.kloss` as a matrix or data frame with proteins in rows
#'           and samples or conditions in columns, containing k_loss estimates.
#'   }
#' @param c
#'   Character string specifying the fill color of the boxplots. Defaults to
#'   `"lightyellow"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reshapes `o$protein.kloss` into long format.
#'   \item Removes `NA` values.
#'   \item Applies a log2 transformation to k_loss values (i.e. plots
#'         `log2(kloss)` on the y-axis).
#'   \item Produces one boxplot per sample/column, showing the distribution
#'         of log2-transformed k_loss values across proteins.
#' }
#'
#' It is assumed that k_loss values are positive; non-positive values would
#' be incompatible with the log2 transformation and should be removed or
#' handled upstream.
#'
#' The function returns a `ggplot2` object, which can be further customized
#' (e.g. additional themes, labels, or facetting).
#'
#' @return
#' A `ggplot` object showing the distribution of protein-level log2(k_loss)
#' values per sample (boxplots).
#'
#' @examples
#' \dontrun{
#'   # Assuming `ps` is a pSILAC object with protein-level k_loss:
#'   p <- plotDistributionKlossProtein(ps)
#'   p
#'
#'   # Custom color:
#'   plotDistributionKlossProtein(ps, c = "gold")
#' }
#'
#' @export
plotDistributionKlossProtein <- function(o, c = "lightyellow") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  if (is.null(o$protein.kloss)) stop("'o$protein.kloss' is NULL; no protein kloss data to plot.")
  if (ncol(o$protein.kloss) == 0L || nrow(o$protein.kloss) == 0L) {
    stop("'o$protein.kloss' has zero rows or columns.")
  }
  
  df_long <- o$protein.kloss %>%
    dplyr::select(-ends_with("stderr")) %>%
    dplyr::select(-contains("source")) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "name",
                        values_to = "value") %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(value = log2(value))
  
  if (nrow(df_long) == 0L) {
    warning("All protein kloss values are NA; returning an empty plot.")
    p_empty <- ggplot2::ggplot() +
      ggplot2::labs(
        title = "Data distribution (protein kloss)",
        x     = "Sample",
        y     = "log2(kloss)"
      ) +
      ggplot2::theme_classic()
    return(p_empty)
  }
  
  p <- df_long %>%
    ggplot2::ggplot(ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_boxplot(
      color         = "black",
      fill          = c,
      outlier.shape = 1,
      outlier.alpha = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "Data distribution (protein kloss)",
      x     = "Sample",
      y     = "log2(kloss)"
    ) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}

#' Plot distribution of protein-level kdeg values per sample
#'
#' @description
#' Creates boxplots of protein-level kdeg values for each sample/column in the
#' `protein.kdeg` matrix of a `pSILAC` object. This is useful for inspecting the
#' distribution, spread, and potential outliers of kdeg at the protein level
#' across samples.
#'
#' @param o
#'   A `pSILAC` object. The function expects:
#'   \itemize{
#'     \item `o$protein.kdeg` as a data frame (or matrix) with proteins in rows
#'           and samples or conditions in columns, containing kdeg estimates.
#'           All columns should be numeric.
#'   }
#' @param c
#'   Character string specifying the fill color of the boxplots. Defaults to
#'   `"lightyellow"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reshapes `o$protein.kdeg` into long format.
#'   \item Removes `NA` values.
#'   \item Removes kdeg values less than or equal to 0 (with a warning).
#'   \item Applies a log2 transformation to remaining kdeg values.
#'   \item Produces one boxplot per sample/column, showing the distribution of
#'         log2-transformed kdeg values across proteins.
#' }
#'
#' @return
#' A `ggplot` object showing the distribution of protein-level log2(kdeg)
#' values per sample (boxplots).
#'
#' @export
plotDistributionKdegProtein <- function(o, c = "lightyellow") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  if (is.null(o$protein.kdeg)) stop("'o$protein.kdeg' is NULL; no protein kdeg data to plot.")
  if (ncol(o$protein.kdeg) == 0L || nrow(o$protein.kdeg) == 0L) {
    stop("'o$protein.kdeg' has zero rows or columns.")
  }
  
  dat <- o$protein.kdeg
  num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  if (length(num_cols) == 0L) stop("'o$protein.kdeg' contains no numeric columns to plot.")
  
  df_long <- dat %>%
    dplyr::select(dplyr::all_of(num_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  if (nrow(df_long) == 0L) {
    warning("All protein kdeg values are NA; returning an empty plot.")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = "Data distribution (protein kdeg)",
          x     = "Sample",
          y     = "log2(kdeg)"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # Remove non-positive kdeg values
  n_removed <- sum(df_long$value <= 0, na.rm = TRUE)
  if (n_removed > 0L) {
    warning(paste0(
      "Removed ", n_removed, " protein kdeg values <= 0 prior to log2 transformation."
    ))
    df_long <- dplyr::filter(df_long, value > 0)
  }
  
  if (nrow(df_long) == 0L) {
    warning("No protein kdeg values remain after filtering; returning an empty plot.")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = "Data distribution (protein kdeg)",
          x     = "Sample",
          y     = "log2(kdeg)"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  df_long <- dplyr::mutate(df_long, value = log2(value))
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_boxplot(
      color         = "black",
      fill          = c,
      outlier.shape = 1,
      outlier.alpha = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "Data distribution (protein kdeg)",
      x     = "Sample",
      y     = "log2(kdeg)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title  = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot distribution of protein-level half-life values per sample
#'
#' @description
#' Creates boxplots of protein-level half-life values for each sample/column in
#' the `protein.halflife` matrix of a `pSILAC` object. This is useful for inspecting
#' the distribution, spread, and potential outliers of protein half-life across
#' samples.
#'
#' @param o
#'   A `pSILAC` object. The function expects:
#'   \itemize{
#'     \item `o$protein.halflife` as a data frame (or matrix) with proteins in rows
#'           and samples or conditions in columns, containing half-life estimates.
#'           All columns should be numeric.
#'   }
#' @param c
#'   Character string specifying the fill color of the boxplots. Defaults to
#'   `"lightyellow"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reshapes `o$protein.halflife` into long format.
#'   \item Removes `NA` values.
#'   \item Removes half-life values less than or equal to 0 (with a warning).
#'   \item Applies a log2 transformation to remaining half-life values.
#'   \item Produces one boxplot per sample/column, showing the distribution of
#'         log2-transformed protein half-life values across proteins.
#' }
#'
#' @return
#' A `ggplot` object showing the distribution of protein-level log2(half-life)
#' values per sample (boxplots).
#'
#' @export
plotDistributionHalflifeProtein <- function(o, c = "lightyellow") {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  if (is.null(o$protein.halflife)) stop("'o$protein.halflife' is NULL; no protein half-life data to plot.")
  if (ncol(o$protein.halflife) == 0L || nrow(o$protein.halflife) == 0L) {
    stop("'o$protein.halflife' has zero rows or columns.")
  }
  
  dat <- o$protein.halflife
  num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  if (length(num_cols) == 0L) stop("'o$protein.halflife' contains no numeric columns to plot.")
  
  df_long <- dat %>%
    dplyr::select(dplyr::all_of(num_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "name",
      values_to = "value"
    ) %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::filter(!is.na(value))
  
  if (nrow(df_long) == 0L) {
    warning("All protein half-life values are NA; returning an empty plot.")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = "Data distribution (protein half-life)",
          x     = "Sample",
          y     = "log2(half-life)"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # Remove non-positive half-life values
  n_removed <- sum(df_long$value <= 0, na.rm = TRUE)
  if (n_removed > 0L) {
    warning(paste0(
      "Removed ", n_removed, " protein half-life values <= 0 prior to log2 transformation."
    ))
    df_long <- dplyr::filter(df_long, value > 0)
  }
  
  if (nrow(df_long) == 0L) {
    warning("No protein half-life values remain after filtering; returning an empty plot.")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = "Data distribution (protein half-life)",
          x     = "Sample",
          y     = "log2(half-life)"
        ) +
        ggplot2::theme_classic()
    )
  }
  
  # Log2 transform
  df_long <- dplyr::mutate(df_long, value = log2(value))
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_boxplot(
      color         = "black",
      fill          = c,
      outlier.shape = 1,
      outlier.alpha = 0.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "Data distribution (protein half-life)",
      x     = "Sample",
      y     = "log2(half-life)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.25),
      plot.title  = ggplot2::element_text(hjust = 0.5)
    )
}

#' Plot correlation matrix of protein-level turnover metrics across samples
#'
#' @description
#' Computes a correlation matrix of protein-level turnover metrics across samples
#' and visualizes it using the \pkg{corrplot} package. Supported metrics are
#' protein-level k_loss, kdeg, and half-life. All metrics are log2-transformed
#' after removal of non-positive values prior to correlation.
#'
#' @param o
#'   A `pSILAC` object. The function expects one of the following slots (depending
#'   on `metric`) to be present as a matrix or data frame with proteins in rows
#'   and samples/conditions in columns:
#'   \itemize{
#'     \item `o$protein.kloss` for `metric = "kloss"`
#'     \item `o$protein.kdeg` for `metric = "kdeg"`
#'     \item `o$protein.halflife` for `metric = "halflife"`
#'   }
#' @param metric
#'   Which protein-level metric to correlate. Must be one of `"kloss"`, `"kdeg"`,
#'   or `"halflife"`. Defaults to `"kloss"`.
#' @param method
#'   Correlation method passed to [stats::cor()]. Defaults to `"spearman"`.
#' @param use
#'   Handling of missing values passed to [stats::cor()]. Defaults to
#'   `"pairwise.complete.obs"`.
#' @param ...
#'   Additional arguments passed to [corrplot::corrplot()], allowing
#'   customization of the correlation plot (e.g. `method`, `tl.cex`, etc.).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Selects the requested protein-level data matrix from the `pSILAC` object.
#'   \item For `metric = "kloss"`, removes auxiliary columns such as standard errors
#'         or `source`.
#'   \item Removes all values less than or equal to zero.
#'   \item Applies a log2 transformation to remaining values.
#'   \item Computes a sample-to-sample correlation matrix and visualizes it
#'         using \pkg{corrplot}.
#' }
#'
#' @return
#' Invisibly returns the correlation matrix (a numeric matrix) used for plotting.
#' The correlation plot is drawn as a side effect.
#'
#' @examples
#' \dontrun{
#'   plotCorProtein(ps, metric = "kloss")
#'   plotCorProtein(ps, metric = "kdeg")
#'   plotCorProtein(ps, metric = "halflife")
#' }
#'
#' @importFrom stats cor
#' @importFrom corrplot corrplot
#' @export
plotCorProtein <- function(o,
                           metric = c("kloss", "kdeg", "halflife"),
                           method = "spearman",
                           use = "pairwise.complete.obs",
                           ...) {
  if (!inherits(o, "pSILAC")) stop("'o' should be a pSILAC object.")
  metric <- match.arg(metric)
  
  # 1. Select data slot
  dat <- switch(
    metric,
    "kloss"    = o$protein.kloss,
    "kdeg"     = o$protein.kdeg,
    "halflife" = o$protein.halflife
  )
  
  slot_name <- switch(
    metric,
    "kloss"    = "protein.kloss",
    "kdeg"     = "protein.kdeg",
    "halflife" = "protein.halflife"
  )
  
  if (is.null(dat)) stop(paste0("'o$", slot_name, "' is NULL; no protein data available."))
  if (ncol(dat) == 0L || nrow(dat) == 0L) stop(paste0("'o$", slot_name, "' has zero rows or columns."))
  if (ncol(dat) == 1L) stop(paste0("'o$", slot_name, "' has only one sample."))
  
  dat <- as.data.frame(dat)
  
  # 2. Column selection
  if (metric == "kloss") {
    dat <- dat %>%
      dplyr::select(-dplyr::ends_with("stderr"), -dplyr::ends_with("source"))
  }
  
  num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  if (length(num_cols) < 2L) {
    stop("Need at least two numeric sample columns to compute a correlation matrix.")
  }
  
  mat <- as.matrix(dat[, num_cols, drop = FALSE])
  
  # 3. Remove non-positive values and log2-transform
  n_removed <- sum(mat <= 0, na.rm = TRUE)
  if (n_removed > 0L) {
    warning(paste0(
      "Removed ", n_removed, " ", metric,
      " values <= 0 prior to log2 transformation."
    ))
    mat[mat <= 0] <- NA_real_
  }
  
  mat <- log2(mat)
  
  # 4. Correlation
  cor_mat <- stats::cor(mat, use = use, method = method)
  
  # 5. Plot
  corrplot::corrplot(
    cor_mat,
    method        = "shade",
    addgrid.col   = "white",
    addCoef.col   = "black",
    tl.col        = "black",
    number.digits = 2,
    ...
  )
  
  invisible(cor_mat)
}



