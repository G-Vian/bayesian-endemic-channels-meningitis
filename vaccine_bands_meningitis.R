# =============================================================================
#  VACCINE-ADJUSTED BANDS FOR MENINGITIS — State of São Paulo, Brazil
# =============================================================================
#
#  WHAT THIS SCRIPT DOES
#  ---------------------------------------------------------------------------
#  It answers a counterfactual question: for a given target year, how many
#  meningitis cases would we have expected WITH the observed vaccination
#  coverage, versus WITHOUT any vaccine-coverage information? For each target and
#  year, two Bayesian models are fitted:
#
#    (1) NO-vaccine baseline  weekly Poisson -> cumulative, with seasonality and
#                             a year effect only. Structurally identical to the
#                             historical-band script.
#    (2) WITH-vaccine model   the SAME model plus vaccine-coverage covariates
#                             (primary series and, where applicable, booster),
#                             entered as standardised z-scores.
#
#  Both produce a full posterior sample matrix (52 weeks x N samples) of the
#  CUMULATIVE trajectory; comparing the two gives the estimated effect of
#  vaccination coverage on expected burden.
#
#  COVERAGE LAG  (config$vaccine$coverage_lag_years)
#  ---------------------------------------------------------------------------
#  Cases in year X are explained by coverage MEASURED in year X - lag. With the
#  default lag = 1, predicting 2022 uses coverage measured in 2021, and each
#  training year is likewise matched to its own previous-year coverage. Setting
#  lag = 0 restores a contemporaneous design (same-year coverage explaining
#  same-year cases) -- a valid attribution design, but NOT a genuine forecast,
#  since it would require knowing the target year's coverage in advance.
#
#  AGGREGATE TARGETS
#  ---------------------------------------------------------------------------
#  For composite targets (Total / Bacterial / Meningococcal) the vaccine effect is
#  not re-estimated from scratch. Component posterior matrices are added and
#  subtracted SAMPLE BY SAMPLE:
#
#      total_WITH = total_baseline - sum(components_baseline) + sum(components_WITH)
#
#  Operating sample-by-sample rather than on quantiles is what keeps uncertainty
#  propagation valid. The 'effect_for_totals' switch controls which components are
#  allowed to contribute to that adjustment.
#
#  DIAGNOSTIC LOG   <<< READ THIS IF A BAND LOOKS IMPLAUSIBLY WIDE >>>
#  ---------------------------------------------------------------------------
#  A vaccine-adjusted band can come out with very large uncertainty for a reason
#  that is invisible in the figure: EXTRAPOLATION in covariate space. If a
#  vaccine's historical coverage was very stable, its training standard deviation
#  is small, so even a moderate real-world drop in the target year becomes an
#  extreme z-score far outside anything the model saw during training. Section 5.5
#  writes a log reporting the training mean/SD, the target-year z-score, an
#  explicit EXTRAPOLATION warning when that z-score is out of range, the vaccine
#  coefficient with its posterior SD, and the WITH/WITHOUT band-width ratio at
#  EW52 -- enough to confirm or rule out that mechanism from real numbers.
#
#  ABOUT THE INPUT DATA   <<< IMPORTANT IF YOU ARE REUSING THIS CODE >>>
#  ---------------------------------------------------------------------------
#  This script reads an ALREADY-CURATED line-list (.RDS): raw SINAN records
#  (Brazil's notifiable disease information system) that have already been
#  cleaned and decoded into readable columns.
#
#  We produced that curated file OURSELVES, starting from RAW SINAN data. Adapting
#  this pipeline to raw SINAN is therefore straightforward: all you need is a
#  preprocessing step that decodes the raw SINAN fields into the column names
#  listed below. Everything downstream consumes those names, not the raw codes.
#
#  EXPECTED COLUMNS in the line-list (one row per notified case):
#    anoepi, sepi     epidemiological year and week; otherwise derived from
#                     data_sintoma (see config$epi_source).
#    data_sintoma     symptom onset date (fallback for year/week).
#    classi_fin       final classification; "Confirmado" marks a confirmed case.
#    grupo_etiologico broad etiological group.
#    etiologia        agent (e.g. "Meningococo", "Pneumococo", "Hemofilo").
#    sorogrupo        meningococcal serogroup.
#    MUNI_RESID       municipality of residence (name).
#    DRS_RESID        health region of residence (name); optional.
#
#  EXPECTED FORMAT of the vaccine-coverage CSVs (one file per vaccine and dose):
#    First column  = municipality, optionally prefixed with its IBGE code.
#    Other columns = one column per YEAR, named with exactly four digits ("2015"),
#                    holding that municipality's coverage for that year.
#    Both ";" and "," separators and decimal commas are handled automatically.
#
#  HOW TO RUN
#  ---------------------------------------------------------------------------
#  1. Replace every "INSERT HERE YOUR DIRECTORY"  placeholder with a real path.
#  2. Set config$vaccine$band_years to your target years, remembering the lag:
#     predicting year X requires coverage measured in year X - lag.
#  3. Source the whole file. Missing packages install automatically.
#
#  REPRODUCIBILITY: all posterior sampling is seeded deterministically from the
#  instance identity (see instance_seed below), so the NO-vaccine baseline yields
#  identical numbers whether it is computed for the figures or re-computed inside
#  the validation module.
#
#  NOTE ON LANGUAGE: some internal variable names are in Portuguese (the working
#  language of the surveillance service); all comments are in English.
# =============================================================================

# =============================================================================
#  1. CONFIGURATION — EDIT HERE
# =============================================================================
config <- list(
  data_path = "INSERT HERE YOUR DIRECTORY" ,
  
  # Output folder exclusively for Vaccine Results
  output_folder = "INSERT HERE YOUR DIRECTORY" ,
  
  # Training years (excluding pandemic 2020/2021). Evaluates targets in 2019/2022.
  train_years = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2022),
  
  epi_source = "auto",
  only_confirmed = TRUE,
  inla_percentiles = c(lower = 0.025, media = 0.50, upper = 0.975),
  n_samples = 1000,
  min_train_cases = 30,
  
  # Epidemiological Targets (Translations maintain strict matching with base data)
  targets = list(
    "Total Meningitis"      = list(tipo = "todas", confirmado = TRUE),
    "Viral Meningitis"      = list(grupo_etiologico = c("Asseptica")),
    "Bacterial Meningitis"  = list(grupo_etiologico = c("Bacteriana", "Nao determinada")),
    "Meningococcal"         = list(etiologia = "Meningococo"),
    "Meningococcal B"       = list(grupo_etiologico = "Bacteriana", etiologia = "Meningococo", sorogrupo = "Sorogrupo B"),
    "Meningococcal C"       = list(grupo_etiologico = "Bacteriana", etiologia = "Meningococo", sorogrupo = "Sorogrupo C"),
    "Pneumococcal"          = list(etiologia = "Pneumococo"),
    "Haemophilus"           = list(etiologia = "Hemofilo")
  ),
  
  inla_models = list(
    State        = list(family = "poisson", seasonal = "rw1", n_harmonics = 2, formula = NULL),
    DRS          = list(family = "poisson", seasonal = "rw1", n_harmonics = 2, formula = NULL),
    CIR          = list(family = "poisson", seasonal = "rw1", n_harmonics = 1, formula = NULL),
    Municipality = list(family = "poisson", seasonal = "rw1", n_harmonics = 1, formula = NULL)
  ),
  
  # Shapefiles for regional mappings
  shapefile_municipalities = "INSERT HERE YOUR DIRECTORY" ,
  col_id_municipalities    = NULL,   # NULL = auto-detects name column
  shapefile_cir            = "INSERT HERE YOUR DIRECTORY" ,
  col_id_cir               = NULL,   # NULL = auto-detects identifier
  shapefile_drs            = "INSERT HERE YOUR DIRECTORY" ,
  col_id_drs               = "X17DRS",
  
  width = 1200, height = 600, dpi = 150,
  
  # ---- VACCINE CONFIGURATION ----------------------------------------------
  vaccine = list(
    folder = "INSERT HERE YOUR DIRECTORY" ,
    files = list(
      penta     = "penta_2012_2022.csv",
      tetra     = "tetra_2002_2016.csv",
      menC      = "menC_2010_2022.csv",
      menC_PR   = "menC_PR_2013_2022.csv",
      pneumo    = "pneumo_2010_2022.csv",
      pneumo_PR = "pneumo_PR_2013_2022.csv"
    ),
    
    cov_by_target = list(
      "Haemophilus"     = list(primary = c("penta", "tetra"), booster = character(0)),
      "Meningococcal C" = list(primary = c("menC"),           booster = c("menC_PR")),
      "Pneumococcal"    = list(primary = c("pneumo"),         booster = c("pneumo_PR"))
    ),
    
    # Target years for which bands are produced. Each must have coverage data
    # available for year (target - coverage_lag_years); see the lag note below.
    band_years = c(2019, 2022),
    aggregation  = "simple_mean",   # "weighted_mean" requires target_pop
    target_pop   = NULL,
    levels       = c("State"),
    standardize  = TRUE,
    skip_municipality = TRUE,       # TRUE = ignores Municipality level to save time
    
    # ---- COVERAGE LAG ------------------------------------------------------
    # coverage_lag_years = N  ->  coverage MEASURED in year Y is used to
    #   explain/predict CASES in year Y + N. Applied consistently in BOTH the
    #   training merge (run_inla_acum_mat_vac) and the target-year covariate
    #   lookup (get_year_covariate) -- see shift_coverage_years() in section 5.
    #   0 = old contemporaneous design (same-year coverage <-> same-year cases).
    #   1 = "coverage of year X-1 predicts cases of year X" (current decision).
    # IMPORTANT: the last usable target year is bounded by your coverage files.
    # With lag = 1, predicting year X requires coverage measured in X-1, so if
    # your CSVs end in year Y then Y+1 is the last feasible target year.
    coverage_lag_years = 1L,
    
    # ---- DIAGNOSTIC LOG -----------------------------------------------------
    # When TRUE, writes a per-instance log (level/region/target/year) reporting:
    # mean and SD of coverage across the training years, the z-score used for the
    # target year, whether that z-score falls OUTSIDE the range seen in training
    # (i.e. EXTRAPOLATION), the estimated vaccine coefficient with its posterior
    # SD, and the band width at EW52 comparing WITH vs WITHOUT vaccine. This lets
    # abnormally wide bands be diagnosed from real numbers rather than guesswork.
    diagnostic_log = TRUE,
    diagnostic_log_file = "vaccine_diagnostic_log.txt",
    
    # ---- BESAG CONFIGURATION ----------------------------------------------
    besag_enable = FALSE,
    besag_levels = c("DRS"),
    
    # ---- AGGREGATE CALCULATION SWITCH -------------------------------------
    # OPTIONS: 
    # "all"       = Uses Hib, MenC, and Pneumo vaccine effects in totals.
    # "menC_only" = Applies ONLY the MenC vaccine effect to the aggregate targets 
    #               (Total Meningitis, Bacterial Meningitis, Meningococcal).
    effect_for_totals = "menC_only"
  )
)

`%||%` <- function(a, b) if (is.null(a)) b else a

# REPRODUCIBILITY: deterministic per-instance seed.
# inla.posterior.sample() draws random samples from the fitted posterior, so
# running the same model twice WITHOUT a fixed seed returns slightly different
# numbers (Monte Carlo noise). Several routines below re-fit the SAME model --
# for example the State baseline is computed once for the band and again inside
# the validation module. We therefore derive the seed from the instance identity
# itself: the same level + region + target + year + tag always maps to the SAME
# seed, no matter where in the script it is called. Re-running the script then
# reproduces every figure and metric exactly.
instance_seed <- function(level, region, target_name, year_target, tag = "") {
  key <- paste(level, region, target_name, year_target, tag, sep = "|")
  as.integer(sum(utf8ToInt(key)) %% 2147483647L)
}

# =============================================================================
#  2. PACKAGE INITIALIZATION
# =============================================================================
message("Checking packages...")
cran_packages <- c("dplyr", "tidyr", "stringr", "lubridate", "ggplot2", "sf", "stringi")

if (isTRUE(config$vaccine$besag_enable)) cran_packages <- c(cran_packages, "spdep")

for (p in cran_packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

if (!requireNamespace("INLA", quietly = TRUE)) {
  message("Installing R-INLA (this may take a while)...")
  tryCatch(install.packages("INLA", repos = c(CRAN = "https://cloud.r-project.org",
                                              INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE),
           error = function(e) message("Failed to install INLA: ", conditionMessage(e)))
}
if (!requireNamespace("INLA", quietly = TRUE))
  stop("The INLA package is not available. Please install it manually.", call. = FALSE)

suppressPackageStartupMessages(library(INLA))
message("Packages OK.\n")

# =============================================================================
#  3. CORE FUNCTIONS (Data Cleaning & Preparation)
# =============================================================================

#' @title Standardize Municipality Names
#' @param x Character vector of municipality names.
#' @return Standardized character vector (uppercase, no accents).
standardize_muni <- function(x) {
  v <- as.character(x)
  v <- tryCatch(stringi::stri_trans_general(v, "Latin-ASCII"), error = function(e) v)
  trimws(toupper(v))
}

#' @title Clean String for Filenames
#' @param x Character string.
#' @return Cleaned string replacing non-alphanumeric chars with underscores.
clean_name <- function(x) {
  gsub("^_|_$", "", gsub("[^A-Z0-9]+", "_", standardize_muni(x)))
}

load_meningitis_base <- function(path) {
  if (!file.exists(path)) stop("File not found:\n  ", path, call. = FALSE)
  d <- readRDS(path)
  if (!is.data.frame(d)) stop("The RDS does not contain a data.frame.", call. = FALSE)
  message("  ", format(nrow(d), big.mark = ","), " rows loaded.")
  as.data.frame(d)
}

prepare_base <- function(d) {
  has_col <- function(x) x %in% names(d)
  source_epi <- config$epi_source %||% "auto"
  used_source <- NULL
  
  if (identical(source_epi, "data_sintoma")) {
    if (!has_col("data_sintoma")) stop("epi_source='data_sintoma' but col is missing.", call. = FALSE)
    d$data_sintoma <- as.Date(d$data_sintoma)
    d$ano_epi <- lubridate::epiyear(d$data_sintoma); d$sem_epi <- lubridate::epiweek(d$data_sintoma)
    used_source <- "data_sintoma"
  } else {
    if (has_col("anoepi") && has_col("sepi")) {
      d$ano_epi <- suppressWarnings(as.integer(as.character(d$anoepi)))
      d$sem_epi <- suppressWarnings(as.integer(as.character(d$sepi)))
      used_source <- "anoepi/sepi (base columns)"
    } else if (has_col("data_sintoma")) {
      d$data_sintoma <- as.Date(d$data_sintoma)
      d$ano_epi <- lubridate::epiyear(d$data_sintoma); d$sem_epi <- lubridate::epiweek(d$data_sintoma)
      used_source <- "data_sintoma [anoepi/sepi missing]"
    } else stop("Base lacks 'anoepi'/'sepi' and 'data_sintoma'.", call. = FALSE)
  }
  
  message("  Year/Week Source: ", used_source, ".")
  col_mun <- intersect(c("MUNI_RESID", "mun", "MUN"), names(d))[1]
  col_drs <- intersect(c("DRS_RESID",  "drs", "DRS"), names(d))[1]
  
  if (is.na(col_mun)) stop("Base lacks municipality column.", call. = FALSE)
  d$mun <- standardize_muni(d[[col_mun]])
  if (!is.na(col_drs)) d$drs <- standardize_muni(d[[col_drs]])
  
  for (cc in c("grupo_etiologico", "etiologia", "sorogrupo", "classi_fin", "evolucao")) {
    if (!has_col(cc)) d[[cc]] <- NA_character_
  }
  return(d)
}

filter_target <- function(d, target, apply_confirmed = TRUE) {
  if (is.character(target)) target <- config$targets[[target]]
  if (is.null(target)) target <- list(tipo = "todas")
  
  if (!isTRUE(target$tipo == "todas")) {
    if (!is.null(target$grupo_etiologico))
      d <- d[!is.na(d$grupo_etiologico) & d$grupo_etiologico %in% target$grupo_etiologico, ]
    if (!is.null(target$etiologia))
      d <- d[!is.na(d$etiologia) & d$etiologia %in% target$etiologia, ]
    if (!is.null(target$sorogrupo))
      d <- d[!is.na(d$sorogrupo) & d$sorogrupo %in% target$sorogrupo, ]
  }
  if (isTRUE(apply_confirmed)) {
    use_conf <- if (!is.null(target$confirmado)) isTRUE(target$confirmado) else isTRUE(config$only_confirmed)
    if (use_conf) d <- d[!is.na(d$classi_fin) & d$classi_fin == "Confirmado", ]
  }
  return(d)
}

list_regions <- function(level) {
  col <- switch(level, "DRS" = "drs", "Municipality" = "mun", "CIR" = "cir",
                stop("Invalid level: ", level))
  if (!col %in% names(base_menin)) stop("Column '", col, "' does not exist in base.")
  sort(unique(na.omit(base_menin[[col]])))
}

aggregate_weekly <- function(base, level, region, target) {
  if (level == "State") {
    d <- base
  } else {
    col <- switch(level, "DRS" = "drs", "Municipality" = "mun", "CIR" = "cir")
    d   <- base[!is.na(base[[col]]) & base[[col]] == region, ]
  }
  d <- filter_target(d, target)
  d %>% filter(!is.na(ano_epi), !is.na(sem_epi)) %>% filter(sem_epi != 53) %>%
    group_by(ano_epi, sem_epi) %>% summarise(cases = n(), .groups = "drop") %>%
    rename(year = ano_epi, week = sem_epi)
}

complete_weeks <- function(df, years) {
  df <- df[df$year %in% years, , drop = FALSE]
  grid <- expand.grid(year = years, week = 1:52)
  merge(grid, df, by = c("year", "week"), all.x = TRUE) %>%
    mutate(cases = ifelse(is.na(cases), 0L, cases)) %>% arrange(year, week)
}

add_harmonics <- function(df, nh) {
  for (k in seq_len(nh)) {
    df[[paste0("sin", k)]] <- sin(2 * pi * k * df$week / 52)
    df[[paste0("cos", k)]] <- cos(2 * pi * k * df$week / 52)
  }
  return(df)
}

build_formula <- function(model) {
  if (!is.null(model$formula)) return(model$formula)
  seasonal <- model$seasonal %||% "rw1"
  year_term <- "f(ano_idx, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  
  if (identical(seasonal, "harmonicos")) {
    nh <- model$n_harmonics %||% 2
    harm <- paste(sprintf("sin%d + cos%d", seq_len(nh), seq_len(nh)), collapse = " + ")
    stats::as.formula(paste("cases ~ 1 +", harm, "+", year_term))
  } else if (identical(seasonal, "nenhum")) {
    stats::as.formula(paste("cases ~ 1 +", year_term))
  } else {
    stats::as.formula(paste("cases ~ 1 +",
                            "f(week, model = 'rw1', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
                            year_term))
  }
}

get_nh_model <- function(model) {
  if (identical(model$seasonal %||% "rw1", "harmonicos")) (model$n_harmonics %||% 2) else 0L
}

simulate_cumulative <- function(samples, target_rows, idx_zip = NA_integer_) {
  vapply(samples, function(s) {
    mu <- exp(s$latent[target_rows, 1]); y <- rpois(length(mu), mu)
    if (!is.na(idx_zip)) {
      p <- as.numeric(s$hyperpar[idx_zip])
      if (is.finite(p) && p > 0) y[runif(length(y)) < p] <- 0L
    }
    cumsum(y)
  }, numeric(length(target_rows)))
}

get_zip_hyper_idx <- function(samples, family) {
  if (!grepl("zeroinflated", family, ignore.case = TRUE)) return(NA_integer_)
  names_hyper <- names(samples[[1]]$hyperpar)
  cand <- grep("zero|prob", names_hyper, ignore.case = TRUE)
  if (length(cand) >= 1) cand[1] else NA_integer_
}

plot_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 7.5), legend.key.size = unit(0.8, "lines"),
          plot.title = element_text(size = base_size + 3, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40"))
}

save_band_plot <- function(plt, filename, region_folder) {
  if (isTRUE(!dir.exists(region_folder))) dir.create(region_folder, recursive = TRUE)
  path <- file.path(region_folder, paste0(filename, ".png"))
  ggsave(path, plt, width = config$width / config$dpi, height = config$height / config$dpi,
         units = "in", dpi = config$dpi, device = "png", bg = "white")
  message("    -> Saved: ", path); invisible(path)
}

# ---- Base Loading Execution -------------------------------------------------
message("Loading meningitis base...")
base_menin <- load_meningitis_base(config$data_path)
base_menin <- prepare_base(base_menin)
message("Base prepared.\n")

# =============================================================================
#  3.5 MUNICIPALITY -> CIR SPATIAL MAPPING
# =============================================================================
get_municipalities_sp <- function() {
  if (is.null(config$shapefile_municipalities) || !nzchar(config$shapefile_municipalities))
    stop("Define config$shapefile_municipalities to map CIR.", call. = FALSE)
  
  mp  <- sf::st_read(config$shapefile_municipalities, quiet = TRUE)
  idc <- config$col_id_municipalities
  
  if (is.null(idc)) {
    gcol  <- attr(mp, "sf_column")
    cands <- setdiff(grep("nm_mun|name_muni|nome|municip|^mun", names(mp),
                          ignore.case = TRUE, value = TRUE), gcol)
    if (!length(cands)) stop("Could not detect municipality name column.", call. = FALSE)
    idc <- cands[1]; message("  [Muni] Name column detected: '", idc, "'.")
  }
  if (!idc %in% names(mp)) stop("Column '", idc, "' missing in municipality shapefile.", call. = FALSE)
  
  mp$chave_mun <- standardize_muni(mp[[idc]])
  mp["chave_mun"]
}

prepare_cir_mapping <- function(base) {
  if ("cir" %in% names(base)) return(base)
  if (is.null(config$shapefile_cir) || !file.exists(config$shapefile_cir))
    stop("CIR Shapefile not found: ", config$shapefile_cir, call. = FALSE)
  
  suppressWarnings(try(sf::sf_use_s2(FALSE), silent = TRUE))
  cir_sf <- sf::st_make_valid(sf::st_read(config$shapefile_cir, quiet = TRUE))
  gcol <- attr(cir_sf, "sf_column"); idc <- config$col_id_cir
  
  if (is.null(idc)) {
    cands <- setdiff(grep("cir|nome|nom|name|regi", names(cir_sf), ignore.case = TRUE, value = TRUE), gcol)
    idc <- if (length(cands)) cands[1] else setdiff(names(cir_sf), gcol)[1]
    message("  [CIR] Identifier column detected: '", idc, "'.")
  }
  
  cir_sf[[idc]] <- standardize_muni(cir_sf[[idc]])
  munpoly <- sf::st_make_valid(get_municipalities_sp())
  cir_sf  <- sf::st_transform(cir_sf, sf::st_crs(munpoly))
  
  cent <- suppressWarnings(sf::st_point_on_surface(munpoly))
  jn   <- suppressWarnings(sf::st_join(cent, cir_sf[idc], join = sf::st_within))
  
  mapa <- data.frame(chave_mun = as.character(jn$chave_mun),
                     cir = as.character(jn[[idc]]), stringsAsFactors = FALSE)
  mapa <- mapa[!is.na(mapa$cir) & !duplicated(mapa$chave_mun), ]
  
  n_miss <- sum(is.na(match(unique(base$mun), mapa$chave_mun)))
  if (n_miss > 0) message("  [CIR] WARNING: ", n_miss, " municipalities without CIR (NA).")
  
  base$cir <- mapa$cir[match(base$mun, mapa$chave_mun)]
  message("  [CIR] ", length(unique(na.omit(base$cir))), " CIRs assigned.")
  return(base)
}

needs_cir <- function() {
  cfg <- config$vaccine
  vac_cir <- "CIR" %in% cfg$levels
  bes_cir <- isTRUE(cfg$besag_enable) && ("CIR" %in% (cfg$besag_levels %||% character(0)))
  vac_cir || bes_cir
}

if (needs_cir()) {
  message("Preparing CIR mapping on base (Municipality -> CIR)...")
  base_menin <- tryCatch(prepare_cir_mapping(base_menin),
                         error = function(e) {
                           message("  [CIR] FAILED: ", e$message, "\n  -> CIR will be SKIPPED."); base_menin
                         })
  if (!"cir" %in% names(base_menin)) {
    config$vaccine$levels <- setdiff(config$vaccine$levels, "CIR")
    config$vaccine$besag_levels <- setdiff(config$vaccine$besag_levels, "CIR")
  }
}

# =============================================================================
#  4. INLA ENGINE
# =============================================================================

matrix_zones <- function(acum_mat) {
  qs <- config$inla_percentiles %||% c(lower = 0.025, media = 0.50, upper = 0.975)
  data.frame(week = 1:nrow(acum_mat),
             lower = apply(acum_mat, 1, quantile, qs[["lower"]], na.rm = TRUE),
             media = apply(acum_mat, 1, quantile, qs[["media"]], na.rm = TRUE),
             upper = apply(acum_mat, 1, quantile, qs[["upper"]], na.rm = TRUE))
}

run_inla_acum_mat <- function(df_train, n_samples = 1000, min_cases = 30, model = NULL, seed = NULL) {
  if (sum(df_train$cases) < min_cases) { message("    [WARNING] Cases < min; skipping."); return(NULL) }
  unique_years <- sort(unique(df_train$year))
  if (length(unique_years) < 3) { message("    [WARNING] < 3 years; skipping."); return(NULL) }
  
  df_train$ano_idx <- match(df_train$year, unique_years)
  if (is.null(model)) model <- list()
  
  family <- model$family %||% "poisson"
  formula_inla <- build_formula(model)
  nh <- get_nh_model(model)
  
  idx_new_year <- length(unique_years) + 1L
  df_pred <- data.frame(cases = NA_integer_, week = 1:52, ano_idx = idx_new_year)
  df_full <- rbind(df_train[, c("cases", "week", "ano_idx")], df_pred)
  
  if (nh > 0) df_full <- add_harmonics(df_full, nh)
  
  fit <- tryCatch(
    inla(formula_inla, family = family, data = df_full,
         control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
         control.predictor = list(link = 1, compute = TRUE), verbose = FALSE),
    error = function(e) { message("    [INLA ERROR] ", e$message); NULL }
  )
  
  if (is.null(fit) || is.null(fit$misc$configs)) {
    message("    [WARNING] INLA configs matrix did not converge. Skipping.")
    return(NULL)
  }
  
  idx_pred <- (nrow(df_train) + 1):nrow(df_full)
  if (!is.null(seed)) set.seed(seed)
  samples <- tryCatch(inla.posterior.sample(n_samples, fit),
                      error = function(e) { message("    [INLA Sample ERROR] ", e$message); NULL })
  if (is.null(samples)) return(NULL)
  
  pred_rows <- grep("^Predictor", rownames(samples[[1]]$latent))
  acum_mat <- simulate_cumulative(samples, pred_rows[idx_pred], get_zip_hyper_idx(samples, family))
  
  list(acum_mat = acum_mat, canal = matrix_zones(acum_mat),
       dic = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
       waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_)
}

run_inla_acum_mat_vac <- function(df_train, cov_by_year, cov_target_year,
                                  n_samples = 1000, min_cases = 30, model = NULL,
                                  standardize = TRUE, seed = NULL) {
  if (sum(df_train$cases) < min_cases) { message("    [WARNING] Cases < min; skipping."); return(NULL) }
  unique_years <- sort(unique(df_train$year))
  if (length(unique_years) < 3) { message("    [WARNING] < 3 years; skipping."); return(NULL) }
  
  df_train$ano_idx <- match(df_train$year, unique_years)
  if (is.null(model)) model <- list()
  
  family <- model$family %||% "poisson"; nh <- get_nh_model(model)
  df_train$year <- suppressWarnings(as.integer(as.character(df_train$year)))
  
  if (!is.null(cov_by_year) && isTRUE(nrow(cov_by_year) > 0)) {
    cov_by_year$year <- suppressWarnings(as.integer(as.character(cov_by_year$ano)))
    df_train <- merge(df_train, cov_by_year, by = "year", all.x = TRUE)
  }
  
  if (!"cov_primaria" %in% names(df_train)) df_train$cov_primaria <- NA_real_
  if (!"cov_reforco" %in% names(df_train)) df_train$cov_reforco <- NA_real_
  
  calc_safe <- function(col, func) { if (isTRUE(any(!is.na(col)))) func(col, na.rm = TRUE) else NA_real_ }
  
  mp <- c(cov_primaria = calc_safe(df_train$cov_primaria, mean),
          cov_reforco  = calc_safe(df_train$cov_reforco, mean))
  sp <- c(cov_primaria = calc_safe(df_train$cov_primaria, stats::sd),
          cov_reforco  = calc_safe(df_train$cov_reforco, stats::sd))
  
  usa_prim <- isTRUE(any(is.finite(df_train$cov_primaria))) && isTRUE(is.finite(sp["cov_primaria"])) && isTRUE(sp["cov_primaria"] > 0)
  usa_ref  <- isTRUE(any(is.finite(df_train$cov_reforco)))  && isTRUE(is.finite(sp["cov_reforco"]))  && isTRUE(sp["cov_reforco"] > 0)
  
  z <- function(x, nm) {
    if (isTRUE(standardize) && isTRUE(is.finite(sp[nm])) && isTRUE(sp[nm] > 0)) {
      (x - mp[nm]) / sp[nm]
    } else { x }
  }
  
  imp <- function(x, nm) { if(isTRUE(length(x) > 0)) { x[is.na(x)] <- mp[nm] }; x }
  
  df_train$zp <- z(imp(df_train$cov_primaria, "cov_primaria"), "cov_primaria")
  df_train$zr <- z(imp(df_train$cov_reforco,  "cov_reforco"),  "cov_reforco")
  
  idx_new_year <- length(unique_years) + 1L
  df_pred <- data.frame(cases = NA_integer_, week = 1:52, ano_idx = idx_new_year)
  
  esc1 <- function(x, nm) {
    x <- suppressWarnings(as.numeric(x))
    if (isTRUE(length(x) >= 1) && isTRUE(is.finite(x[1]))) x[1] else mp[[nm]]
  }
  cp <- esc1(cov_target_year$cov_primaria, "cov_primaria")
  cr <- esc1(cov_target_year$cov_reforco,  "cov_reforco")
  
  df_pred$zp <- z(cp, "cov_primaria")
  df_pred$zr <- z(cr, "cov_reforco")
  
  cols_keep <- c("cases", "week", "ano_idx", "zp", "zr")
  df_full <- rbind(df_train[, cols_keep], df_pred[, cols_keep])
  if (nh > 0) df_full <- add_harmonics(df_full, nh)
  
  base_form  <- build_formula(model)
  termos_vac <- c(if (usa_prim) "zp", if (usa_ref) "zr")
  
  if (length(termos_vac) > 0) {
    rhs <- paste(deparse(base_form[[3]]), collapse = " ")
    formula_inla <- stats::as.formula(paste("cases ~", rhs, "+", paste(termos_vac, collapse = " + ")))
  } else {
    formula_inla <- base_form
  }
  
  fit <- tryCatch(
    inla(formula_inla, family = family, data = df_full,
         control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
         control.predictor = list(link = 1, compute = TRUE),
         control.fixed = list(prec = list(zp = 4, zr = 4)), verbose = FALSE),
    error = function(e) { message("    [INLA vac ERROR] ", e$message); NULL }
  )
  
  if (is.null(fit) || is.null(fit$misc$configs)) {
    message("    [WARNING] INLA vac configs did not converge. Skipping.")
    return(NULL)
  }
  
  idx_pred <- (nrow(df_train) + 1):nrow(df_full)
  if (!is.null(seed)) set.seed(seed)
  samples <- tryCatch(inla.posterior.sample(n_samples, fit),
                      error = function(e) { message("    [INLA Sample vac ERROR] ", e$message); NULL })
  if (is.null(samples)) return(NULL)
  
  pred_rows <- grep("^Predictor", rownames(samples[[1]]$latent))
  acum_mat <- simulate_cumulative(samples, pred_rows[idx_pred], get_zip_hyper_idx(samples, family))
  coefs <- tryCatch(fit$summary.fixed, error = function(e) NULL)
  
  # Diagnostic bundle: everything log_vaccine_diagnostics() needs to detect
  # extrapolation in covariate space (target-year z-score compared with the range
  # seen during training) and to report the vaccine coefficient's uncertainty.
  diag <- list(
    mp = mp, sp = sp,
    z_train_range_primaria = if (isTRUE(any(is.finite(df_train$zp)))) range(df_train$zp, na.rm = TRUE) else c(NA_real_, NA_real_),
    z_train_range_reforco  = if (isTRUE(any(is.finite(df_train$zr)))) range(df_train$zr, na.rm = TRUE) else c(NA_real_, NA_real_),
    z_target_primaria = df_pred$zp[1], z_target_reforco = df_pred$zr[1],
    cov_target_primaria_raw = cp, cov_target_reforco_raw = cr,
    coef_primaria_sd = if (isTRUE(usa_prim) && isTRUE("zp" %in% rownames(coefs))) coefs["zp","sd"] else NA_real_,
    coef_reforco_sd  = if (isTRUE(usa_ref)  && isTRUE("zr" %in% rownames(coefs))) coefs["zr","sd"] else NA_real_
  )
  
  list(acum_mat = acum_mat, canal = matrix_zones(acum_mat),
       dic = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
       waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_,
       coef_primaria = if (isTRUE(usa_prim) && isTRUE("zp" %in% rownames(coefs))) coefs["zp","mean"] else NA_real_,
       coef_reforco  = if (isTRUE(usa_ref)  && isTRUE("zr" %in% rownames(coefs))) coefs["zr","mean"] else NA_real_,
       diag = diag)
}

# =============================================================================
#  5. VACCINE DATA PROCESSING
# =============================================================================

read_vaccine_csv <- function(path) {
  if (!file.exists(path)) { message("  [Vaccine] File not found: ", path); return(NULL) }
  l1  <- readLines(path, n = 1, warn = FALSE, encoding = "latin1")
  sep <- if (lengths(regmatches(l1, gregexpr(";", l1))) >= lengths(regmatches(l1, gregexpr(",", l1)))) ";" else ","
  
  df <- tryCatch(utils::read.table(path, sep = sep, header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE, fileEncoding = "latin1",
                                   colClasses = "character", quote = "\""),
                 error = function(e) { message("  [Vaccine] Error reading ", path, ": ", e$message); NULL })
  if (is.null(df)) return(NULL)
  
  names(df) <- trimws(names(df)); col_muni <- names(df)[1]
  anos_cols <- grep("^[0-9]{4}$", names(df), value = TRUE)
  if (!length(anos_cols)) { message("  [Vaccine] No year columns in ", path); return(NULL) }
  
  muni_raw <- trimws(df[[col_muni]])
  manter <- !is.na(muni_raw) & nzchar(muni_raw) & !grepl("^total", muni_raw, ignore.case = TRUE)
  df <- df[manter, , drop = FALSE]; muni_raw <- muni_raw[manter]
  
  cod <- sub("^\\s*([0-9]{6,7}).*$", "\\1", muni_raw); cod[cod == muni_raw] <- NA_character_
  nome <- trimws(sub("^\\s*[0-9]{6,7}\\s*", "", muni_raw))
  code6 <- ifelse(is.na(cod), NA_character_, substr(cod, 1, 6))
  
  to_num <- function(x) { 
    x <- trimws(x); x[x %in% c("", "-", "NA", "ND")] <- NA
    x <- gsub("\\.", "", x); x <- gsub(",", ".", x)
    suppressWarnings(as.numeric(x)) 
  }
  
  long <- do.call(rbind, lapply(anos_cols, function(a)
    data.frame(code_muni_6 = code6, mun = standardize_muni(nome), ano = as.integer(a),
               cobertura = to_num(df[[a]]), stringsAsFactors = FALSE)))
  
  long[!is.na(long$mun) & long$mun != "", , drop = FALSE]
}

load_vaccines <- function() {
  cfg <- config$vaccine; out <- list()
  for (rot in names(cfg$files)) {
    path <- file.path(cfg$folder, cfg$files[[rot]])
    message("  [Vaccine] Reading ", rot, ": ", path)
    d <- read_vaccine_csv(path)
    if (!is.null(d)) { names(d)[names(d) == "cobertura"] <- rot; out[[rot]] <- d }
  }
  out
}

map_muni_region <- function(level) {
  if (level == "State") return(NULL)
  col <- switch(level, "DRS" = "drs", "CIR" = "cir", "Municipality" = "mun", stop("Invalid level: ", level))
  if (!col %in% names(base_menin)) stop("Column '", col, "' missing for aggregation in ", level, ".")
  unique(data.frame(mun = base_menin$mun, regiao = base_menin[[col]], stringsAsFactors = FALSE))
}

aggregate_vaccine_level <- function(vac_long, col_val, level) {
  cfg <- config$vaccine; v <- vac_long; v$val <- v[[col_val]]
  v <- v[!is.na(v$val), , drop = FALSE]
  if (!nrow(v)) return(data.frame(regiao = character(0), ano = integer(0)))
  
  if (level == "Municipality") {
    agg <- aggregate(val ~ mun + ano, data = v, FUN = mean, na.rm = TRUE)
    names(agg) <- c("regiao", "ano", col_val); return(agg)
  }
  
  if (level == "State") {
    if (identical(cfg$aggregation, "weighted_mean") && !is.null(cfg$target_pop)) {
      v <- merge(v, cfg$target_pop, by = c("code_muni_6", "ano"), all.x = TRUE)
      agg <- do.call(rbind, lapply(split(v, v$ano), function(d) {
        w <- d$pop; if (all(is.na(w))) w <- rep(1, nrow(d))
        data.frame(regiao = "SP", ano = d$ano[1], val = stats::weighted.mean(d$val, w, na.rm = TRUE)) 
      }))
    } else { 
      agg <- aggregate(val ~ ano, data = v, FUN = mean, na.rm = TRUE)
      agg$regiao <- "SP"; agg <- agg[, c("regiao", "ano", "val")] 
    }
    names(agg)[names(agg) == "val"] <- col_val; return(agg)
  }
  
  mm <- map_muni_region(level); v2 <- merge(v, mm, by = "mun", all.x = TRUE)
  v2 <- v2[!is.na(v2$regiao), , drop = FALSE]
  if (!nrow(v2)) return(data.frame(regiao = character(0), ano = integer(0)))
  
  if (identical(cfg$aggregation, "weighted_mean") && !is.null(cfg$target_pop)) {
    v2 <- merge(v2, cfg$target_pop, by = c("code_muni_6", "ano"), all.x = TRUE)
    agg <- do.call(rbind, lapply(split(v2, list(v2$regiao, v2$ano), drop = TRUE), function(d) {
      w <- d$pop; if (all(is.na(w))) w <- rep(1, nrow(d))
      data.frame(regiao = d$regiao[1], ano = d$ano[1], val = stats::weighted.mean(d$val, w, na.rm = TRUE)) 
    }))
  } else {
    agg <- aggregate(val ~ regiao + ano, data = v2, FUN = mean, na.rm = TRUE)
  }
  
  names(agg)[names(agg) == "val"] <- col_val; return(agg)
}

covariate_target <- function(vacs_long, target_name, level, region) {
  spec <- config$vaccine$cov_by_target[[target_name]]; if (is.null(spec)) return(NULL)
  
  build_group <- function(labels) {
    if (isTRUE(is.null(labels)) || isTRUE(length(labels) == 0)) return(NULL)
    parts <- list()
    for (rot in labels) {
      vl <- vacs_long[[rot]]; if (isTRUE(is.null(vl))) next
      ag <- aggregate_vaccine_level(vl, rot, level)
      if (level != "State") ag <- ag[ag$regiao == region, , drop = FALSE]
      if (isTRUE(is.null(ag)) || isTRUE(nrow(ag) == 0)) next
      parts[[rot]] <- ag[, c("ano", rot)]
    }
    if (isTRUE(length(parts) == 0)) return(NULL)
    m <- Reduce(function(a, b) merge(a, b, by = "ano", all = TRUE), parts)
    valcols <- setdiff(names(m), "ano")
    data.frame(ano = m$ano, val = rowMeans(m[, valcols, drop = FALSE], na.rm = TRUE))
  }
  
  prim <- build_group(spec$primary); ref <- build_group(spec$booster)
  if (is.null(prim) && is.null(ref)) return(NULL)
  
  anos <- sort(unique(c(prim$ano, ref$ano))); out <- data.frame(ano = anos)
  out$cov_primaria <- if (!is.null(prim)) prim$val[match(anos, prim$ano)] else NA_real_
  out$cov_reforco  <- if (!is.null(ref))  ref$val[match(anos, ref$ano)]   else NA_real_
  out
}

# Shifts the coverage table forward by 'lag' years: a row originally labelled
# year Y now applies to year Y+lag. With lag = 1, coverage MEASURED in Y is used
# to explain/predict CASES in Y+1. Shifting the table ONCE here propagates the
# lag to BOTH the training merge (inside run_inla_acum_mat_vac) and the
# target-year lookup (get_year_covariate), because both read this same shifted
# table -- so the two can never fall out of sync. lag = 0 restores the
# contemporaneous design (same-year coverage explaining same-year cases).
shift_coverage_years <- function(cov_tab, lag = 0L) {
  if (is.null(cov_tab) || !isTRUE(nrow(cov_tab) > 0) || identical(as.integer(lag), 0L)) return(cov_tab)
  cov_tab$ano <- cov_tab$ano + as.integer(lag)
  cov_tab
}

get_year_covariate <- function(cov_tab, year_target) {
  if (is.null(cov_tab) || !isTRUE(nrow(cov_tab) > 0)) {
    return(list(cov_primaria = NA_real_, cov_reforco = NA_real_))
  }
  anos_tab <- suppressWarnings(as.integer(as.character(cov_tab$ano)))
  alvo <- suppressWarnings(as.integer(year_target))
  row_data <- cov_tab[!is.na(anos_tab) & anos_tab == alvo, , drop = FALSE]
  
  fetch <- function(col) {
    if (isTRUE(nrow(row_data) > 0) && col %in% names(row_data) && length(row_data[[col]]) >= 1) {
      val <- row_data[[col]][1]
      if (isTRUE(is.finite(val))) return(val)
    }
    NA_real_
  }
  list(cov_primaria = fetch("cov_primaria"), cov_reforco = fetch("cov_reforco"))
}

# =============================================================================
#  5.5 VACCINE DIAGNOSTIC LOG -- diagnosing abnormally wide uncertainty bands
# =============================================================================
# Accumulates log lines in memory and writes a .txt at the end of the routine.
# WHY THIS EXISTS: a Bayesian band can come out with enormous uncertainty for a
# reason that is invisible in the final figure. The most common cause here is
# EXTRAPOLATION in covariate space: if a vaccine's historical coverage was very
# stable, its training standard deviation is small, so even a moderate real-world
# drop in the target year becomes an extreme z-score far outside anything the
# model saw. Multiplying an uncertain coefficient by an extreme z-score and then
# exponentiating (Poisson log link) makes the interval explode. This log reports
# exactly the numbers needed to confirm or rule out that mechanism.
.VACLOG <- new.env()
.VACLOG$linhas <- character(0)

vaclog <- function(...) {
  txt <- paste0(...)
  .VACLOG$linhas <- c(.VACLOG$linhas, txt)
  invisible(txt)
}

vaclog_flush <- function(output_folder, filename) {
  if (!length(.VACLOG$linhas)) return(invisible(NULL))
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  path <- file.path(output_folder, filename)
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines(.VACLOG$linhas, con)
  close(con)
  message("\n  [Vaccine Diagnostics] Log saved to:\n    -> ", path)
  invisible(path)
}

.fmt <- function(x, d = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (isTRUE(length(x) == 0) || isTRUE(!is.finite(x))) return("NA")
  formatC(x, format = "f", digits = d)
}

log_vaccine_diagnostics <- function(level, region, target_name, year_target, sem, com) {
  if (!isTRUE(config$vaccine$diagnostic_log %||% TRUE)) return(invisible(NULL))
  vaclog(strrep("=", 90))
  vaclog(sprintf("INSTANCE: %s | %s | %s | Year %s", level, region, target_name, as.character(year_target)))
  vaclog(strrep("-", 90))
  
  if (is.null(com) || is.null(com$diag)) {
    vaclog("  [!] WITH-vaccine model produced no diagnostic info (model may have failed/skipped).")
  } else {
    d <- com$diag
    vaclog(sprintf("  Coverage lag applied: %s year(s)", as.character(config$vaccine$coverage_lag_years %||% 0L)))
    vaclog(sprintf("  PRIMARY coverage -> train mean=%s sd=%s | target(raw)=%s | target z=%s | train z range=[%s, %s]",
                   .fmt(d$mp[["cov_primaria"]]), .fmt(d$sp[["cov_primaria"]]), .fmt(d$cov_target_primaria_raw),
                   .fmt(d$z_target_primaria), .fmt(d$z_train_range_primaria[1]), .fmt(d$z_train_range_primaria[2])))
    vaclog(sprintf("  BOOSTER coverage -> train mean=%s sd=%s | target(raw)=%s | target z=%s | train z range=[%s, %s]",
                   .fmt(d$mp[["cov_reforco"]]), .fmt(d$sp[["cov_reforco"]]), .fmt(d$cov_target_reforco_raw),
                   .fmt(d$z_target_reforco), .fmt(d$z_train_range_reforco[1]), .fmt(d$z_train_range_reforco[2])))
    
    zp <- suppressWarnings(as.numeric(d$z_target_primaria)); rp <- d$z_train_range_primaria
    zr <- suppressWarnings(as.numeric(d$z_target_reforco));  rr <- d$z_train_range_reforco
    extra_p <- isTRUE(is.finite(zp)) && isTRUE(all(is.finite(rp))) && (zp < rp[1] || zp > rp[2])
    extra_r <- isTRUE(is.finite(zr)) && isTRUE(all(is.finite(rr))) && (zr < rr[1] || zr > rr[2])
    if (extra_p) vaclog("  [WARNING] target-year PRIMARY coverage z-score is OUTSIDE the training range -> EXTRAPOLATION.")
    if (extra_r) vaclog("  [WARNING] target-year BOOSTER  coverage z-score is OUTSIDE the training range -> EXTRAPOLATION.")
    
    vaclog(sprintf("  Coefficient (zp, primary) posterior: mean=%s  sd=%s",
                   .fmt(com$coef_primaria, 4), .fmt(d$coef_primaria_sd, 4)))
    vaclog(sprintf("  Coefficient (zr, booster) posterior: mean=%s  sd=%s",
                   .fmt(com$coef_reforco, 4), .fmt(d$coef_reforco_sd, 4)))
  }
  
  if (!is.null(sem) && !is.null(com)) {
    w_sem <- suppressWarnings(sem$canal$upper[52] - sem$canal$lower[52])
    w_com <- suppressWarnings(com$canal$upper[52] - com$canal$lower[52])
    ratio <- if (isTRUE(is.finite(w_sem)) && isTRUE(w_sem > 0)) w_com / w_sem else NA_real_
    vaclog(sprintf("  Band width at EW52 -> NO vaccine=%s | WITH vaccine=%s | ratio(WITH/NO)=%s",
                   .fmt(w_sem, 1), .fmt(w_com, 1), .fmt(ratio, 2)))
    if (isTRUE(is.finite(ratio)) && ratio > 3)
      vaclog("  [WARNING] WITH-vaccine band is >3x WIDER than NO-vaccine band at EW52 -- see extrapolation flags above.")
  } else {
    vaclog("  (band-width comparison skipped: sem or com model missing)")
  }
  vaclog("")
}

# =============================================================================
#  6. COMPARISON (WITH / WITHOUT VACCINE)
# =============================================================================

plot_predicted_vs_observed <- function(series, title, year_target, level, rot_reg, rot_tgt) {
  df_lin <- rbind(
    data.frame(week = series$week, value = series$pred_sem, type = "Predicted NO vaccine"),
    data.frame(week = series$week, value = series$pred_com, type = "Predicted WITH vaccine"))
  df_obs <- data.frame(week = series$week, value = series$obs_acum)
  
  cores  <- c("Predicted NO vaccine" = "#1565C0", "Predicted WITH vaccine" = "#C62828")
  p <- ggplot() +
    geom_ribbon(data = series, aes(x = week, ymin = li_sem, ymax = ls_sem), fill = "#1565C0", alpha = 0.12) +
    geom_ribbon(data = series, aes(x = week, ymin = li_com, ymax = ls_com), fill = "#C62828", alpha = 0.12) +
    geom_line(data = df_lin, aes(x = week, y = value, color = type), linewidth = 0.8) +
    geom_point(data = df_obs, aes(x = week, y = value), color = "black", size = 1.1) +
    geom_line(data = df_obs, aes(x = week, y = value), color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_color_manual(values = cores, name = NULL) +
    scale_x_continuous(breaks = seq(1, 52, 4), labels = function(x) paste0("EW", sprintf("%02d", x))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(title = paste0(title, " — Predicted vs Observed"),
         subtitle = "Black = observed; Blue = predicted NO vaccine; Red = predicted WITH vaccine",
         x = "Epidemiological week", y = "Cumulative cases in the year") +
    plot_theme()
  
  folder <- file.path(config$output_folder, as.character(year_target), level, "predito_vs_observado")
  save_band_plot(p, paste0("predobs_", tolower(level), "_", rot_reg, "_", rot_tgt, "_", year_target), folder)
}

plot_simple_zones <- function(canal, obs_acum, title, year_target) {
  ymax <- suppressWarnings(max(c(canal$upper, if (isTRUE(nrow(obs_acum) > 0)) obs_acum$acum else 0), na.rm = TRUE)) * 1.10
  if (!isTRUE(is.finite(ymax)) || isTRUE(ymax <= 0)) ymax <- suppressWarnings(max(c(canal$upper, 1), na.rm = TRUE))
  canal$teto <- if (isTRUE(is.finite(ymax))) max(ymax, 1) else 1
  
  lab_obs    <- paste0("Cumulative observed ", year_target)
  lab_centro <- "Median"
  
  p <- ggplot(canal, aes(x = week)) +
    geom_ribbon(aes(ymin = upper, ymax = teto,  fill = "Epidemic")) +
    geom_ribbon(aes(ymin = media, ymax = upper, fill = "Alert")) +
    geom_ribbon(aes(ymin = lower, ymax = media, fill = "Expected")) +
    geom_ribbon(aes(ymin = 0,     ymax = lower, fill = "Below expected")) +
    geom_line(aes(y = lower), color = "black", linewidth = 0.3) +
    geom_line(aes(y = upper), color = "black", linewidth = 0.3) +
    geom_line(aes(y = media, color = lab_centro), linewidth = 0.7)
  
  if (isTRUE(!is.null(obs_acum)) && isTRUE(nrow(obs_acum) > 0)) {
    p <- p +
      geom_line(data = obs_acum, aes(x = week, y = acum, color = lab_obs), linewidth = 0.8) +
      geom_point(data = obs_acum, aes(x = week, y = acum, color = lab_obs), size = 1.6)
  }
  
  cores_lin <- setNames("#1565C0", lab_centro); cores_lin[lab_obs] <- "#000000"
  
  p +
    scale_fill_manual(values = c("Below expected" = "#0bc82c", "Expected" = "#edea44",
                                 "Alert" = "#fe3d3d", "Epidemic" = "#ba0606"),
                      breaks = c("Below expected","Expected","Alert","Epidemic"), name = NULL) +
    scale_color_manual(values = cores_lin, name = NULL) +
    scale_x_continuous(breaks = seq(1, 52, 4), labels = function(x) paste0("EW", sprintf("%02d", x))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    guides(fill = guide_legend(order = 1, nrow = 1), color = guide_legend(order = 2, nrow = 1)) +
    labs(title = title, subtitle = paste0("Band ", year_target),
         x = "Epidemiological week", y = "Cumulative cases in the year") +
    plot_theme()
}
compare_vaccine_instance <- function(level, region, target_name, year_target, vacs_long) {
  spec <- config$vaccine$cov_by_target[[target_name]]; if (is.null(spec)) return(NULL)
  
  rot_reg <- if (level == "State") "ESP" else clean_name(region)
  rot_tgt <- clean_name(target_name)
  message("  [Vac cmp] ", level, ":", region, " | ", target_name, " | Year ", year_target)
  
  df_all <- aggregate_weekly(base_menin, level, if (level == "State") NULL else region, target_name)
  train_years <- setdiff(config$train_years, year_target)
  df_train <- df_all %>% filter(year %in% train_years) %>% complete_weeks(train_years)
  
  if (isTRUE(sum(df_train$cases, na.rm = TRUE) < config$min_train_cases)) { message("    Few cases; skipping."); return(NULL) }
  
  cov_tab <- covariate_target(vacs_long, target_name, level, region)
  if (is.null(cov_tab)) { message("    No coverage; skipping."); return(NULL) }
  # Apply the configured lag ONCE here; both the training merge and the
  # target-year lookup below read this same shifted table.
  cov_tab <- shift_coverage_years(cov_tab, config$vaccine$coverage_lag_years %||% 0L)
  
  cov_by_year <- cov_tab[, c("ano", "cov_primaria", "cov_reforco")]
  cov_year_target <- get_year_covariate(cov_tab, year_target)
  
  if (isTRUE(all(is.na(unlist(cov_year_target))))) {
    message("    No coverage in target year ", year_target, "; skipping.") ; return(NULL) 
  }
  
  obs <- df_all %>% filter(year == year_target) %>% select(week, cases)
  obs <- data.frame(week = 1:52) %>% left_join(obs, by = "week") %>%
    mutate(cases = tidyr::replace_na(cases, 0L)) %>% arrange(week) %>% mutate(acum = cumsum(cases))
  y_final <- suppressWarnings(max(obs$acum, na.rm = TRUE))
  
  sem <- run_inla_acum_mat(df_train, config$n_samples, config$min_train_cases, config$inla_models[[level]],
                            seed = instance_seed(level, region, target_name, year_target, "sem"))
  com <- run_inla_acum_mat_vac(df_train, cov_by_year, cov_year_target, config$n_samples,
                               config$min_train_cases, config$inla_models[[level]], config$vaccine$standardize,
                               seed = instance_seed(level, region, target_name, year_target, "com"))
  
  # Diagnostic log -- deliberately runs even if one of the two models failed, so
  # that the log records WHY it failed. See section 5.5.
  tryCatch(log_vaccine_diagnostics(level, region, target_name, year_target, sem, com),
           error = function(e) message("    [Vaccine Diagnostics ERROR] ", e$message))
  
  if (is.null(sem) || is.null(com)) { message("    A model failed; skipping."); return(NULL) }
  
  title <- paste0(target_name, " — Band ", year_target)
  
  folder <- file.path(config$output_folder, as.character(year_target), level)
  
  save_band_plot(plot_simple_zones(sem$canal, obs, paste0(title, " (NO Vaccine)"), year_target),
                 paste0("band_", tolower(level), "_", rot_reg, "_", rot_tgt, "_", year_target, "_NOvac"), folder)
  save_band_plot(plot_simple_zones(com$canal, obs, paste0(title, " (WITH Vaccine)"), year_target),
                 paste0("band_", tolower(level), "_", rot_reg, "_", rot_tgt, "_", year_target, "_WITHvac"), folder)
  
  series <- data.frame(level = level, region = region, target = target_name, year = year_target, week = 1:52,
                       obs_acum = obs$acum,
                       pred_sem = sem$canal$media, li_sem = sem$canal$lower, ls_sem = sem$canal$upper,
                       pred_com = com$canal$media, li_com = com$canal$lower, ls_com = com$canal$upper, stringsAsFactors = FALSE)
  
  tryCatch(plot_predicted_vs_observed(series, title, year_target, level, rot_reg, rot_tgt),
           error = function(e) message("    [Plot Error] ", e$message))
  
  safenum <- function(x) if (isTRUE(length(x) > 0)) as.numeric(x[1]) else NA_real_
  
  metrics <- data.frame(
    level = level, region = region, target = target_name, year = year_target,
    obs_final = y_final,
    pred_final_sem = safenum(sem$canal$media[52]),
    pred_final_com = safenum(com$canal$media[52]),
    erro_final_sem = abs(safenum(sem$canal$media[52]) - y_final),
    erro_final_com = abs(safenum(com$canal$media[52]) - y_final),
    largura_banda_sem = safenum(sem$canal$upper[52]) - safenum(sem$canal$lower[52]),
    largura_banda_com = safenum(com$canal$upper[52]) - safenum(com$canal$lower[52]),
    coef_vac_primaria = safenum(com$coef_primaria),
    coef_vac_reforco  = safenum(com$coef_reforco),
    DIC_com  = safenum(com$dic),
    WAIC_com = safenum(com$waic),
    stringsAsFactors = FALSE
  )
  
  list(metrics = metrics, series = series, acum_sem = sem$acum_mat, acum_com = com$acum_mat)
}

# =============================================================================
#  7. AGGREGATE MENINGITIS CALCULATIONS (COMPLETE SUM + SWITCH LOGIC)
# =============================================================================

#' @title Compute Aggregate Impacts for Targets
#' @description Sums the MCMC matrices for specific components (with and without vaccine)
#'   to derive a precise total predicted effect. Integrates the 'menC_only' config switch.
compute_aggregate_target <- function(level, region, year_target, target_name, saved_models, vac_long) {
  ch <- function(al) paste(level, region, al, year_target)
  
  req_comps <- character(0)
  if (target_name %in% c("Total Meningitis", "Bacterial Meningitis")) {
    req_comps <- c("Haemophilus", "Meningococcal C", "Pneumococcal")
  } else if (target_name == "Meningococcal") {
    req_comps <- c("Meningococcal C")
  } else {
    return(NULL)
  }
  
  # === NOVO SWITCH DE EFEITO ISOLADO MENINGO C ===
  if (isTRUE(config$vaccine$effect_for_totals == "menC_only")) {
    req_comps <- c("Meningococcal C")
  }
  
  comps_list <- lapply(req_comps, function(c_name) saved_models[[ch(c_name)]])
  if (any(vapply(comps_list, is.null, logical(1)))) return(NULL)
  
  n_samples_min <- min(vapply(comps_list, function(x) ncol(x$acum_sem), integer(1)))
  clip_mat <- function(m) m[, seq_len(n_samples_min), drop = FALSE]
  
  sum_sem <- Reduce(`+`, lapply(comps_list, function(x) clip_mat(x$acum_sem)))
  sum_com <- Reduce(`+`, lapply(comps_list, function(x) clip_mat(x$acum_com)))
  
  df_all <- aggregate_weekly(base_menin, level, if (level == "State") NULL else region, target_name)
  train_years <- setdiff(config$train_years, year_target)
  df_train <- df_all %>% filter(year %in% train_years) %>% complete_weeks(train_years)
  
  if (sum(df_train$cases, na.rm = TRUE) < config$min_train_cases) return(NULL)
  
  tot_sem <- run_inla_acum_mat(df_train, config$n_samples, config$min_train_cases, config$inla_models[[level]],
                                seed = instance_seed(level, region, target_name, year_target, "sem"))
  if (is.null(tot_sem)) return(NULL)
  
  tot_mat <- clip_mat(tot_sem$acum_mat)
  
  # Arithmetic: Subtract Baseline sum, Add Vaccine sum (incorporates only active components)
  remainder <- tot_mat - sum_sem
  total_com <- remainder + sum_com
  total_sem <- tot_mat
  
  canal_sem <- matrix_zones(total_sem)
  canal_com <- matrix_zones(total_com)
  
  obs <- df_all %>% filter(year == year_target) %>% select(week, cases)
  obs <- data.frame(week = 1:52) %>% left_join(obs, by = "week") %>%
    mutate(cases = tidyr::replace_na(cases, 0L)) %>% arrange(week) %>% mutate(acum = cumsum(cases))
  y_final <- max(obs$acum, na.rm = TRUE)
  
  series <- data.frame(level = level, region = region, target = target_name, year = year_target, week = 1:52,
                       obs_acum = obs$acum,
                       pred_sem = canal_sem$media, li_sem = canal_sem$lower, ls_sem = canal_sem$upper,
                       pred_com = canal_com$media, li_com = canal_com$lower, ls_com = canal_com$upper, stringsAsFactors = FALSE)
  
  title <- paste0(target_name, " — Band ", year_target)
  rr <- if (level == "State") "ESP" else clean_name(region)
  
  tryCatch(plot_predicted_vs_observed(series, title, year_target, level, rr, clean_name(target_name)),
           error = function(e) message("    [Total Plot Error] ", e$message))
  
  metrics <- data.frame(level = level, region = region, target = target_name, year = year_target,
                        obs_final = y_final, pred_final_sem = canal_sem$media[52], pred_final_com = canal_com$media[52],
                        erro_final_sem = abs(canal_sem$media[52] - y_final), erro_final_com = abs(canal_com$media[52] - y_final),
                        largura_banda_sem = canal_sem$upper[52] - canal_sem$lower[52],
                        largura_banda_com = canal_com$upper[52] - canal_com$lower[52],
                        coef_vac_primaria = NA_real_, coef_vac_reforco = NA_real_,
                        DIC_com = NA_real_, WAIC_com = NA_real_, stringsAsFactors = FALSE)
  
  list(metrics = metrics, series = series)
}

# =============================================================================
#  8. TEMPORAL METRICS EXPORTER (4-week block TXT generator)
# =============================================================================

#' @title Export Temporal Block Metrics
#' @description Saves cumulative block (4-week intervals) metrics to a TXT file 
#'   specifically for the State Level, mapped out for every analyzed Target.
export_temporal_blocks <- function(series_all, output_folder) {
  if (nrow(series_all) == 0) return(invisible(NULL))
  
  # Filter only for the State (ESP)
  df <- series_all[series_all$level == "State", ]
  if (nrow(df) == 0) return(invisible(NULL))
  
  # Group weeks into Blocks of 4 (e.g. Block 1: W1-4, Block 2: W5-8...)
  df$block <- ceiling(df$week / 4)
  
  blocks_df <- df %>%
    dplyr::group_by(target, year, block) %>%
    dplyr::summarise(
      week_start   = min(week),
      week_end     = max(week),
      obs_acum_end = max(obs_acum),
      pred_sem_end = max(pred_sem),
      pred_com_end = max(pred_com),
      .groups      = "drop"
    ) %>%
    dplyr::group_by(target, year) %>%
    dplyr::arrange(block) %>%
    dplyr::mutate(
      obs_cases      = obs_acum_end - dplyr::lag(obs_acum_end, default = 0),
      pred_sem_cases = pred_sem_end - dplyr::lag(pred_sem_end, default = 0),
      pred_com_cases = pred_com_end - dplyr::lag(pred_com_end, default = 0),
      vaccine_effect = pred_sem_cases - pred_com_cases
    ) %>%
    dplyr::ungroup()
  
  file_path <- file.path(output_folder, "temporal_block_metrics_ESP.txt")
  sink(file_path)
  
  cat("===============================================================================\n")
  cat("           TEMPORAL BLOCK METRICS (STATE LEVEL - ESP)\n")
  cat("===============================================================================\n\n")
  
  targets <- unique(blocks_df$target)
  for (tgt in targets) {
    cat(sprintf("TARGET: %s\n", tgt))
    cat(strrep("-", 79), "\n")
    tgt_data <- blocks_df[blocks_df$target == tgt, ]
    
    years <- unique(tgt_data$year)
    for (yr in years) {
      cat(sprintf("  YEAR: %s\n", yr))
      cat(sprintf("  %-6s %-12s | %-12s | %-12s | %-12s | %-12s\n",
                  "Block", "Weeks", "Observed", "No Vaccine", "With Vaccine", "Vac Effect"))
      cat(paste0("  ", strrep("-", 75), "\n"))
      
      yr_data <- tgt_data[tgt_data$year == yr, ]
      for (i in 1:nrow(yr_data)) {
        row <- yr_data[i,]
        cat(sprintf("  %-6d %02d to %02d    | %-12.1f | %-12.1f | %-12.1f | %-12.1f\n",
                    row$block, row$week_start, row$week_end,
                    row$obs_cases, row$pred_sem_cases, row$pred_com_cases, row$vaccine_effect))
      }
      cat("\n")
    }
    cat("\n")
  }
  sink()
  message("\n  [Metrics] Temporal block metrics saved flawlessly to:\n    -> ", file_path)
}

# =============================================================================
#  9. EXECUTION ROUTINE
# =============================================================================

execute_vaccine_comparison <- function() {
  cfg <- config$vaccine
  message("\n===== BANDS WITH/WITHOUT VACCINE — COMPLETE SUM =====\n")
  vac_long <- load_vaccines()
  if (!length(vac_long)) { message("  [Vaccine] No CSV loaded; aborting."); return(invisible(NULL)) }
  
  targets_vac <- intersect(names(cfg$cov_by_target), names(config$targets))
  levels <- cfg$levels; if (isTRUE(cfg$skip_municipality)) levels <- setdiff(levels, "Municipality")
  
  rows <- list()
  series_all <- list()
  saved_models <- list()
  
  # Iterate over Standard Vaccine Targets
  for (level in levels) {
    regions <- switch(level, "State" = "SP", "DRS" = list_regions("DRS"),
                      "CIR" = list_regions("CIR"), "Municipality" = list_regions("Municipality"))
    for (region in regions) {
      for (target_name in targets_vac) {
        for (year_target in cfg$band_years) {
          res <- tryCatch(
            compare_vaccine_instance(level, region, target_name, year_target, vac_long),
            error = function(e) { message("    [ERROR vac cmp] ", level, "/", region, "/", target_name, ": ", e$message); NULL }
          )
          if (!is.null(res)) {
            ch <- paste(level, region, target_name, year_target)
            rows[[ch]] <- res$metrics
            series_all[[ch]] <- res$series
            saved_models[[ch]] <- list(acum_sem = res$acum_sem, acum_com = res$acum_com)
          }
        }
      }
    }
  }
  
  # AGGREGATE TOTALS
  aggregate_targets <- c("Total Meningitis", "Bacterial Meningitis", "Meningococcal")
  active_aggregates <- intersect(aggregate_targets, names(config$targets))
  
  if (length(active_aggregates) > 0) {
    for (level in levels) {
      regions <- switch(level, "State" = "SP", "DRS" = list_regions("DRS"),
                        "CIR" = list_regions("CIR"), "Municipality" = list_regions("Municipality"))
      for (region in regions) {
        for (year_target in cfg$band_years) {
          for (tgt in active_aggregates) {
            tot <- tryCatch(
              compute_aggregate_target(level, region, year_target, tgt, saved_models, vac_long),
              error = function(e) { message("    [Total Vac ERROR] ", e$message); NULL }
            )
            if (!is.null(tot)) {
              ch <- paste(level, region, tgt, year_target)
              rows[[ch]] <- tot$metrics
              series_all[[ch]] <- tot$series
            }
          }
        }
      }
    }
  }
  
  # Saving all logic
  if (length(rows)) {
    folder <- config$output_folder
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    
    utils::write.csv(do.call(rbind, rows), file.path(folder, "comparison_with_without_vaccine.csv"),
                     row.names = FALSE, fileEncoding = "UTF-8")
    
    if (length(series_all)) {
      df_series <- do.call(rbind, series_all)
      utils::write.csv(df_series, file.path(folder, "series_predicted_observed.csv"),
                       row.names = FALSE, fileEncoding = "UTF-8")
      
      # Generate the 4-week temporal block metrics TXT.
      export_temporal_blocks(df_series, folder)
    }
    
    message("\n=== Vaccine comparison generated and saved to:\n    ", normalizePath(folder, mustWork = FALSE), " ===")
  } else {
    message("\n[Vaccine] No instances produced comparisons.\n")
  }
  
  # Write out the diagnostic log accumulated during this routine.
  tryCatch(vaclog_flush(config$output_folder, config$vaccine$diagnostic_log_file %||% "vaccine_diagnostic_log.txt"),
           error = function(e) message("    [Vaccine Diagnostics] failed to save log: ", e$message))
}

# NOTE: the Besag (joint spatial) routine is omitted from this execution block to
# keep it readable. It follows the same logic and can be appended at the end.

execute_vaccine_comparison()


# =============================================================================
#  10. VALIDATION METRICS MODULE  (self-contained — paste at the END of script)
# -----------------------------------------------------------------------------
#  Computes, for every instance (level / region / target / year) and for BOTH
#  models (NO vaccine vs WITH vaccine):
#
#     MAE | RMSE | Coverage | WIS | Relative WIS | Final Error | MAPE | DIC | WAIC
#
#  Math is identical to metricas_acum() from the historical-band script.
#  Results are written to a TXT (plus a CSV for further analysis).
#
#  NOTE: the module re-fits the INLA models, because WIS and Coverage require the
#  full posterior sample matrix, which is not kept after execute_vaccine_comparison().
#  Use metrics_cfg below to restrict levels/targets/years and keep runtime sane.
# =============================================================================

metrics_cfg <- list(
  coverage_level = 0.90,          # nominal central interval used for Coverage
  levels         = c("State"),    # c("State", "DRS") — DRS multiplies runtime a lot
  targets        = "all",         # "all" or e.g. c("Meningococcal C", "Pneumococcal")
  years          = NULL,          # NULL = uses config$vaccine$band_years
  include_aggregates = TRUE,      # Total / Bacterial / Meningococcal (derived matrices)
  n_samples      = NULL,          # NULL = uses config$n_samples
  file_txt       = "validation_metrics_vaccine.txt",
  file_csv       = "validation_metrics_vaccine.csv"
)

# ---------------------------------------------------------------------------
#  10.1 Core metric engine
# ---------------------------------------------------------------------------
metrics_acum <- function(acum_mat, obs_cases, dic = NA_real_, waic = NA_real_,
                         cov_level = 0.90) {
  if (is.null(acum_mat)) return(NULL)
  if (is.null(dim(acum_mat))) acum_mat <- matrix(acum_mat, nrow = length(obs_cases))
  
  n <- min(nrow(acum_mat), length(obs_cases))
  if (n < 1) return(NULL)
  acum_mat <- acum_mat[seq_len(n), , drop = FALSE]
  
  y_acum    <- cumsum(obs_cases[seq_len(n)])
  pred_mean <- rowMeans(acum_mat)
  
  mae  <- mean(abs(y_acum - pred_mean))
  rmse <- sqrt(mean((y_acum - pred_mean)^2))
  
  a     <- (1 - cov_level) / 2
  qlow  <- apply(acum_mat, 1, quantile, probs = a,     names = FALSE, na.rm = TRUE)
  qhigh <- apply(acum_mat, 1, quantile, probs = 1 - a, names = FALSE, na.rm = TRUE)
  coverage <- mean(y_acum >= qlow & y_acum <= qhigh) * 100
  
  alphas <- 1 - c(0.50, 0.80, 0.90, 0.95); K <- length(alphas)
  wis_i <- vapply(seq_len(n), function(i) {
    am   <- acum_mat[i, ]
    soma <- 0.5 * abs(y_acum[i] - stats::median(am, na.rm = TRUE))
    for (al in alphas) {
      l <- quantile(am, al / 2,     names = FALSE, na.rm = TRUE)
      u <- quantile(am, 1 - al / 2, names = FALSE, na.rm = TRUE)
      is_a <- (u - l) +
        (2 / al) * (l - y_acum[i]) * (y_acum[i] < l) +
        (2 / al) * (y_acum[i] - u) * (y_acum[i] > u)
      soma <- soma + (al / 2) * is_a
    }
    soma / (K + 0.5)
  }, numeric(1))
  
  wis     <- mean(wis_i, na.rm = TRUE)
  escala  <- mean(y_acum, na.rm = TRUE)
  wis_rel <- if (is.finite(escala) && escala > 0) wis / escala else NA_real_
  
  erro_final <- abs(pred_mean[n] - y_acum[n])
  mape_final <- if (y_acum[n] > 0) erro_final / y_acum[n] * 100 else NA_real_
  
  data.frame(MAE = mae, RMSE = rmse, Coverage = coverage,
             WIS = wis, WIS_rel = wis_rel,
             FinalError = erro_final, MAPE_final = mape_final,
             DIC = as.numeric(dic %||% NA_real_), WAIC = as.numeric(waic %||% NA_real_),
             n_weeks = n, obs_final = y_acum[n], pred_final = pred_mean[n],
             stringsAsFactors = FALSE)
}

metrics_row <- function(level, region, target, year, model, m) {
  if (is.null(m)) return(NULL)
  cbind(data.frame(level = level, region = region, target = target,
                   year = year, model = model, stringsAsFactors = FALSE), m)
}

# ---------------------------------------------------------------------------
#  10.2 Helpers: observed weekly cases and training frame for one instance
# ---------------------------------------------------------------------------
metrics_obs_cases <- function(level, region, target_name, year_target) {
  df_all <- aggregate_weekly(base_menin, level,
                             if (level == "State") NULL else region, target_name)
  o <- df_all %>% filter(year == year_target) %>% select(week, cases)
  o <- data.frame(week = 1:52) %>% left_join(o, by = "week") %>%
    mutate(cases = tidyr::replace_na(cases, 0L)) %>% arrange(week)
  o$cases
}

metrics_train_df <- function(level, region, target_name, year_target) {
  df_all <- aggregate_weekly(base_menin, level,
                             if (level == "State") NULL else region, target_name)
  ty <- setdiff(config$train_years, year_target)
  df_all %>% filter(year %in% ty) %>% complete_weeks(ty)
}

# ---------------------------------------------------------------------------
#  10.3 Main routine
# ---------------------------------------------------------------------------
run_metrics_module <- function() {
  message("\n===== VALIDATION METRICS MODULE =====\n")
  
  cfg   <- config$vaccine
  n_smp <- metrics_cfg$n_samples %||% config$n_samples
  cov_l <- metrics_cfg$coverage_level %||% 0.90
  years <- metrics_cfg$years %||% cfg$band_years
  levels_run <- intersect(metrics_cfg$levels, c("State", "DRS", "CIR", "Municipality"))
  
  vac_long <- load_vaccines()
  if (!length(vac_long)) { message("  [Metrics] No vaccine CSV loaded; aborting."); return(invisible(NULL)) }
  
  targets_vac <- intersect(names(cfg$cov_by_target), names(config$targets))
  if (!identical(metrics_cfg$targets, "all"))
    targets_vac <- intersect(targets_vac, metrics_cfg$targets)
  
  rows <- list(); store <- list()
  
  # ---- (A) direct vaccine targets ------------------------------------------
  for (level in levels_run) {
    regions <- switch(level, "State" = "SP", "DRS" = list_regions("DRS"),
                      "CIR" = list_regions("CIR"), "Municipality" = list_regions("Municipality"))
    for (region in regions) for (target_name in targets_vac) for (yt in years) {
      
      key <- paste(level, region, target_name, yt)
      message("  [Metrics] ", level, ":", region, " | ", target_name, " | ", yt)
      
      res <- tryCatch({
        df_train <- metrics_train_df(level, region, target_name, yt)
        if (sum(df_train$cases, na.rm = TRUE) < config$min_train_cases) stop("few cases")
        
        cov_tab <- covariate_target(vac_long, target_name, level, region)
        if (is.null(cov_tab)) stop("no coverage data")
        # Same lag as in compare_vaccine_instance(), so the metrics module
        # evaluates EXACTLY the design that produced the published figures.
        cov_tab <- shift_coverage_years(cov_tab, config$vaccine$coverage_lag_years %||% 0L)
        cov_yt <- get_year_covariate(cov_tab, yt)
        if (all(is.na(unlist(cov_yt)))) stop("no coverage in target year")
        
        sem <- run_inla_acum_mat(df_train, n_smp, config$min_train_cases,
                                 config$inla_models[[level]],
                                 seed = instance_seed(level, region, target_name, yt, "sem"))
        com <- run_inla_acum_mat_vac(df_train, cov_tab[, c("ano", "cov_primaria", "cov_reforco")],
                                     cov_yt, n_smp, config$min_train_cases,
                                     config$inla_models[[level]], cfg$standardize,
                                     seed = instance_seed(level, region, target_name, yt, "com"))
        if (is.null(sem) || is.null(com)) stop("a model failed")
        list(sem = sem, com = com)
      }, error = function(e) { message("    skipped: ", conditionMessage(e)); NULL })
      
      if (is.null(res)) next
      
      obs <- metrics_obs_cases(level, region, target_name, yt)
      store[[key]] <- list(acum_sem = res$sem$acum_mat, acum_com = res$com$acum_mat)
      
      rows[[paste0(key, "|sem")]] <- metrics_row(level, region, target_name, yt, "NO vaccine",
                                                 metrics_acum(res$sem$acum_mat, obs, res$sem$dic, res$sem$waic, cov_l))
      rows[[paste0(key, "|com")]] <- metrics_row(level, region, target_name, yt, "WITH vaccine",
                                                 metrics_acum(res$com$acum_mat, obs, res$com$dic, res$com$waic, cov_l))
    }
  }
  
  # ---- (B) aggregate targets (same arithmetic as compute_aggregate_target) --
  if (isTRUE(metrics_cfg$include_aggregates)) {
    aggs <- intersect(c("Total Meningitis", "Bacterial Meningitis", "Meningococcal"),
                      names(config$targets))
    for (level in levels_run) {
      regions <- switch(level, "State" = "SP", "DRS" = list_regions("DRS"),
                        "CIR" = list_regions("CIR"), "Municipality" = list_regions("Municipality"))
      for (region in regions) for (tgt in aggs) for (yt in years) {
        
        comps <- if (tgt == "Meningococcal") "Meningococcal C"
        else c("Haemophilus", "Meningococcal C", "Pneumococcal")
        if (identical(cfg$effect_for_totals, "menC_only")) comps <- "Meningococcal C"
        
        got <- lapply(comps, function(cn) store[[paste(level, region, cn, yt)]])
        if (any(vapply(got, is.null, logical(1)))) next
        
        message("  [Metrics] ", level, ":", region, " | ", tgt, " (aggregate) | ", yt)
        
        res <- tryCatch({
          df_train <- metrics_train_df(level, region, tgt, yt)
          if (sum(df_train$cases, na.rm = TRUE) < config$min_train_cases) stop("few cases")
          tot <- run_inla_acum_mat(df_train, n_smp, config$min_train_cases,
                                   config$inla_models[[level]],
                                   seed = instance_seed(level, region, tgt, yt, "sem"))
          if (is.null(tot)) stop("total model failed")
          
          k    <- min(c(ncol(tot$acum_mat), vapply(got, function(x) ncol(x$acum_sem), integer(1))))
          clip <- function(m) m[, seq_len(k), drop = FALSE]
          
          sum_sem <- Reduce(`+`, lapply(got, function(x) clip(x$acum_sem)))
          sum_com <- Reduce(`+`, lapply(got, function(x) clip(x$acum_com)))
          tot_mat <- clip(tot$acum_mat)
          
          list(sem = tot_mat, com = (tot_mat - sum_sem) + sum_com,
               dic = tot$dic, waic = tot$waic)
        }, error = function(e) { message("    skipped: ", conditionMessage(e)); NULL })
        
        if (is.null(res)) next
        obs <- metrics_obs_cases(level, region, tgt, yt)
        key <- paste(level, region, tgt, yt)
        
        rows[[paste0(key, "|sem")]] <- metrics_row(level, region, tgt, yt, "NO vaccine",
                                                   metrics_acum(res$sem, obs, res$dic, res$waic, cov_l))
        rows[[paste0(key, "|com")]] <- metrics_row(level, region, tgt, yt, "WITH vaccine",
                                                   metrics_acum(res$com, obs, NA_real_, NA_real_, cov_l))
      }
    }
  }
  
  if (!length(rows)) { message("\n  [Metrics] No instance produced metrics.\n"); return(invisible(NULL)) }
  export_validation_metrics(do.call(rbind, rows), config$output_folder)
}

# ---------------------------------------------------------------------------
#  10.4 TXT / CSV exporter
# ---------------------------------------------------------------------------
export_validation_metrics <- function(tab, output_folder) {
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  tab <- tab[order(tab$level, tab$region, tab$target, tab$year, tab$model), ]
  
  utils::write.csv(tab, file.path(output_folder, metrics_cfg$file_csv),
                   row.names = FALSE, fileEncoding = "UTF-8")
  
  path <- file.path(output_folder, metrics_cfg$file_txt)
  con  <- file(path, open = "w", encoding = "UTF-8"); on.exit(close(con), add = TRUE)
  cl   <- round(100 * (metrics_cfg$coverage_level %||% 0.90))
  
  writeLines(c(
    strrep("=", 112),
    " MENINGITIS SP — VACCINE BANDS | PREDICTIVE VALIDATION METRICS (CUMULATIVE)",
    paste0(" Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
    paste0(" Coverage evaluated at the ", cl, "% central posterior interval."),
    " MAE / RMSE : error of the cumulative curve across the 52 epidemiological weeks.",
    " WIS        : Weighted Interval Score (50/80/90/95%); WIS_rel = WIS / mean(observed).",
    " FinalError : |predicted - observed| cumulative total at EW52; MAPE_final in %.",
    " DIC / WAIC : INLA information criteria (NA for derived aggregate WITH-vaccine curves).",
    strrep("=", 112), ""), con)
  
  f <- function(x, d = 3) if (!is.finite(x)) "      NA" else formatC(x, format = "f", digits = d)
  hdr <- sprintf("    %-5s %-14s %8s %8s %9s %8s %9s %10s %9s %9s %9s",
                 "Year", "Model", "MAE", "RMSE", "Coverage", "WIS", "WIS_rel",
                 "FinalErr", "MAPE(%)", "DIC", "WAIC")
  
  for (lv in unique(tab$level)) {
    writeLines(c(paste0("LEVEL: ", lv), strrep("-", 112)), con)
    d_lv <- tab[tab$level == lv, ]
    for (rg in unique(d_lv$region)) {
      writeLines(paste0("  REGION: ", rg), con)
      d_rg <- d_lv[d_lv$region == rg, ]
      for (tg in unique(d_rg$target)) {
        writeLines(c(paste0("    TARGET: ", tg), hdr, paste0("    ", strrep("-", 108))), con)
        d_tg <- d_rg[d_rg$target == tg, ]
        for (i in seq_len(nrow(d_tg))) {
          r <- d_tg[i, ]
          writeLines(sprintf("    %-5s %-14s %8s %8s %9s %8s %9s %10s %9s %9s %9s",
                             r$year, r$model, f(r$MAE, 2), f(r$RMSE, 2), f(r$Coverage, 1),
                             f(r$WIS, 3), f(r$WIS_rel, 4), f(r$FinalError, 2),
                             f(r$MAPE_final, 2), f(r$DIC, 1), f(r$WAIC, 1)), con)
        }
        writeLines("", con)
      }
    }
    writeLines("", con)
  }
  
  writeLines(c(strrep("=", 112),
               " OVERALL SUMMARY (mean across all instances, by model)",
               strrep("=", 112),
               sprintf(" %-14s %8s %8s %9s %8s %9s %10s %9s", "Model", "MAE", "RMSE",
                       "Coverage", "WIS", "WIS_rel", "FinalErr", "MAPE(%)")), con)
  for (md in unique(tab$model)) {
    d <- tab[tab$model == md, ]
    mn <- function(cc) mean(d[[cc]], na.rm = TRUE)
    writeLines(sprintf(" %-14s %8s %8s %9s %8s %9s %10s %9s", md,
                       f(mn("MAE"), 2), f(mn("RMSE"), 2), f(mn("Coverage"), 1),
                       f(mn("WIS"), 3), f(mn("WIS_rel"), 4),
                       f(mn("FinalError"), 2), f(mn("MAPE_final"), 2)), con)
  }
  writeLines("", con)
  
  message("\n  [Metrics] Saved to:\n    -> ", path,
         "\n    -> ", file.path(output_folder, metrics_cfg$file_csv))
  invisible(path)
}

# ---------------------------------------------------------------------------
#  10.5 RUN
# ---------------------------------------------------------------------------
run_metrics_module()

# =============================================================================
#  11. COMBINED FIGURE (VACCINE) — side-by-side panels, single shared legend
# =============================================================================
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))

# Etiologies shown in the combined figure (only those with coverage data).
combo_targets_vac <- c("Meningococcal C", "Pneumococcal")

gerar_figura_combinada_vacina <- function(year_target, cenario = c("sem", "com")) {
  cenario  <- match.arg(cenario)
  vac_long <- load_vaccines()
  letras   <- letters[seq_along(combo_targets_vac)]
  
  paineis <- Map(function(al, lt) {
    res <- tryCatch(compare_vaccine_instance("State", "SP", al, year_target, vac_long), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    s <- res$series
    canal <- if (cenario == "sem")
      data.frame(week = s$week, lower = s$li_sem, media = s$pred_sem, upper = s$ls_sem)
    else
      data.frame(week = s$week, lower = s$li_com, media = s$pred_com, upper = s$ls_com)
    obs <- data.frame(week = s$week, acum = s$obs_acum)
    plot_simple_zones(canal, obs, paste0(lt, ") ", al), year_target) + theme(plot.subtitle = element_blank())
  }, combo_targets_vac, letras)
  
  paineis <- Filter(Negate(is.null), paineis)
  if (!length(paineis)) { message("  [Combined vaccine figure] no panel was generated."); return(invisible(NULL)) }
  
  combinada <- patchwork::wrap_plots(paineis, ncol = 2) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  pasta <- file.path(config$output_folder, as.character(year_target), "Figura_Combinada")
  if (!dir.exists(pasta)) dir.create(pasta, recursive = TRUE)
  rot <- if (cenario == "sem") "NOvac" else "WITHvac"
  arq <- file.path(pasta, paste0("figura_combinada_", rot, "_", year_target, ".png"))
  ggsave(arq, combinada, width = 2000/150, height = 1200/150, units = "in", dpi = 150, bg = "white")
  message("  -> Combined figure (", rot, ") saved to: ", arq)
}

for (yt in config$vaccine$band_years) {
  gerar_figura_combinada_vacina(yt, "sem")
  gerar_figura_combinada_vacina(yt, "com")
}
