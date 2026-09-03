# =============================================================================
#  HISTORICAL BAND (ENDEMIC CHANNEL) FOR MENINGITIS — State of São Paulo, Brazil
#  Author: adapted from a SRAG script (Gabriel Vinicius Vian) for meningitis
# =============================================================================
#
#  WHAT THIS SCRIPT DOES
#  ---------------------------------------------------------------------------
#  It builds an "endemic channel" (historical band) for meningitis surveillance:
#  a four-zone reference band (Below expected / Expected / Alert / Epidemic)
#  showing how many CUMULATIVE cases would be expected across the 52
#  epidemiological weeks of a year, given the historical pattern. The observed
#  trajectory of a chosen target year is overlaid on that band, so an
#  epidemiologist can see at a glance whether the year is running low, normal,
#  or above historical expectation.
#
#  TWO ENGINES ARE AVAILABLE (config$metodo_banda)
#    "inla"             Bayesian model (weekly Poisson -> cumulative) fitted with
#                       R-INLA; the four zones are posterior percentiles. This is
#                       the method used in our analysis.
#    "canal_geometrico" Classic geometric endemic channel (Bortman 1999),
#                       reproducing the endemic.channels package exactly. Kept as
#                       a benchmark for side-by-side comparison.
#
#  THE MODEL, IN ONE PARAGRAPH
#  ---------------------------------------------------------------------------
#  Weekly case counts are modelled as Poisson. The linear predictor has an
#  intercept, a seasonal term across the 52 weeks (first-order random walk by
#  default, or Fourier harmonics), and an i.i.d. year effect. To predict a new
#  year, a block of rows carrying a NEW year index is appended, so the target
#  year's effect is drawn from its prior instead of being fitted -- this is what
#  makes the band a genuine PREDICTIVE interval rather than a description of the
#  training years. We then draw N samples from the posterior, simulate weekly
#  counts from each sample, and take the CUMULATIVE SUM within each sample. The
#  four zones are quantiles of the resulting 52 x N matrix.
#
#  OPTIONAL COMPONENTS
#  ---------------------------------------------------------------------------
#    * Besag spatial model  : fits all regions of a level jointly, with a spatial
#                             random effect from a shapefile adjacency graph.
#    * State reconstruction : fits each sub-region and sums the cumulative
#                             trajectories SAMPLE BY SAMPLE. Summing by sample is
#                             what propagates uncertainty correctly; summing
#                             quantiles would be statistically invalid.
#    * Contingency analysis : overlays how an already-closed year LOOKED at
#                             several moments during that year, using partial
#                             database snapshots -- i.e. visualises reporting delay.
#    * Validation           : leave-one-year-out and forward-chaining (block)
#                             validation with MAE, RMSE, interval coverage, WIS,
#                             final cumulative error, MAPE, DIC and WAIC.
#
#  ABOUT THE INPUT DATA   <<< IMPORTANT IF YOU ARE REUSING THIS CODE >>>
#  ---------------------------------------------------------------------------
#  This script reads an ALREADY-CURATED line-list (.RDS): a dataset in which raw
#  SINAN records (Brazil's notifiable disease information system) have already
#  been cleaned and decoded into readable columns.
#
#  We produced that curated file OURSELVES, starting from RAW SINAN data. Adapting
#  this pipeline to raw SINAN is therefore straightforward: all you need is a
#  preprocessing step that decodes the raw SINAN fields into the column names
#  listed below. Nothing downstream has to change, because every routine here
#  consumes those column names -- not the raw SINAN numeric codes.
#
#  EXPECTED COLUMNS (one row per notified case):
#    anoepi, sepi     epidemiological year and week (integers). If missing, they
#                     are derived from data_sintoma (see config$fonte_epi).
#    data_sintoma     symptom onset date (fallback source for year/week).
#    classi_fin       final classification; "Confirmado" marks a confirmed case.
#    grupo_etiologico broad etiological group (e.g. "Asseptica", "Bacteriana").
#    etiologia        agent (e.g. "Meningococo", "Pneumococo", "Hemofilo").
#    sorogrupo        meningococcal serogroup (e.g. "Sorogrupo B", "Sorogrupo C").
#    MUNI_RESID       municipality of residence (name).
#    DRS_RESID        health region of residence (name); optional.
#
#  HOW TO RUN
#  ---------------------------------------------------------------------------
#  1. Replace every "INSERT HERE YOUR DIRECTORY" placeholder below with a real
#     path on your machine.
#  2. Set config$anos_treino (training years) and config$ano_atual (the year whose
#     observed trajectory is drawn on top of the band).
#  3. Source the whole file. Missing CRAN packages install automatically; R-INLA
#     is fetched from its own repository on first run.
#
#  REPRODUCIBILITY: every posterior sampling call is seeded deterministically
#  from the instance identity (see instance_seed below), so re-running the script
#  reproduces all figures and metrics exactly.
#
#  KNOWN LIMITATION: years with 53 epidemiological weeks have their week-53 cases
#  dropped, so that every year contributes exactly 52 points and years can be
#  stacked in a single matrix. This follows the endemic.channels convention, but
#  it does mean annual totals are marginally underestimated in those years.
# =============================================================================

# =============================================================================
#  1) CONFIGURATIONS  —  EDIT HERE ✍️
# =============================================================================
config <- list(
  
  # Meningitis line-list database (SINAN / CVE-SP):
  caminho_dados = "INSERT HERE YOUR DIRECTORY",
  pasta_saida   = "INSERT HERE YOUR DIRECTORY",
  
  # Training: 2013 to 2024, EXCLUDING 2020-2022 (your epidemiological decision).
  anos_treino = c(2013, 2014, 2015, 2016, 2017, 2018, 2023, 2024, 2025),
  ano_atual   = 2019,
  
  # Filter for CONFIRMED cases (classi_fin == "Confirmado"). Now it's PER TARGET:
  #   - each target can have 'confirmado = TRUE/FALSE' (see targets list);
  #   - if the target DOES NOT define it, this global DEFAULT is used.
  # Default FALSE = mirrors the channel doc (B/C without filter). "Total Meningitis"
  # is marked as confirmado = TRUE (as requested). The diagnostic shows
  # ALL vs CONFIRMED per year for verification.
  apenas_confirmados = TRUE,
  
  # =========================================================================
  # 1.0 BAND METHOD 🚫
  # =========================================================================
  #   "inla"             = Bayesian band (Option A: Weekly Poisson -> cumulative).
  #                        This is OUR method (kept). The ZONES/percentiles are the
  #                        same as the endemic channel (see inla_percentis below).
  #   "canal_geometrico" = Pure geometric endemic channel (endemic.channels).
  metodo_banda = "inla", ## <<< LEAVE ONLY INLA HERE!
  
  # Percentiles extracted from the INLA POSTERIOR to form the 4 channel zones
  # (below/expected/alert/epidemic), MIRRORING the endemic channel (95% interval):
  #   lower = 0.025, media = 0.50, upper = 0.975.
  # (change here for another level; e.g., c(lower=0.05, media=0.50, upper=0.95) for 90%.)
  inla_percentis = c(lower = 0.025, media = 0.50, upper = 0.975),
  
  # Pure geometric channel parameters (only if metodo_banda = "canal_geometrico").
  canal = list(
    nivel_confianca = 0.95,
    tipo_intervalo  = "ic_media",  # "ic_media" (± t*s/sqrt(n), package) | "preditivo" (± t*s)
    offset_log      = 1
  ),
  
  # =========================================================================
  # 1.0.1 EPIDEMIOLOGICAL YEAR/WEEK SOURCE 🚫
  # =========================================================================
  #   "auto"         = uses anoepi/sepi from the database if they exist (RECOMMENDED; matches
  #                    the endemic channels code), otherwise derives from data_sintoma.
  #   "data_sintoma" = forces derivation from data_sintoma (epiyear/epiweek).
  fonte_epi = "auto",
  
  # Print case count diagnostics (ALL vs CONFIRMED per year)?
  diagnostico_casos = TRUE,
  
  # =========================================================================
  # 1.1 EXECUTION MODE 🚫
  # =========================================================================
  # rodar_besag (only applies to metodo_banda = "inla"; the geometric channel is always
  # calculated region by region, without spatial pooling):
  #   TRUE  = Joint Spatial Fit (Besag).  FALSE = Independent Fit.
  rodar_besag = FALSE,
  
  # =============================== ==========================================
  # 1.2 GEOGRAPHIC SCOPE 🚫  
  # =========================================================================
  niveis = c("Estado"),   # "Estado" | "DRS" | "CIR" | "Municipio"
  
  drs_selecionadas = NULL,
  cir_selecionadas = NULL,
  mun_selecionados = NULL,
  
  # =========================================================================
  # 1.3 EPIDEMIOLOGICAL TARGETS (definition of "meningitis") 🚫
  # =========================================================================
  # Rules per target (AND between fields):
  #   list(tipo="todas") | list(grupo_etiologico=c(...)) | list(etiologia=c(...))
  #   | list(sorogrupo=c(...))  (combine fields)
  # Optional field 'confirmado = TRUE/FALSE' overrides the global default apenas_confirmados.
  alvos = list(
    "Total Meningitis"     = list(tipo = "todas", confirmado = TRUE),  # <<< as requested: only confirmed
    "Viral Meningitis"     = list(grupo_etiologico = c("Asseptica")),
    "Bacterial Meningitis" = list(grupo_etiologico = c("Bacteriana", "Nao determinada")),
    "Meningococcus"        = list(etiologia = "Meningococo"),
    "Meningococcus B"      = list(grupo_etiologico = "Bacteriana", etiologia = "Meningococo", sorogrupo = "Sorogrupo B"),
    "Meningococcus C"      = list(grupo_etiologico = "Bacteriana", etiologia = "Meningococo", sorogrupo = "Sorogrupo C"),
    "Pneumococcus"         = list(etiologia = "Pneumococo"),
    "Haemophilus"          = list(etiologia = "Hemofilo")
  ), 
  
  n_amostras = 1000,   # (only used in the "inla" method)
  
  # =========================================================================
  # 1.4 Minimum threshold of cases in training 🚫
  # =========================================================================
  min_casos_treino = 30,
  
  # =========================================================================
  # 1.5 RECONSTRUCTION OF THE STATE FROM A GRANULAR LEVEL (DRS/CIR/Municipality) 🚫
  # =========================================================================
  # Instead of fitting a single model for the State, it fits a model PER REGION
  # of the granular level and SUMS the cumulative trajectories SAMPLE BY SAMPLE from
  # the posterior (sum of cumsum = cumsum of sum; summing by sample preserves uncertainty,
  # summing quantiles would NOT be correct). Generates the same 4 zones.
  reconstrucao = list(
    ativar     = FALSE,
    nivel_grao = "DRS",          # "DRS" | "CIR" | "Municipio"
    # Which versions to generate: FALSE = independent fits per region;
    # TRUE = joint spatial fit (Besag). You can request both.
    besag      = c(FALSE, TRUE),
    alvos      = c("Total Meningitis"),   # or "Todos"
    # Region of the granular level with few cases: 0 = fits all regardless (recommended,
    # otherwise skipped ones enter as ZERO and the total gets underestimated).
    min_casos_regiao = 0,
    
    # BLOCK VALIDATION (forward-chaining) of the reconstructed State
    validacao_block = FALSE,
    val_n_amostras  = 500,
    val_min_treino  = 3,          # min number of years in training for 1st fold
    val_nivel_cobertura = 0.90
  ),
  
  # =========================================================================
  # 1.5.1 CONTINGENCY ANALYSIS (evolution of a consolidated year) ✍️
  # =========================================================================
  # Overlays the accumulated trajectories of an already CONSOLIDATED year (e.g., 2025)
  # onto the HISTORICAL BAND (INLA), seen at various moments of the year via
  # PARTIAL BASES by EW (13/26/39/50) + the consolidated one. Same logic as the
  # endemic channels code, but the BACKGROUND is our historical band (4 zones), NOT the channel.
  # The background band is trained WITHOUT the consolidated year.
  contingencia = list(
    ativar    = TRUE,                 # <<< SWITCH turns on/off
    ano       = 2025,                  # consolidated year to overlay
    # Partial bases (line-list, SAME format as the database). Label -> .RDS path.
    # The order defines the line order (increasing EW). (These partial bases are not available in this repository)
    pasta_parciais = "INSERT HERE YOUR DIRECTORY",
    parciais = list(
      "EW 13" = "base_parcial_SE13_2025.RDS",
      "EW 26" = "base_parcial_SE26_2025.RDS",
      "EW 39" = "base_parcial_SE39_2025.RDS",
      "EW 50" = "base_parcial_SE50_2025.RDS"
    ),
    # up to which EW each partial is considered "known" (NULL = infers from label)
    corte_por_parcial = NULL,
    incluir_consolidada = TRUE,        # also overlays the closed year (current base)
    alvos     = "Todos",               # "Todos" or vector of target names
    nivel     = "Estado",              # "Estado" (uses the whole database)
    # partial line colors (the consolidated one is always black)
    cores_parciais = c("EW 13" = "#9ECAE1", "EW 26" = "#4292C6",
                       "EW 39" = "#08519C", "EW 50" = "#082E54")
  ),
      
  carinhas       = c(typical = "🙂", high = "😟", very_high = "😢"),
      
  # "inla" method
  # 4 geometric channel zones (mirrors the package: Below/Expected/Alert/Epidemic)
  carinhas_canal = c(below = "🙂", expected = "😐", alert = "😟", epidemic = "😱"),
  
  # =========================================================================
  # 1.7 INLA MODELS BY LEVEL (only used if metodo_banda = "inla")   🚫
  # =========================================================================
  modelos_inla = list(
    Estado    = list(familia = "poisson", sazonal = "rw1", n_harmonicos = 2, formula = NULL),
    DRS       = list(familia = "poisson", sazonal = "rw1", n_harmonicos = 2, formula = NULL),
    CIR       = list(familia = "poisson", sazonal = "rw1", n_harmonicos = 1, formula = NULL),
    Municipio = list(familia = "poisson", sazonal = "rw1", n_harmonicos = 1, formula = NULL)
  ),
  
  # =========================================================================
  # 1.8 TOPOLOGY AND SPATIAL EFFECT (Besag) — only if rodar_besag = TRUE and method is "inla"  🚫
  # =========================================================================
  besag = list(
    familia = NULL, sazonal = NULL, n_harmonicos = NULL, formula = NULL,  # NULL = inherits from level
    # YOUR SHAPEFILES:
    usar_geobr = FALSE, ano_geobr = 2020,
    shapefile_municipios = "INSERT HERE YOUR DIRECTORY",
    coluna_id_municipios = NULL,
    shapefile_drs = "INSERT HERE YOUR DIRECTORY",
    coluna_id_drs = "X17DRS",
    shapefile_cir = "INSERT HERE YOUR DIRECTORY",
    coluna_id_cir = NULL
  ),
  grafos        = list(Estado = NULL, DRS = NULL, CIR = NULL, Municipio = NULL),
  regioes_grafo = list(Estado = NULL, DRS = NULL, CIR = NULL, Municipio = NULL),
  
  cir = list(
    shapefile = "INSERT HERE YOUR DIRECTORY",
    coluna_id = NULL, lookup = NULL
  ),
  
  # =========================================================================
  # 1.9 VALIDATION — uses the INLA model (Option A) as a predictive benchmark,
  #     INDEPENDENT of the method used to draw the band.✍️
  # =========================================================================
  validacao = list(
    ativar         = TRUE,
    incluir_estado = TRUE,
    drs            =  character(0),
    cir            = character(0),
    municipios     = character(0),
    alvos          =   c("Total Meningitis", "Meningococcus", "Viral Meningitis",  "Bacterial Meningitis", "Meningococcus", "Meningococcus B",  "Meningococcus C" , "Pneumococcus" ), 
    n_amostras     = 500,
    nivel_cobertura = 0.90,
    min_casos_treino = 30,
    besag_ativar  = FALSE,
    besag_niveis  = c("DRS"),
    besag_regioes = list(Estado = NULL, DRS = "TAUBATE", CIR = NULL, Municipio = NULL),
    block_ativar       = TRUE,
    block_besag_ativar = FALSE,
    block_min_treino   = 3
  ),
  
  largura = 1200, altura = 600, dpi = 150
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
instance_seed <- function(nivel, regiao, alvo, ano, tag = "") {
  key <- paste(nivel, regiao, alvo, ano, tag, sep = "|")
  as.integer(sum(utf8ToInt(key)) %% 2147483647L)
}

if (is.null(config$ano_atual) || !is.numeric(config$ano_atual) || is.na(config$ano_atual))
  stop("config$ano_atual must be a numeric year (e.g., 2026).", call. = FALSE)

# =============================================================================
#  2) PACKAGES
# =============================================================================
message("Checking packages...")
pacotes_cran <- c("dplyr", "tidyr", "stringr", "lubridate", "ggplot2")
for (p in pacotes_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
# INLA is only necessary for metodo_banda = "inla" and for validation.
precisa_inla <- identical(config$metodo_banda, "inla") || isTRUE(config$validacao$ativar) ||
  isTRUE(config$validacao$block_ativar) || isTRUE(config$validacao$besag_ativar) ||
  isTRUE(config$validacao$block_besag_ativar)
if (precisa_inla) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    message("Installing R-INLA (may take a while the first time)...")
    tryCatch(
      install.packages("INLA",
                       repos = c(CRAN = "https://cloud.r-project.org",
                                 INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE),
      error = function(e) message("Failed to install INLA: ", conditionMessage(e)))
  }
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("The INLA package is not available. Install it manually once, or use ",
         "metodo_banda = 'canal_geometrico' and disable validation.", call. = FALSE)
  suppressPackageStartupMessages(library(INLA))
}
message("Packages OK.\n")

# =============================================================================
#  3) DATABASE LOADING AND PREPARATION
# =============================================================================
message("Loading meningitis database...")

padronizar_mun <- function(x) {
  v <- as.character(x)
  v <- tryCatch(stringi::stri_trans_general(v, "Latin-ASCII"), error = function(e) v)
  trimws(toupper(v))
}
limpar_nome <- function(x) gsub("^_|_$", "", gsub("[^A-Z0-9]+", "_", padronizar_mun(x)))

carregar_base_meningite <- function(caminho) {
  if (!file.exists(caminho)) stop("File not found:\n  ", caminho, call. = FALSE)
  d <- readRDS(caminho)
  if (!is.data.frame(d)) stop("The RDS does not contain a data.frame.", call. = FALSE)
  message("  ", format(nrow(d), big.mark = ","), " rows.")
  as.data.frame(d)
}

# NOTE: Does NOT filter confirmed cases here (the filter is applied later, after the diagnostic).
preparar_base <- function(d) {
  tem <- function(x) x %in% names(d)
  fonte <- config$fonte_epi %||% "auto"
  
  usou <- NULL
  if (identical(fonte, "data_sintoma")) {
    if (!tem("data_sintoma")) stop("fonte_epi='data_sintoma' but the database lacks 'data_sintoma'.", call. = FALSE)
    d$data_sintoma <- as.Date(d$data_sintoma)
    d$ano_epi <- lubridate::epiyear(d$data_sintoma)
    d$sem_epi <- lubridate::epiweek(d$data_sintoma)
    usou <- "data_sintoma (epiyear/epiweek)"
  } else {  # "auto": prioritizes anoepi/sepi (matches the endemic channels code)
    if (tem("anoepi") && tem("sepi")) {
      d$ano_epi <- suppressWarnings(as.integer(as.character(d$anoepi)))
      d$sem_epi <- suppressWarnings(as.integer(as.character(d$sepi)))
      usou <- "anoepi/sepi (database columns)"
    } else if (tem("data_sintoma")) {
      d$data_sintoma <- as.Date(d$data_sintoma)
      d$ano_epi <- lubridate::epiyear(d$data_sintoma)
      d$sem_epi <- lubridate::epiweek(d$data_sintoma)
      usou <- "data_sintoma (epiyear/epiweek) [anoepi/sepi missing]"
    } else stop("Database missing 'anoepi'/'sepi' and 'data_sintoma'.", call. = FALSE)
  }
  message("  Epidemiological year/week source: ", usou, ".")
  n_na <- sum(is.na(d$ano_epi) | is.na(d$sem_epi))
  if (n_na > 0)
    message("  WARNING: ", format(n_na, big.mark = ","),
            " row(s) with missing year/week (will be ignored in the count).")
  
  # standardized municipality (name) and DRS (name)
  col_mun <- intersect(c("MUNI_RESID", "mun", "MUN"),  names(d))[1]
  col_drs <- intersect(c("DRS_RESID",  "drs", "DRS"),  names(d))[1]
  if (is.na(col_mun)) stop("Database missing municipality column (MUNI_RESID).", call. = FALSE)
  d$mun <- padronizar_mun(d[[col_mun]])
  if (!is.na(col_drs)) d$drs <- padronizar_mun(d[[col_drs]])
  
  for (cc in c("grupo_etiologico", "etiologia", "sorogrupo", "classi_fin", "evolucao"))
    if (!tem(cc)) d[[cc]] <- NA_character_
  d
}

base_menin <- carregar_base_meningite(config$caminho_dados)
base_menin <- preparar_base(base_menin)   # still WITHOUT confirmed filter
message("Database prepared (without confirmed filter yet).\n")

# =============================================================================
#  4) AUXILIARY FUNCTIONS (geography, target, aggregation)
# =============================================================================
listar_regioes <- function(nivel) {
  col <- switch(nivel, "DRS" = "drs", "Municipio" = "mun", "CIR" = "cir",
                stop("Invalid level: ", nivel))
  if (!col %in% names(base_menin)) stop("Column '", col, "' does not exist in the database.")
  sort(unique(na.omit(base_menin[[col]])))
}

# applies the filter for an epidemiological target (AND between defined fields).
# aplicar_confirmado = FALSE turns off the confirmed filter (used in diagnostics).
filtrar_alvo <- function(d, alvo, aplicar_confirmado = TRUE) {
  nome_alvo <- if (is.character(alvo)) alvo else NULL
  if (is.character(alvo)) alvo <- config$alvos[[alvo]]
  if (is.null(alvo)) alvo <- list(tipo = "todas")
  # etiological filters (skipped when tipo == "todas")
  if (!isTRUE(alvo$tipo == "todas")) {
    if (!is.null(alvo$grupo_etiologico))
      d <- d[!is.na(d$grupo_etiologico) & d$grupo_etiologico %in% alvo$grupo_etiologico, ]
    if (!is.null(alvo$etiologia))
      d <- d[!is.na(d$etiologia) & d$etiologia %in% alvo$etiologia, ]
    if (!is.null(alvo$sorogrupo))
      d <- d[!is.na(d$sorogrupo) & d$sorogrupo %in% alvo$sorogrupo, ]
  }
  # CONFIRMED filter (per target; default = config$apenas_confirmados)
  if (isTRUE(aplicar_confirmado)) {
    usar_conf <- if (!is.null(alvo$confirmado)) isTRUE(alvo$confirmado)
    else isTRUE(config$apenas_confirmados)
    if (usar_conf) d <- d[!is.na(d$classi_fin) & d$classi_fin == "Confirmado", ]
  }
  d
}

# weekly count (year, week, cases) for a region and a target
agregar_semanal <- function(base, nivel, regiao, alvo) {
  if (nivel == "Estado") d <- base
  else {
    col <- switch(nivel, "DRS" = "drs", "Municipio" = "mun", "CIR" = "cir")
    d   <- base[!is.na(base[[col]]) & base[[col]] == regiao, ]
  }
  d <- filtrar_alvo(d, alvo)
  d %>%
    filter(!is.na(ano_epi), !is.na(sem_epi)) %>%
    filter(sem_epi != 53) %>%      # removes EW53 (same as doc 7: remove_53 / SEPI != 53)
    group_by(ano_epi, sem_epi) %>%
    summarise(casos = n(), .groups = "drop") %>%
    rename(ano = ano_epi, semana = sem_epi)
}

# fills weeks 1:52 with zero for each training year
completar_semanas <- function(df, anos) {
  df <- df[df$ano %in% anos, , drop = FALSE]
  grade <- expand.grid(ano = anos, semana = 1:52)
  merge(grade, df, by = c("ano", "semana"), all.x = TRUE) %>%
    mutate(casos = ifelse(is.na(casos), 0L, casos)) %>%
    arrange(ano, semana)
}

# adds harmonic columns sin_k/cos_k (for sazonal = "harmonicos")
adicionar_harmonicos <- function(df, nh) {
  for (k in seq_len(nh)) {
    df[[paste0("sin", k)]] <- sin(2 * pi * k * df$semana / 52)
    df[[paste0("cos", k)]] <- cos(2 * pi * k * df$semana / 52)
  }
  df
}

# builds the INLA formula from the level spec (seasonal + family)
construir_formula <- function(modelo) {
  if (!is.null(modelo$formula)) return(modelo$formula)
  sazonal <- modelo$sazonal %||% "rw1"
  ano_term <- "f(ano_idx, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  if (identical(sazonal, "harmonicos")) {
    nh <- modelo$n_harmonicos %||% 2
    harm <- paste(sprintf("sin%d + cos%d", seq_len(nh), seq_len(nh)), collapse = " + ")
    stats::as.formula(paste("casos ~ 1 +", harm, "+", ano_term))
  } else if (identical(sazonal, "nenhum")) {
    stats::as.formula(paste("casos ~ 1 +", ano_term))
  } else {
    stats::as.formula(paste(
      "casos ~ 1 +",
      "f(semana, model = 'rw1', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
      ano_term))
  }
}
nh_do_modelo <- function(modelo)
  if (identical(modelo$sazonal %||% "rw1", "harmonicos")) (modelo$n_harmonicos %||% 2) else 0L

# =============================================================================
#  4.0.1) TARGET-YEAR WEEK CUTOFF
# =============================================================================
# Decides how many epidemiological weeks of config$ano_atual should be treated
# as "known" when building the observed trajectory overlay:
#   - if config$ano_atual IS the real current year (matches Sys.Date()), use the
#     current (partial) epidemiological week -> mirrors a live/running bulletin;
#   - otherwise (retrospective year, e.g. 2019), treat it as a CLOSED year and
#     use all 52 weeks -> completes the full historical series.
# Returns the LAST week considered "known" (0 if none).
semana_final_do_ano <- function(ano) {
  if (isTRUE(ano == lubridate::epiyear(Sys.Date()))) {
    max(lubridate::epiweek(Sys.Date()) - 1L, 0L)
  } else {
    52L
  }
}

# =============================================================================
#  4.1) CASE DIAGNOSTICS + application of the 'confirmed' filter
# =============================================================================
diagnostico_casos <- function(d) {
  message("\n================= CASE DIAGNOSTICS (sanity check) =================")
  message("Total rows in the database: ", format(nrow(d), big.mark = ","))
  
  # 1) classi_fin values (this is where 'Confirmado' might be written differently!)
  cf <- table(d$classi_fin, useNA = "ifany")
  message("\nValues of classi_fin (the filter requires EXACTLY 'Confirmado'):")
  for (nm in names(cf)) message("   '", nm, "' -> ", format(as.integer(cf[[nm]]), big.mark = ","))
  
  # 2) per year: ALL vs CONFIRMED, for some key targets
  alvos_diag <- intersect(c("Total Meningitis", "Meningococcus B", "Meningococcus C", "Meningococcus"),
                          names(config$alvos))
  for (al in alvos_diag) {
    sub  <- filtrar_alvo(d, al, aplicar_confirmado = FALSE)
    subc <- sub[!is.na(sub$classi_fin) & sub$classi_fin == "Confirmado", ]
    tt_all  <- tapply(rep(1L, nrow(sub)),  sub$ano_epi,  sum)
    tt_conf <- tapply(rep(1L, nrow(subc)), subc$ano_epi, sum)
    anos <- sort(unique(c(as.integer(names(tt_all)), as.integer(names(tt_conf)))))
    anos <- anos[!is.na(anos)]
    message("\n[", al, "] cases per year (ALL | CONFIRMED):")
    for (a in anos) {
      va <- if (as.character(a) %in% names(tt_all))  tt_all[[as.character(a)]]  else 0L
      vc <- if (as.character(a) %in% names(tt_conf)) tt_conf[[as.character(a)]] else 0L
      message("   ", a, ": ", va, " | ", vc)
    }
    message("   TOTAL: ", nrow(sub), " | ", nrow(subc))
  }
  message("\n(If 'CONFIRMED' is much lower than 'ALL', check the classi_fin filter.")
  message(" If any year disappears completely, check NA in anoepi/sepi/data_sintoma.)")
  message("==================================================================\n")
}

if (isTRUE(config$diagnostico_casos)) diagnostico_casos(base_menin)

# The confirmed filter is NO LONGER global: it is applied PER TARGET inside
# filtrar_alvo (target's 'confirmado' field, or the default config$apenas_confirmados).
# Thus, "Total Meningitis" comes out confirmed and B/C come out with all notifications.
message("Confirmed filter is now PER TARGET (default = ", config$apenas_confirmados,
        "; 'Total Meningitis' = confirmed).\n")

# =============================================================================
#  4.5) CIR — creates the 'cir' column in the base (spatial join municipality -> CIR)
# =============================================================================
obter_municipios_sp <- function() {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("The 'sf' package is required (install.packages('sf')).", call. = FALSE)
  if (!is.null(config$besag$shapefile_municipios)) {
    mp  <- sf::st_read(config$besag$shapefile_municipios, quiet = TRUE)
    idc <- config$besag$coluna_id_municipios
    if (is.null(idc)) {
      gcol  <- attr(mp, "sf_column")
      cands <- setdiff(grep("nm_mun|name_muni|nome|municip|^mun", names(mp),
                            ignore.case = TRUE, value = TRUE), gcol)
      if (!length(cands))
        stop("Could not detect the municipality name column; define config$besag$coluna_id_municipios.", call. = FALSE)
      idc <- cands[1]
      message("  [Muni] name column detected: '", idc, "'.")
    }
    if (!idc %in% names(mp)) stop("Column '", idc, "' does not exist in the municipalities shapefile.", call. = FALSE)
    mp$chave_mun <- padronizar_mun(mp[[idc]])
    return(mp["chave_mun"])
  }
  if (isTRUE(config$besag$usar_geobr)) {
    if (!requireNamespace("geobr", quietly = TRUE))
      install.packages("geobr", repos = "https://cloud.r-project.org")
    ano_geo <- config$besag$ano_geobr %||% 2020
    mp <- geobr::read_municipality(code_muni = "SP", year = ano_geo, showProgress = FALSE)
    mp$chave_mun <- padronizar_mun(mp$name_muni)
    return(mp["chave_mun"])
  }
  stop("No polygon source: define shapefile_municipios or usar_geobr = TRUE.", call. = FALSE)
}

preparar_cir_na_base <- function(base) {
  if ("cir" %in% names(base)) return(base)
  if (!is.null(config$cir$lookup)) {
    lk <- config$cir$lookup
    if (is.character(lk) && length(lk) == 1) lk <- utils::read.csv(lk, stringsAsFactors = FALSE)
    names(lk) <- tolower(names(lk))
    if (!all(c("chave_mun", "cir") %in% names(lk)))
      stop("config$cir$lookup must have 'chave_mun' and 'cir' columns.", call. = FALSE)
    lk$chave_mun <- padronizar_mun(lk$chave_mun)
    base$cir <- lk$cir[match(base$mun, lk$chave_mun)]
    return(base)
  }
  if (!requireNamespace("sf", quietly = TRUE))
    stop("The 'sf' package is required to map municipalities to CIR.", call. = FALSE)
  if (is.null(config$cir$shapefile) || !file.exists(config$cir$shapefile))
    stop("CIR shapefile not found: ", config$cir$shapefile, call. = FALSE)
  suppressWarnings(try(sf::sf_use_s2(FALSE), silent = TRUE))
  cir_sf <- sf::st_make_valid(sf::st_read(config$cir$shapefile, quiet = TRUE))
  gcol <- attr(cir_sf, "sf_column")
  idc  <- config$cir$coluna_id
  if (is.null(idc)) {
    cands <- setdiff(grep("cir|nome|nom|name|regi", names(cir_sf), ignore.case = TRUE, value = TRUE), gcol)
    idc <- if (length(cands)) cands[1] else setdiff(names(cir_sf), gcol)[1]
    message("  [CIR] identifier column detected: '", idc, "'.")
  }
  if (!idc %in% names(cir_sf)) stop("Column '", idc, "' does not exist in the CIR shapefile.", call. = FALSE)
  cir_sf[[idc]] <- padronizar_mun(cir_sf[[idc]])
  munpoly <- sf::st_make_valid(obter_municipios_sp())
  cir_sf  <- sf::st_transform(cir_sf, sf::st_crs(munpoly))
  cent <- suppressWarnings(sf::st_point_on_surface(munpoly))
  jn   <- suppressWarnings(sf::st_join(cent, cir_sf[idc], join = sf::st_within))
  mapa <- data.frame(chave_mun = as.character(jn$chave_mun),
                     cir = as.character(jn[[idc]]), stringsAsFactors = FALSE)
  mapa <- mapa[!is.na(mapa$cir) & !duplicated(mapa$chave_mun), ]
  n_sem <- sum(is.na(match(unique(base$mun), mapa$chave_mun)))
  if (n_sem > 0) message("  [CIR] WARNING: ", n_sem, " municipality(ies) without CIR (will be NA).")
  base$cir <- mapa$cir[match(base$mun, mapa$chave_mun)]
  message("  [CIR] ", length(unique(na.omit(base$cir))), " CIR(s) assigned.")
  base
}

precisa_cir <- ("CIR" %in% config$niveis) ||
  ("CIR" %in% config$validacao$besag_niveis) ||
  (length(config$validacao$cir) > 0)
if (isTRUE(precisa_cir)) {
  message("Preparing CIR in the database (municipality -> CIR)...")
  base_menin <- preparar_cir_na_base(base_menin)
}

# =============================================================================
#  5) ENGINE (A) — GEOMETRIC ENDEMIC CHANNEL  (MIRRORS endemic.channels)
# =============================================================================
# Reproduces EXACTLY endemic.channels::endemic_channel(type='geometric') (Bortman 1999),
# verified in the package source:
#   inc = CUMULATIVE per year; per epidemiological week:
#     log_inc  = log(inc + 1)
#     log_mean  = mean(log_inc)
#     log_upper = log_mean + qt(0.975, n_anos-1) * sd(log_inc) / sqrt(n_anos)
#     log_lower = log_mean - qt(0.975, n_anos-1) * sd(log_inc) / sqrt(n_anos)
#     mean/upper/lower = exp(log_...) - 1      (the package DOES NOT truncate at 0)
# Zones (package): Below expected [0,lower] | Expected [lower,mean] |
#                  Alert [mean,upper] | Epidemic [upper, Inf].
banda_canal_geometrico <- function(df_treino, anos, conf = 0.95,
                                   tipo = "ic_media", offset = 1, min_casos = 30) {
  total <- sum(df_treino$casos)
  if (total < min_casos) {
    message("    [WARNING] Cases in training (", total, ") < minimum (", min_casos, "); skipping.")
    return(NULL)
  }
  anos <- sort(unique(anos)); N <- length(anos)   # n_anos (same as package)
  if (N < 2) { message("    [WARNING] < 2 years; skipping."); return(NULL) }
  
  # 52 x N matrix with the CUMULATIVE of each year (df_treino comes complete in 1..52)
  M <- sapply(anos, function(yr) {
    v <- df_treino$casos[df_treino$ano == yr]
    if (length(v) != 52) return(rep(NA_real_, 52))
    cumsum(v)
  })
  if (is.null(dim(M))) return(NULL)
  
  tq <- stats::qt(1 - (1 - conf) / 2, df = N - 1)   # package: qt(0.975, n_anos-1)
  res <- t(apply(M, 1, function(x) {
    L <- log(x + offset)                            # package: log(inc + 1)
    m <- mean(L); s <- stats::sd(L)                 # sample sd (denom. n-1)
    if (!is.finite(s)) s <- 0
    half <- if (identical(tipo, "preditivo")) tq * s else tq * s / sqrt(N)
    c(media = exp(m) - offset,                      # package: exp(log_mean) - 1
      lower = exp(m - half) - offset,               # WITHOUT truncating at 0 (same as package)
      upper = exp(m + half) - offset)
  }))
  data.frame(semana = 1:52, media = res[, "media"],
             lower = res[, "lower"], upper = res[, "upper"])
}

# =============================================================================
#  5.1) ENGINE (B) — OPTION A (INLA)  [only used if metodo_banda = "inla"]
# =============================================================================
simular_acumulado <- function(amostras, linhas_alvo, idx_zip = NA_integer_) {
  vapply(amostras, function(s) {
    mu <- exp(s$latent[linhas_alvo, 1]); y <- rpois(length(mu), mu)
    if (!is.na(idx_zip)) {
      p <- as.numeric(s$hyperpar[idx_zip])
      if (is.finite(p) && p > 0) y[runif(length(y)) < p] <- 0L
    }
    cumsum(y)
  }, numeric(length(linhas_alvo)))
}
idx_zip_hyper <- function(amostras, familia) {
  if (!grepl("zeroinflated", familia, ignore.case = TRUE)) return(NA_integer_)
  nomes <- names(amostras[[1]]$hyperpar); cand <- grep("zero|prob", nomes, ignore.case = TRUE)
  if (length(cand) >= 1) cand[1] else NA_integer_
}
rodar_inla_acumulado <- function(df_treino, n_amostras = 1000, min_casos = 30, modelo = NULL, seed = NULL) {
  total_casos <- sum(df_treino$casos)
  if (total_casos < min_casos) {
    message("    [WARNING] Cases in training (", total_casos, ") < minimum (", min_casos, "); skipping."); return(NULL)
  }
  anos_u <- sort(unique(df_treino$ano))
  if (length(anos_u) < 3) { message("    [WARNING] Less than 3 years; skipping."); return(NULL) }
  df_treino$ano_idx <- match(df_treino$ano, anos_u)
  if (is.null(modelo)) modelo <- list()
  familia <- modelo$familia %||% "poisson"; formula_inla <- construir_formula(modelo); nh <- nh_do_modelo(modelo)
  idx_ano_novo <- length(anos_u) + 1L
  df_pred <- data.frame(casos = NA_integer_, semana = 1:52, ano_idx = idx_ano_novo)
  df_full <- rbind(df_treino[, c("casos", "semana", "ano_idx")], df_pred)
  if (nh > 0) df_full <- adicionar_harmonicos(df_full, nh)
  fit <- tryCatch(
    inla(formula_inla, family = familia, data = df_full,
         control.compute = list(config = TRUE),
         control.predictor = list(link = 1, compute = TRUE), verbose = FALSE),
    error = function(e) { message("    [INLA ERROR] ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  idx_pred <- (nrow(df_treino) + 1):nrow(df_full)
  if (!is.null(seed)) set.seed(seed)
  amostras <- inla.posterior.sample(n_amostras, fit)
  pred_rows <- grep("^Predictor", rownames(amostras[[1]]$latent))
  linhas_alvo <- pred_rows[idx_pred]; idx_zip <- idx_zip_hyper(amostras, familia)
  acum_mat <- simular_acumulado(amostras, linhas_alvo, idx_zip)
  # Channel ZONES from the POSTERIOR: same percentiles as the endemic channel.
  qs <- config$inla_percentis %||% c(lower = 0.025, media = 0.50, upper = 0.975)
  data.frame(semana = 1:52,
             lower = apply(acum_mat, 1, quantile, qs[["lower"]], na.rm = TRUE),
             media = apply(acum_mat, 1, quantile, qs[["media"]], na.rm = TRUE),
             upper = apply(acum_mat, 1, quantile, qs[["upper"]], na.rm = TRUE))
}

# =============================================================================
#  6) PLOTS
# =============================================================================
tema_banda <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title         = element_text(size = base_size + 1),
      axis.text.x        = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y        = element_text(size = base_size),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = 7.5),
      legend.key.size    = unit(0.8, "lines"),
      legend.box.spacing = unit(2, "pt"),
      legend.box.margin  = margin(t = -4),
      plot.title         = element_text(size = base_size + 3, hjust = 0.5, face = "bold"),
      plot.subtitle      = element_text(size = base_size,     hjust = 0.5, color = "grey40"),
      plot.margin        = margin(8, 12, 8, 8)
    )
}

# --- 4-ZONE CHANNEL FIGURE (used for both INLA and geometric channel) ---
# Receives a data.frame with columns: semana, lower, media, upper.
# Draws the 4 zones (Below/Expected/Alert/Epidemic) with the SAME colors as
# the endemic.channels package, the lower/upper lines, and the center line (media),
# and overlays the cumulative observed data for the current year.
grafico_canal_zonas <- function(canal, obs_acum, titulo, subtitulo, lab_obs,
                                lab_centro = "Median") {
  ymax <- suppressWarnings(max(c(canal$upper, obs_acum$acum), na.rm = TRUE)) * 1.10
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(c(canal$upper, 1), na.rm = TRUE)
  ymax <- max(ymax, 1); canal$teto <- ymax
  
  rot <- c("Below expected", "Expected", "Alert", "Epidemic")
  if (isTRUE(config$usar_carinhas))
    rot <- paste(config$carinhas_canal[c("below", "expected", "alert", "epidemic")], rot)
  
  g <- ggplot(canal, aes(x = semana)) +
    geom_ribbon(aes(ymin = upper, ymax = teto,  fill = "Epidemic")) +
    geom_ribbon(aes(ymin = media, ymax = upper, fill = "Alert")) +
    geom_ribbon(aes(ymin = lower, ymax = media, fill = "Expected")) +
    geom_ribbon(aes(ymin = 0,     ymax = lower, fill = "Below expected")) +
    geom_line(aes(y = lower), color = "black", linewidth = 0.3) +
    geom_line(aes(y = upper), color = "black", linewidth = 0.3) +
    geom_line(aes(y = media, color = lab_centro), linewidth = 0.7)
  
  if (!is.null(obs_acum) && nrow(obs_acum) > 0) {
    g <- g +
      geom_line(data = obs_acum, aes(x = semana, y = acum, color = lab_obs), linewidth = 0.8) +
      geom_point(data = obs_acum, aes(x = semana, y = acum, color = lab_obs), size = 2.0)
  }
  
  cores_lin <- setNames("#1565C0", lab_centro); cores_lin[lab_obs] <- "#000000"
  
  g +
    scale_fill_manual(
      values = c("Below expected" = "#0bc82c", "Expected" = "#edea44",
                 "Alert" = "#fe3d3d", "Epidemic" = "#ba0606"),   # exact package colors
      breaks = c("Below expected", "Expected", "Alert", "Epidemic"),
      labels = rot) +
    scale_color_manual(values = cores_lin) +
    scale_x_continuous(breaks = 1:52, labels = function(x) paste0("EW", sprintf("%02d", x)),
                       limits = c(1, 52), expand = expansion(add = 0)) +
    scale_y_continuous(breaks = function(l) pretty(l, n = 6),
                       expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(title = titulo, subtitle = subtitulo,
         x = "Epidemiological week", y = "Cumulative cases in the year") +
    guides(fill = guide_legend(order = 1, nrow = 1), color = guide_legend(order = 2, nrow = 1)) +
    tema_banda()
}

# =============================================================================
#  7) SAVE FIGURE
# =============================================================================
salvar_banda <- function(plt, nome_arquivo, pasta_regiao) {
  if (!dir.exists(pasta_regiao)) dir.create(pasta_regiao, recursive = TRUE)
  arq <- file.path(pasta_regiao, paste0(nome_arquivo, ".png"))
  ggsave(arq, plt, width = config$largura / config$dpi, height = config$altura / config$dpi,
         units = "in", dpi = config$dpi, device = "png", bg = "white")
  message("    -> saved: ", arq); invisible(arq)
}

# =============================================================================
#  8) MAIN ENGINE — band + figure for ONE region + ONE target
# =============================================================================
processar_regiao_alvo <- function(nivel, regiao, alvo_nome) {
  rot_reg <- if (nivel == "Estado") "SP_State" else limpar_nome(regiao)
  rot_alv <- limpar_nome(alvo_nome)
  message("  [", nivel, "] ", regiao, " | ", alvo_nome)
  
  df_todos  <- agregar_semanal(base_menin, nivel, regiao, alvo_nome)
  df_treino <- df_todos %>% filter(ano %in% config$anos_treino) %>% completar_semanas(config$anos_treino)
  
  # CUMULATIVE observed for the target year (full 52 weeks if retrospective; partial if current)
  sem_hoje <- semana_final_do_ano(config$ano_atual) + 1L
  obs_atual <- df_todos %>% filter(ano == config$ano_atual, semana < sem_hoje) %>% select(semana, casos)
  if (sem_hoje > 1) {
    obs_atual <- data.frame(semana = 1:(sem_hoje - 1)) %>%
      left_join(obs_atual, by = "semana") %>%
      mutate(casos = tidyr::replace_na(casos, 0L)) %>% arrange(semana) %>% mutate(acum = cumsum(casos))
  } else obs_atual <- data.frame(semana = integer(0), casos = integer(0), acum = numeric(0))
  
  titulo_reg <- if (nivel == "Estado") "State of São Paulo" else paste0(nivel, ": ", regiao)
  titulo_alv <- paste0("Meningitis (cumulative) — ", alvo_nome)
  lab_obs    <- paste0("Cumulative observed ", config$ano_atual)
  metodo     <- config$metodo_banda %||% "canal_geometrico"
  
  if (identical(metodo, "canal_geometrico")) {
    message("    Geometric endemic channel...")
    canal <- banda_canal_geometrico(df_treino, config$anos_treino,
                                    conf   = config$canal$nivel_confianca %||% 0.95,
                                    tipo   = config$canal$tipo_intervalo  %||% "ic_media",
                                    offset = config$canal$offset_log      %||% 1,
                                    min_casos = config$min_casos_treino)
    if (is.null(canal)) return(invisible(NULL))
    subt <- paste0("Geometric endemic channel (endemic.channels, Bortman; CI ",
                   round(100 * (config$canal$nivel_confianca %||% 0.95)), "%) | Training: ",
                   paste(range(config$anos_treino), collapse = "-"), " (excl. 2020-2022)")
    plt <- grafico_canal_zonas(canal, if (nrow(obs_atual) > 0) obs_atual else NULL,
                               titulo_alv, subt, lab_obs,
                               lab_centro = "Geometric mean")
  } else {
    message("    Fitting INLA model (Weekly Poisson -> cumulative)...")
    canal <- rodar_inla_acumulado(df_treino, config$n_amostras, config$min_casos_treino,
                                  config$modelos_inla[[nivel]],
                                  seed = instance_seed(nivel, rot_reg, alvo_nome, config$ano_atual, "banda"))
    if (is.null(canal)) return(invisible(NULL))
    qs <- config$inla_percentis %||% c(lower = 0.025, media = 0.50, upper = 0.975)
    subt <- paste0(" ")
    plt <- grafico_canal_zonas(canal, if (nrow(obs_atual) > 0) obs_atual else NULL,
                               titulo_alv, subt, lab_obs,
                               lab_centro = "Median")
  }
  
  pasta_regiao <- file.path(config$pasta_saida, paste0("bandas_", nivel), rot_reg)
  salvar_banda(plt, paste0("banda_", tolower(nivel), "_", rot_reg, "_", rot_alv), pasta_regiao)
}

# =============================================================================
#  8.5) BESAG — JOINT SPATIAL FIT (only "inla" method)
# =============================================================================
grafo_de_shapefile <- function(nivel, path, coluna_id) {
  for (pk in c("sf", "spdep")) if (!requireNamespace(pk, quietly = TRUE))
    stop("The '", pk, "' package is required for Besag.", call. = FALSE)
  suppressWarnings(try(sf::sf_use_s2(FALSE), silent = TRUE))
  mapa <- sf::st_make_valid(sf::st_read(path, quiet = TRUE))
  gcol <- attr(mapa, "sf_column")
  if (is.null(coluna_id)) {
    cands <- setdiff(grep("drs|cir|nome|nom|name|regi", names(mapa), ignore.case = TRUE, value = TRUE), gcol)
    coluna_id <- if (length(cands)) cands[1] else setdiff(names(mapa), gcol)[1]
    message("  [Besag ", nivel, "] region column detected: '", coluna_id, "'.")
  }
  if (!coluna_id %in% names(mapa)) stop("Column '", coluna_id, "' does not exist in the ", nivel, " shapefile.", call. = FALSE)
  mapa$regiao <- padronizar_mun(mapa[[coluna_id]])
  regpoly <- aggregate(mapa[gcol], by = list(regiao = mapa$regiao), FUN = function(x) x[1], do_union = TRUE)
  regpoly <- regpoly[order(regpoly$regiao), ]
  regpoly <- suppressWarnings(sf::st_buffer(sf::st_make_valid(regpoly), 0))
  viz <- spdep::poly2nb(regpoly); n_ilhas <- sum(spdep::card(viz) == 0)
  if (n_ilhas > 0) message("  [Besag ", nivel, "] WARNING: ", n_ilhas, " region(s) without neighbors.")
  arq <- tempfile(fileext = ".adj"); spdep::nb2INLA(arq, viz); G <- INLA::inla.read.graph(arq)
  regioes <- as.character(regpoly$regiao)
  reg_base <- tryCatch(listar_regioes(nivel), error = function(e) character(0))
  falta <- setdiff(reg_base, regioes)
  if (length(falta)) message("  [Besag ", nivel, "] ", length(falta), " region(s) in the database without polygon (e.g.: ",
                             paste(utils::head(falta, 3), collapse = ", "), ").")
  message("  [Besag ", nivel, "] SHAPEFILE graph: ", length(regioes), " regions.")
  list(grafo = G, regioes = regioes)
}
preparar_grafo_besag <- function(nivel) {
  for (pk in c("sf", "spdep")) if (!requireNamespace(pk, quietly = TRUE))
    stop("The '", pk, "' package is required for Besag.", call. = FALSE)
  if (!is.null(config$grafos[[nivel]]) && !is.null(config$regioes_grafo[[nivel]]))
    return(list(grafo = config$grafos[[nivel]], regioes = config$regioes_grafo[[nivel]]))
  shp_path <- switch(nivel, "DRS" = config$besag$shapefile_drs, "CIR" = config$besag$shapefile_cir, NULL)
  shp_id   <- switch(nivel, "DRS" = config$besag$coluna_id_drs, "CIR" = config$besag$coluna_id_cir, NULL)
  if (!is.null(shp_path)) return(grafo_de_shapefile(nivel, shp_path, shp_id))
  munpoly <- obter_municipios_sp()
  col <- switch(nivel, "DRS" = "drs", "CIR" = "cir", "Municipio" = "mun", stop("Besag does not apply to level ", nivel))
  mapa <- unique(base_menin[!is.na(base_menin[[col]]) & !is.na(base_menin$mun), c("mun", col)])
  names(mapa) <- c("chave_mun", "regiao"); mapa$chave_mun <- as.character(mapa$chave_mun)
  munpoly <- merge(munpoly, mapa, by = "chave_mun")
  if (nrow(munpoly) == 0) stop("Could not match database municipalities with polygons.", call. = FALSE)
  suppressWarnings(try(sf::sf_use_s2(FALSE), silent = TRUE)); munpoly <- sf::st_make_valid(munpoly)
  regpoly <- aggregate(munpoly["chave_mun"], by = list(regiao = munpoly$regiao), FUN = function(x) x[1], do_union = TRUE)
  regpoly <- regpoly[order(regpoly$regiao), ]
  regpoly <- suppressWarnings(sf::st_buffer(sf::st_make_valid(regpoly), 0))
  viz <- spdep::poly2nb(regpoly); n_ilhas <- sum(spdep::card(viz) == 0)
  if (n_ilhas > 0) message("  [Besag] WARNING: ", n_ilhas, " region(s) without neighbors (islands).")
  arq <- tempfile(fileext = ".adj"); spdep::nb2INLA(arq, viz); G <- INLA::inla.read.graph(arq)
  regioes <- as.character(regpoly$regiao)
  message("  [Besag] level ", nivel, ": graph with ", length(regioes), " regions.")
  list(grafo = G, regioes = regioes)
}
construir_formula_besag <- function(modelo, G) {
  if (!is.null(modelo$formula)) { f <- modelo$formula; environment(f) <- environment(); return(f) }
  sazonal <- modelo$sazonal %||% "rw1"
  ano_term <- "f(ano_idx, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  esp_term <- "f(regiao_idx, model = 'besag', graph = G, scale.model = TRUE, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  if (identical(sazonal, "harmonicos")) {
    nh <- modelo$n_harmonicos %||% 2
    harm <- paste(sprintf("sin%d + cos%d", seq_len(nh), seq_len(nh)), collapse = " + ")
    f <- stats::as.formula(paste("casos ~ 1 +", harm, "+", ano_term, "+", esp_term))
  } else if (identical(sazonal, "nenhum")) {
    f <- stats::as.formula(paste("casos ~ 1 +", ano_term, "+", esp_term))
  } else {
    f <- stats::as.formula(paste("casos ~ 1 +",
                                 "f(semana, model = 'rw1', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
                                 ano_term, "+", esp_term))
  }
  environment(f) <- environment(); f
}
modelo_besag_efetivo <- function(nivel) {
  base_m <- config$modelos_inla[[nivel]] %||% list(); ovr <- config$besag
  list(familia = ovr$familia %||% (base_m$familia %||% "poisson"),
       sazonal = ovr$sazonal %||% (base_m$sazonal %||% "rw1"),
       n_harmonicos = ovr$n_harmonicos %||% (base_m$n_harmonicos %||% 2),
       formula = ovr$formula)
}
processar_nivel_besag <- function(nivel, alvos_rodar) {
  message("\n>>> [BESAG] JOINT SPATIAL fit for level ", nivel, " ...")
  gb <- preparar_grafo_besag(nivel); regioes <- gb$regioes; G <- gb$grafo
  if (length(regioes) < 2) { message("  [Besag] < 2 regions; skipping."); return(invisible(NULL)) }
  modelo <- modelo_besag_efetivo(nivel); familia <- modelo$familia %||% "poisson"; nh <- nh_do_modelo(modelo)
  for (alvo_nome in alvos_rodar) {
    message("  [BESAG ", nivel, "] target: ", alvo_nome)
    lst <- lapply(regioes, function(r) {
      cont <- agregar_semanal(base_menin, nivel, r, alvo_nome)
      tr <- completar_semanas(cont[cont$ano %in% config$anos_treino, , drop = FALSE], config$anos_treino)
      tr$regiao <- r; tr
    })
    treino <- do.call(rbind, lst)
    if (sum(treino$casos) < config$min_casos_treino) { message("    [WARNING] few cases in level; skipping."); next }
    anos_u <- sort(unique(treino$ano)); treino$ano_idx <- match(treino$ano, anos_u)
    treino$regiao_idx <- match(treino$regiao, regioes); idx_ano_novo <- length(anos_u) + 1L
    pred <- expand.grid(semana = 1:52, regiao = regioes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    pred$casos <- NA_integer_; pred$ano_idx <- idx_ano_novo; pred$regiao_idx <- match(pred$regiao, regioes)
    dados <- rbind(treino[, c("casos", "semana", "ano_idx", "regiao_idx")], pred[, c("casos", "semana", "ano_idx", "regiao_idx")])
    if (nh > 0) dados <- adicionar_harmonicos(dados, nh)
    formula_besag <- construir_formula_besag(modelo, G)
    fit <- tryCatch(inla(formula_besag, family = familia, data = dados,
                         control.compute = list(config = TRUE), control.predictor = list(link = 1, compute = TRUE),
                         verbose = FALSE), error = function(e) { message("    [INLA ERROR besag] ", e$message); NULL })
    if (is.null(fit)) next
    set.seed(instance_seed(nivel, "ALL", alvo_nome, config$ano_atual, "besag"))
    amostras <- inla.posterior.sample(config$n_amostras, fit)
    pred_rows <- grep("^Predictor", rownames(amostras[[1]]$latent)); idx_zip <- idx_zip_hyper(amostras, familia)
    n_tr <- nrow(treino)
    for (r in regioes) {
      pos_pred <- which(pred$regiao == r); linhas <- pred_rows[n_tr + pos_pred]
      acum_mat <- simular_acumulado(amostras, linhas, idx_zip)
      qs <- config$inla_percentis %||% c(lower = 0.025, media = 0.50, upper = 0.975)
      canal <- data.frame(semana = 1:52,
                          lower = apply(acum_mat, 1, quantile, qs[["lower"]], na.rm = TRUE),
                          media = apply(acum_mat, 1, quantile, qs[["media"]], na.rm = TRUE),
                          upper = apply(acum_mat, 1, quantile, qs[["upper"]], na.rm = TRUE))
      df_todos_r <- agregar_semanal(base_menin, nivel, r, alvo_nome)
      sem_hoje <- semana_final_do_ano(config$ano_atual) + 1L
      obs_atual <- df_todos_r %>% filter(ano == config$ano_atual, semana < sem_hoje) %>% select(semana, casos)
      if (sem_hoje > 1) {
        obs_atual <- data.frame(semana = 1:(sem_hoje - 1)) %>% left_join(obs_atual, by = "semana") %>%
          mutate(casos = tidyr::replace_na(casos, 0L)) %>% arrange(semana) %>% mutate(acum = cumsum(casos))
      } else obs_atual <- data.frame(semana = integer(0), acum = numeric(0))
      rot_reg <- limpar_nome(r); rot_alv <- limpar_nome(alvo_nome)
      plt <- grafico_canal_zonas(canal, if (nrow(obs_atual) > 0) obs_atual else NULL,
                                 paste0("Meningitis (cumulative) — ", alvo_nome),
                                 paste0("SPATIAL INLA (Besag); zones = posterior percentiles | Training: ",
                                        paste(range(config$anos_treino), collapse = "-"), " (excl. 2020-2022)"),
                                 paste0("Cumulative observed ", config$ano_atual), lab_centro = "Median")
      pasta_regiao <- file.path(config$pasta_saida, paste0("bandas_", nivel), rot_reg)
      salvar_banda(plt, paste0("banda_", tolower(nivel), "_", rot_reg, "_", rot_alv), pasta_regiao)
    }
  }
}

# =============================================================================
#  8.6) RECONSTRUCTION OF THE STATE FROM A GRANULAR LEVEL (DRS/CIR/Municipality)
# =============================================================================
#  PRINCIPLE: fits the model at the GRANULAR level (one region at a time, or all together via
#  Besag) and sums the cumulative trajectories SAMPLE BY SAMPLE from the posterior.
#    - summing by sample   -> correctly propagates uncertainty  ✔
#    - summing quantiles   -> WRONG (invalid)
#  Since cumsum is linear, summing cumulatives = cumulating the weekly sum.
#
#  Without Besag: fits are independent; when summing sample i from each region,
#  it assumes INDEPENDENCE between regions (the band tends to be NARROWER
#  than the directly fitted State band, due to error cancellation).
#  With Besag: all regions come from the SAME joint fit, so spatial correlation
#  is already embedded in the sample.
# =============================================================================

# complete grid (52 weeks x years) of counts per region of the granular level
grade_grao <- function(nivel_grao, alvo_nome, anos, regioes) {
  do.call(rbind, lapply(regioes, function(r) {
    cont <- agregar_semanal(base_menin, nivel_grao, r, alvo_nome)
    tr <- completar_semanas(cont, anos)
    tr$regiao <- r
    tr
  }))
}

# list of regions in the granular level (uses shapefile for Besag, to match graph)
regioes_do_grao <- function(nivel_grao, usar_besag) {
  if (isTRUE(usar_besag)) preparar_grafo_besag(nivel_grao)$regioes
  else listar_regioes(nivel_grao)
}

# --- (A) INDEPENDENT fit of ONE region -> cumulative matrix (n_weeks x n_samples)
matriz_acum_indep <- function(df_treino, semanas_pred, modelo = NULL,
                              n_amostras = 1000, min_casos = 0, seed = NULL) {
  if (sum(df_treino$casos) < min_casos) return(NULL)
  anos_u <- sort(unique(df_treino$ano))
  if (length(anos_u) < 3) return(NULL)
  df_treino$ano_idx <- match(df_treino$ano, anos_u)
  if (is.null(modelo)) modelo <- list()
  familia <- modelo$familia %||% "poisson"
  formula_inla <- construir_formula(modelo); nh <- nh_do_modelo(modelo)
  idx_novo <- length(anos_u) + 1L       # target year = NEW index (drawn year effect)
  df_pred <- data.frame(casos = NA_integer_, semana = semanas_pred, ano_idx = idx_novo)
  df_full <- rbind(df_treino[, c("casos", "semana", "ano_idx")], df_pred)
  if (nh > 0) df_full <- adicionar_harmonicos(df_full, nh)
  fit <- tryCatch(
    inla(formula_inla, family = familia, data = df_full,
         control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
         control.predictor = list(link = 1, compute = TRUE), verbose = FALSE),
    error = function(e) { message("      [INLA ERROR region] ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  idx_pred  <- (nrow(df_treino) + 1):nrow(df_full)
  if (!is.null(seed)) set.seed(seed)
  amostras  <- inla.posterior.sample(n_amostras, fit)
  pred_rows <- grep("^Predictor", rownames(amostras[[1]]$latent))
  m <- simular_acumulado(amostras, pred_rows[idx_pred], idx_zip_hyper(amostras, familia))
  if (is.null(dim(m))) m <- matrix(m, nrow = length(semanas_pred))
  list(acum_mat = m,
       dic  = if (!is.null(fit$dic$dic))   fit$dic$dic   else NA_real_,
       waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_)
}

# --- (B) sums SAMPLE BY SAMPLE across regions, with INDEPENDENT fits
soma_grao_indep <- function(nivel_grao, alvo_nome, regioes, anos_treino,
                            semanas_pred, n_amostras, min_casos_regiao = 0) {
  modelo <- config$modelos_inla[[nivel_grao]]
  soma <- NULL; n_ok <- 0L; falhas <- character(0); dics <- c(); waics <- c()
  for (r in regioes) {
    cont <- agregar_semanal(base_menin, nivel_grao, r, alvo_nome)
    tr   <- completar_semanas(cont, anos_treino)
    res  <- matriz_acum_indep(tr, semanas_pred, modelo, n_amostras, min_casos_regiao,
                              seed = instance_seed(nivel_grao, r, alvo_nome,
                                                   paste(range(anos_treino), collapse = "-"), "recon_indep"))
    if (is.null(res)) { falhas <- c(falhas, r); next }
    soma <- if (is.null(soma)) res$acum_mat else soma + res$acum_mat
    dics <- c(dics, res$dic); waics <- c(waics, res$waic); n_ok <- n_ok + 1L
  }
  if (is.null(soma)) return(NULL)
  if (length(falhas))
    message("    [Reconstruction] WARNING: ", length(falhas), " region(s) without fit ",
            "(entered as ZERO -> total UNDERESTIMATED): ",
            paste(utils::head(falhas, 5), collapse = ", "))
  list(acum_mat = soma, n_regioes = n_ok, falhas = falhas,
       dic = mean(dics, na.rm = TRUE), waic = mean(waics, na.rm = TRUE))
}

# --- (C) sums SAMPLE BY SAMPLE across regions, with JOINT fit (Besag)
soma_grao_besag <- function(nivel_grao, alvo_nome, anos_treino, semanas_pred, n_amostras) {
  gb <- preparar_grafo_besag(nivel_grao); regioes <- gb$regioes; G <- gb$grafo
  if (length(regioes) < 2) { message("    [Besag Reconstruction] < 2 regions; skipping."); return(NULL) }
  modelo  <- modelo_besag_efetivo(nivel_grao)
  familia <- modelo$familia %||% "poisson"; nh <- nh_do_modelo(modelo)
  
  treino <- grade_grao(nivel_grao, alvo_nome, anos_treino, regioes)
  if (is.null(treino) || !nrow(treino)) return(NULL)
  anos_u <- sort(unique(treino$ano))
  if (length(anos_u) < 3) return(NULL)
  treino$ano_idx    <- match(treino$ano, anos_u)
  treino$regiao_idx <- match(treino$regiao, regioes)
  
  pred <- expand.grid(semana = semanas_pred, regiao = regioes,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  pred$casos <- NA_integer_; pred$ano_idx <- length(anos_u) + 1L
  pred$regiao_idx <- match(pred$regiao, regioes)
  
  dados <- rbind(treino[, c("casos", "semana", "ano_idx", "regiao_idx")],
                 pred[,   c("casos", "semana", "ano_idx", "regiao_idx")])
  if (nh > 0) dados <- adicionar_harmonicos(dados, nh)
  formula_besag <- construir_formula_besag(modelo, G)
  
  fit <- tryCatch(
    inla(formula_besag, family = familia, data = dados,
         control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
         control.predictor = list(link = 1, compute = TRUE), verbose = FALSE),
    error = function(e) { message("    [INLA ERROR besag reconstruction] ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  
  set.seed(instance_seed(nivel_grao, "ALL", alvo_nome, paste(range(anos_treino), collapse = "-"), "recon_besag"))
  amostras  <- inla.posterior.sample(n_amostras, fit)
  pred_rows <- grep("^Predictor", rownames(amostras[[1]]$latent))
  idx_zip   <- idx_zip_hyper(amostras, familia); n_tr <- nrow(treino)
  
  soma <- NULL
  for (r in regioes) {
    pos <- which(pred$regiao == r)
    pos <- pos[order(pred$semana[pos])]
    m <- simular_acumulado(amostras, pred_rows[n_tr + pos], idx_zip)
    if (is.null(dim(m))) m <- matrix(m, nrow = length(pos))
    soma <- if (is.null(soma)) m else soma + m
  }
  list(acum_mat = soma, n_regioes = length(regioes), falhas = character(0),
       dic  = if (!is.null(fit$dic$dic))   fit$dic$dic   else NA_real_,
       waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_)
}

# matrix -> 4 zones (same percentiles as channel)
zonas_de_matriz <- function(acum_mat, semanas = 1:52) {
  qs <- config$inla_percentis %||% c(lower = 0.025, media = 0.50, upper = 0.975)
  data.frame(semana = semanas,
             lower = apply(acum_mat, 1, quantile, qs[["lower"]], na.rm = TRUE),
             media = apply(acum_mat, 1, quantile, qs[["media"]], na.rm = TRUE),
             upper = apply(acum_mat, 1, quantile, qs[["upper"]], na.rm = TRUE))
}

# observed cumulative of the STATE SUMMING the granular level (excludes records missing DRS/CIR/mun!)
obs_acum_estado_por_grao <- function(nivel_grao, alvo_nome, regioes, ano, ate_semana = NULL) {
  col <- switch(nivel_grao, "DRS" = "drs", "CIR" = "cir", "Municipio" = "mun")
  d <- base_menin[!is.na(base_menin[[col]]) & base_menin[[col]] %in% regioes, ]
  cont <- agregar_semanal(d, "Estado", NULL, alvo_nome)
  cont <- cont[cont$ano == ano, c("semana", "casos"), drop = FALSE]
  smax <- if (is.null(ate_semana)) 52L else ate_semana
  if (smax < 1) return(data.frame(semana = integer(0), casos = integer(0), acum = numeric(0)))
  out <- data.frame(semana = 1:smax)
  out <- merge(out, cont, by = "semana", all.x = TRUE)
  out$casos[is.na(out$casos)] <- 0L
  out <- out[order(out$semana), ]
  out$acum <- cumsum(out$casos)
  out
}

# --- ENGINE: generates the reconstructed STATE band (one version: with or without Besag)
reconstruir_estado <- function(nivel_grao, alvo_nome, usar_besag) {
  tag <- if (isTRUE(usar_besag)) "Besag" else "independent"
  message("  [Reconstruction ", tag, "] ", nivel_grao, " -> State | ", alvo_nome)
  regioes <- regioes_do_grao(nivel_grao, usar_besag)
  if (!length(regioes)) { message("    no regions; skipping."); return(invisible(NULL)) }
  
  res <- if (isTRUE(usar_besag))
    soma_grao_besag(nivel_grao, alvo_nome, config$anos_treino, 1:52, config$n_amostras)
  else
    soma_grao_indep(nivel_grao, alvo_nome, regioes, config$anos_treino, 1:52,
                    config$n_amostras, config$reconstrucao$min_casos_regiao %||% 0)
  if (is.null(res)) { message("    could not reconstruct; skipping."); return(invisible(NULL)) }
  
  canal <- zonas_de_matriz(res$acum_mat, 1:52)
  
  # observed for the target year (summing granular level) + comparison with direct total
  sem_hoje <- semana_final_do_ano(config$ano_atual) + 1L
  obs <- obs_acum_estado_por_grao(nivel_grao, alvo_nome, regioes, config$ano_atual,
                                  ate_semana = max(sem_hoje - 1L, 0L))
  tot_direto <- sum(agregar_semanal(base_menin, "Estado", NULL, alvo_nome) %>%
                      filter(ano == config$ano_atual) %>% pull(casos))
  tot_grao   <- if (nrow(obs)) max(obs$acum) else 0
  if (tot_direto != tot_grao)
    message("    [Reconstruction] Total ", config$ano_atual, " summing ", nivel_grao, " = ",
            tot_grao, " vs direct State total = ", tot_direto,
            " (difference = records without ", nivel_grao, ").")
  
  qs <- config$inla_percentis %||% c(lower = 0.025, media = 0.50, upper = 0.975)
  subt <- paste0("Reconstructed State summing ", res$n_regioes, " ", nivel_grao,
                 " (", tag, ", sum sample by sample); zones = percentiles (",
                 paste0(round(100 * qs[c("lower","media","upper")]), collapse = "/"),
                 ") | Training: ", paste(range(config$anos_treino), collapse = "-"))
  plt <- grafico_canal_zonas(canal, if (nrow(obs) > 0) obs else NULL,
                             paste0("Meningitis (cumulative) — ", alvo_nome),
                             subt, paste0("Cumulative observed ", config$ano_atual),
                             lab_centro = "Median")
  pasta <- file.path(config$pasta_saida, "bandas_Estado_reconstruido", nivel_grao)
  salvar_banda(plt, paste0("banda_estado_recon_", tolower(nivel_grao), "_",
                           if (isTRUE(usar_besag)) "besag" else "indep", "_",
                           limpar_nome(alvo_nome)), pasta)
  invisible(canal)
}

# =============================================================================
#  8.7) CONTINGENCY ANALYSIS — evolution of a consolidated year over the band
# =============================================================================
#  Same logic as the endemic channels code (contingency loop), but the
#  BACKGROUND is OUR INLA historical band (4 zones), not the geometric channel.
#  For an already consolidated year (e.g., 2025):
#    - background = INLA historical band trained WITHOUT the consolidated year;
#    - lines = cumulative of the year seen in each PARTIAL base (up to that base's EW)
#              + the CONSOLIDATED cumulative (current base), all layered on top.
# =============================================================================

# infers the EW cutoff from the label ("EW 13" -> 13); NULL/no number -> 52
corte_da_parcial <- function(rotulo, override = NULL) {
  if (!is.null(override)) return(as.integer(override))
  m <- regmatches(rotulo, regexpr("[0-9]+", rotulo))
  if (length(m) && nzchar(m)) as.integer(m) else 52L
}

# loads and prepares a partial line-list base reusing preparar_base()
carregar_parcial <- function(caminho) {
  if (!file.exists(caminho)) { message("    [Conting] file not found: ", caminho); return(NULL) }
  d <- tryCatch(readRDS(caminho), error = function(e) { message("    [Conting] error reading ", caminho, ": ", e$message); NULL })
  if (is.null(d)) return(NULL)
  d <- as.data.frame(d)
  tryCatch(preparar_base(d), error = function(e) { message("    [Conting] error preparing ", caminho, ": ", e$message); NULL })
}

# cumulative 1..ate_semana of a year/target from an already prepared line-list base
acum_ano_de_base <- function(base_prep, nivel, alvo_nome, ano, ate_semana = 52L) {
  cont <- agregar_semanal(base_prep, nivel, NULL, alvo_nome)
  cont <- cont[cont$ano == ano & cont$semana <= ate_semana, c("semana", "casos"), drop = FALSE]
  if (ate_semana < 1) return(data.frame(semana = integer(0), acum = numeric(0)))
  out <- data.frame(semana = 1:ate_semana)
  out <- merge(out, cont, by = "semana", all.x = TRUE)
  out$casos[is.na(out$casos)] <- 0L
  out <- out[order(out$semana), ]
  out$acum <- cumsum(out$casos)
  out[, c("semana", "acum")]
}

# draws the background band + the partial/consolidated trajectories on top
grafico_contingencia <- function(canal, linhas_df, titulo, subtitulo, cores_linhas,
                                 lab_centro = "Median") {
  ymax <- suppressWarnings(max(c(canal$upper, linhas_df$acum), na.rm = TRUE)) * 1.10
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(c(canal$upper, 1), na.rm = TRUE)
  ymax <- max(ymax, 1); canal$teto <- ymax
  
  rot <- c("Below expected", "Expected", "Alert", "Epidemic")
  if (isTRUE(config$usar_carinhas))
    rot <- paste(config$carinhas_canal[c("below", "expected", "alert", "epidemic")], rot)
  
  ordem <- names(cores_linhas)
  linhas_df$serie <- factor(linhas_df$serie, levels = ordem)
  
  g <- ggplot(canal, aes(x = semana)) +
    geom_ribbon(aes(ymin = upper, ymax = teto,  fill = "Epidemic")) +
    geom_ribbon(aes(ymin = media, ymax = upper, fill = "Alert")) +
    geom_ribbon(aes(ymin = lower, ymax = media, fill = "Expected")) +
    geom_ribbon(aes(ymin = 0,     ymax = lower, fill = "Below expected")) +
    geom_line(aes(y = lower), color = "black", linewidth = 0.3) +
    geom_line(aes(y = upper), color = "black", linewidth = 0.3) +
    geom_line(aes(y = media), color = "#1565C0", linewidth = 0.5, linetype = "dotted") +
    geom_line(data = linhas_df, aes(x = semana, y = acum, color = serie, linewidth = serie)) +
    scale_fill_manual(
      values = c("Below expected" = "#0bc82c", "Expected" = "#edea44",
                 "Alert" = "#fe3d3d", "Epidemic" = "#ba0606"),
      breaks = c("Below expected", "Expected", "Alert", "Epidemic"),
      labels = rot, name = NULL) +
    scale_color_manual(values = cores_linhas, breaks = ordem, name = NULL) +
    scale_linewidth_manual(values = setNames(
      ifelse(ordem == "Consolidado", 1.6, 1.1), ordem), guide = "none") +
    scale_x_continuous(breaks = 1:52, labels = function(x) paste0("EW", sprintf("%02d", x)),
                       limits = c(1, 52), expand = expansion(add = 0)) +
    scale_y_continuous(breaks = function(l) pretty(l, n = 6),
                       expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(title = titulo, subtitle = subtitulo,
         x = "Epidemiological week", y = "Cumulative cases in the year") +
    guides(fill  = guide_legend(order = 1, nrow = 1),
           color = guide_legend(order = 2, nrow = 1)) +
    tema_banda()
  g
}

# generates the contingency of ONE target
contingencia_alvo <- function(alvo_nome, nivel, ano, parciais_prep, cortes,
                              incluir_consolidada, cores_parciais) {
  message("  [Contingency] ", nivel, " | ", alvo_nome, " | year ", ano)
  
  # 1) BACKGROUND: INLA historical band trained WITHOUT the consolidated year
  df_todos  <- agregar_semanal(base_menin, nivel, NULL, alvo_nome)
  anos_fundo <- setdiff(config$anos_treino, ano)
  df_treino <- df_todos %>% filter(ano %in% anos_fundo) %>% completar_semanas(anos_fundo)
  canal <- rodar_inla_acumulado(df_treino, config$n_amostras, config$min_casos_treino,
                                config$modelos_inla[[nivel]],
                                seed = instance_seed(nivel, "SP", alvo_nome, ano, "contingencia"))
  if (is.null(canal)) { message("    null background band; skipping."); return(invisible(NULL)) }
  
  # 2) LINES: each partial up to its EW + consolidated (entire year from the current base)
  linhas <- list(); cores <- c()
  for (rot in names(parciais_prep)) {
    bp <- parciais_prep[[rot]]; if (is.null(bp)) next
    ate <- cortes[[rot]]
    a <- acum_ano_de_base(bp, nivel, alvo_nome, ano, ate)
    if (nrow(a)) { a$serie <- rot; linhas[[rot]] <- a; cores[rot] <- cores_parciais[[rot]] %||% "#4292C6" }
  }
  if (isTRUE(incluir_consolidada)) {
    a <- acum_ano_de_base(base_menin, nivel, alvo_nome, ano, 52L)
    if (nrow(a)) { a$serie <- "Consolidado"; linhas[["Consolidado"]] <- a; cores["Consolidado"] <- "#000000" }
  }
  if (!length(linhas)) { message("    no trajectories; skipping."); return(invisible(NULL)) }
  linhas_df <- do.call(rbind, linhas)
  
  qs <- config$inla_percentis %||% c(lower = 0.025, media = 0.50, upper = 0.975)
  subt <- paste0(" ")
  titulo <- paste0("Contingency ", ano, "\nMeningitis (cumulative) — ", alvo_nome)
  
  plt <- grafico_contingencia(canal, linhas_df, titulo, subt, cores)
  pasta <- file.path(config$pasta_saida, "Contingencia", as.character(ano))
  salvar_banda(plt, paste0("contingencia_", tolower(nivel), "_", ano, "_", limpar_nome(alvo_nome)), pasta)
}

executar_contingencia <- function() {
  cfg <- config$contingencia
  message("\n===== CONTINGENCY ANALYSIS — year ", cfg$ano, " =====\n")
  
  # loads the partial bases ONCE (reused by all targets)
  parciais_prep <- list(); cortes <- list()
  for (rot in names(cfg$parciais)) {
    caminho <- file.path(cfg$pasta_parciais %||% ".", cfg$parciais[[rot]])
    message("  [Contingency] loading partial ", rot, ": ", caminho)
    parciais_prep[[rot]] <- carregar_parcial(caminho)
    cortes[[rot]] <- corte_da_parcial(rot, cfg$corte_por_parcial[[rot]])
  }
  parciais_prep <- parciais_prep[!vapply(parciais_prep, is.null, logical(1))]
  if (!length(parciais_prep) && !isTRUE(cfg$incluir_consolidada)) {
    message("  [Contingency] no valid partial base and consolidated is off; aborting.\n"); return(invisible(NULL))
  }
  
  alvos_c <- if (identical(cfg$alvos, "Todos")) names(config$alvos) else cfg$alvos
  nivel   <- cfg$nivel %||% "Estado"
  for (alvo_nome in alvos_c)
    tryCatch(contingencia_alvo(alvo_nome, nivel, cfg$ano, parciais_prep, cortes,
                               cfg$incluir_consolidada, cfg$cores_parciais),
             error = function(e) message("    [ERROR contingency / ", alvo_nome, "] ", e$message))
  message("\n=== Contingency completed (folder 'Contingencia/", cfg$ano, "'). ===\n")
}

# =============================================================================
message("\n============= HISTORICAL BAND (CUMULATIVE) — MENINGITIS =============\n")
alvos_rodar <- names(config$alvos)
usar_besag  <- isTRUE(config$rodar_besag) && identical(config$metodo_banda, "inla")
if (isTRUE(config$rodar_besag) && !identical(config$metodo_banda, "inla"))
  message("(Warning: rodar_besag only applies to metodo_banda='inla'; in geometric channel ",
          "each region is calculated independently.)\n")

regioes_por_nivel <- setNames(lapply(config$niveis, function(nivel) {
  switch(nivel,
         "Estado"    = "SP",
         "DRS"       = if (!is.null(config$drs_selecionadas)) config$drs_selecionadas else listar_regioes("DRS"),
         "CIR"       = if (!is.null(config$cir_selecionadas)) config$cir_selecionadas else listar_regioes("CIR"),
         "Municipio" = if (!is.null(config$mun_selecionados)) config$mun_selecionados else listar_regioes("Municipio"))
}), config$niveis)

n_fig <- sum(vapply(regioes_por_nivel, length, integer(1))) * length(alvos_rodar)
message("Will attempt to generate ~", n_fig, " figures.\n")

for (nivel in config$niveis) {
  message("\n>>> Level: ", nivel, "\n")
  if (usar_besag && nivel != "Estado") {
    tryCatch(processar_nivel_besag(nivel, alvos_rodar),
             error = function(e) message("    [BESAG ERROR ", nivel, "] ", e$message))
  } else {
    for (regiao in regioes_por_nivel[[nivel]])
      for (alvo_nome in alvos_rodar)
        tryCatch(processar_regiao_alvo(nivel, regiao, alvo_nome),
                 error = function(e) message("    [ERROR] ", nivel, " / ", regiao, " / ", alvo_nome, ": ", e$message))
  }
}
message("\n=== Bands completed! Figures in: ", normalizePath(config$pasta_saida, mustWork = FALSE), " ===\n")

# =============================================================================
#  9.5) EXECUTION — STATE RECONSTRUCTED FROM GRANULAR LEVEL
# =============================================================================
if (isTRUE(config$reconstrucao$ativar)) {
  message("\n===== STATE RECONSTRUCTED FROM ",
          config$reconstrucao$nivel_grao, " =====\n")
  alvos_rec <- if (identical(config$reconstrucao$alvos, "Todos")) names(config$alvos)
  else config$reconstrucao$alvos
  for (alvo_nome in alvos_rec)
    for (ub in config$reconstrucao$besag)
      tryCatch(reconstruir_estado(config$reconstrucao$nivel_grao, alvo_nome, ub),
               error = function(e)
                 message("    [ERROR reconstruction ", if (isTRUE(ub)) "besag" else "indep",
                         " / ", alvo_nome, "] ", e$message))
}

# =============================================================================
#  9.6) EXECUTION — CONTINGENCY ANALYSIS (switch)
# =============================================================================
if (isTRUE(config$contingencia$ativar)) executar_contingencia()

# =============================================================================
#  10) VALIDATION — CUMULATIVE (LOYO and blocks), via INLA model (Option A)
#      It's a predictive benchmark, INDEPENDENT of the band drawing method.
# =============================================================================
serie_validacao <- function(contagem, anos, ano_atual) {
  partes <- list()
  for (yr in anos) {
    if (yr == ano_atual) {
      w <- contagem$semana[contagem$ano == yr]; if (length(w) == 0) next; smax <- max(w)
    } else smax <- 52L
    g <- data.frame(ano = yr, semana = 1:smax)
    g <- merge(g, contagem[contagem$ano == yr, c("semana", "casos")], by = "semana", all.x = TRUE)
    g$casos[is.na(g$casos)] <- 0L
    partes[[as.character(yr)]] <- g[order(g$semana), c("ano", "semana", "casos")]
  }
  if (length(partes) == 0) return(NULL)
  do.call(rbind, partes)
}
ajustar_fold_acum <- function(treino, teste, modelo, n_amostras, seed = NULL) {
  anos_tr <- sort(unique(treino$ano)); treino$ano_idx <- match(treino$ano, anos_tr)
  idx_novo <- length(anos_tr) + 1L; n_test <- nrow(teste)
  teste2 <- data.frame(casos = NA_integer_, semana = teste$semana, ano_idx = idx_novo)
  dados  <- rbind(treino[, c("casos", "semana", "ano_idx")], teste2)
  if (is.null(modelo)) modelo <- list()
  familia <- modelo$familia %||% "poisson"; formula_inla <- construir_formula(modelo); nh <- nh_do_modelo(modelo)
  if (nh > 0) dados <- adicionar_harmonicos(dados, nh)
  fit <- tryCatch(inla(formula_inla, family = familia, data = dados,
                       control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                       control.predictor = list(link = 1, compute = TRUE), verbose = FALSE),
                  error = function(e) { message("      [ERROR fold] ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  n_tr <- nrow(treino)
  if (!is.null(seed)) set.seed(seed)
  amostras <- inla.posterior.sample(n_amostras, fit)
  pred_rows <- grep("^Predictor", rownames(amostras[[1]]$latent)); linhas <- pred_rows[(n_tr + 1):nrow(dados)]
  idx_zip <- idx_zip_hyper(amostras, familia); acum_mat <- simular_acumulado(amostras, linhas, idx_zip)
  if (is.null(dim(acum_mat))) acum_mat <- matrix(acum_mat, nrow = n_test)
  list(acum_mat = acum_mat,
       dic  = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
       waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_)
}
metricas_acum <- function(fold, teste, nivel_cob) {
  acum_mat <- fold$acum_mat; y_acum <- cumsum(teste$casos)
  n <- length(y_acum); pred_mean <- rowMeans(acum_mat)
  mae <- mean(abs(y_acum - pred_mean)); rmse <- sqrt(mean((y_acum - pred_mean)^2))
  a <- (1 - nivel_cob) / 2
  qlow  <- apply(acum_mat, 1, quantile, probs = a, names = FALSE)
  qhigh <- apply(acum_mat, 1, quantile, probs = 1 - a, names = FALSE)
  cobertura <- mean(y_acum >= qlow & y_acum <= qhigh) * 100
  alphas <- 1 - c(0.50, 0.80, 0.90, 0.95); K <- length(alphas)
  wis_i <- vapply(seq_len(n), function(i) {
    am <- acum_mat[i, ]; med <- median(am); soma <- 0.5 * abs(y_acum[i] - med)
    for (al in alphas) {
      l <- quantile(am, al / 2, names = FALSE); u <- quantile(am, 1 - al / 2, names = FALSE)
      is_a <- (u - l) + (2 / al) * (l - y_acum[i]) * (y_acum[i] < l) + (2 / al) * (y_acum[i] - u) * (y_acum[i] > u)
      soma <- soma + (al / 2) * is_a
    }
    soma / (K + 0.5)
  }, numeric(1))
  wis <- mean(wis_i, na.rm = TRUE); escala <- mean(y_acum, na.rm = TRUE)
  wis_rel <- if (is.finite(escala) && escala > 0) wis / escala else NA_real_
  erro_total <- abs(pred_mean[n] - y_acum[n])
  mape_final <- if (y_acum[n] > 0) abs(pred_mean[n] - y_acum[n]) / y_acum[n] * 100 else NA_real_
  metricas <- data.frame(MAE = mae, RMSE = rmse, Cobertura = cobertura, WIS = wis, WIS_rel = wis_rel,
                         ErroTotalFinal = erro_total, MAPE_final = mape_final, DIC = fold$dic, WAIC = fold$waic)
  pred <- data.frame(ano = teste$ano, semana = teste$semana, obs_acum = y_acum,
                     pred_acum = pred_mean, li = qlow, ls = qhigh)
  list(metricas = metricas, pred = pred)
}
validar_instancia_acum <- function(nivel, regiao, alvo_nome, modelo, anos_val,
                                   n_amostras, min_casos, nivel_cob, forward = FALSE, min_treino = 3) {
  contagem <- agregar_semanal(base_menin, nivel, regiao, alvo_nome)
  series <- serie_validacao(contagem, anos_val, config$ano_atual)
  if (is.null(series)) return(NULL)
  anos_disp <- sort(unique(series$ano))
  if (forward) {
    if (length(anos_disp) <= min_treino) return(NULL)
    ks <- (min_treino + 1):length(anos_disp)
    folds_def <- lapply(ks, function(k) list(yr = anos_disp[k], tr = anos_disp[seq_len(k - 1)]))
  } else folds_def <- lapply(anos_disp, function(yr) list(yr = yr, tr = setdiff(anos_disp, yr)))
  folds_m <- list(); preds <- list()
  for (fd in folds_def) {
    yr <- fd$yr; treino <- series[series$ano %in% fd$tr, ]; teste <- series[series$ano == yr, ]
    if (nrow(teste) == 0 || sum(treino$casos) < min_casos || length(unique(treino$ano)) < 3) next
    fold <- ajustar_fold_acum(treino, teste, modelo, n_amostras,
                              seed = instance_seed(nivel, regiao, alvo_nome, yr, "fold_loyo"))
    if (is.null(fold)) next
    mp <- metricas_acum(fold, teste, nivel_cob); mp$metricas$ano <- yr
    folds_m[[as.character(yr)]] <- mp$metricas; preds[[as.character(yr)]] <- mp$pred
  }
  if (length(folds_m) == 0) return(NULL)
  folds_df <- do.call(rbind, folds_m)
  medias <- sapply(c("MAE","RMSE","Cobertura","WIS","WIS_rel","ErroTotalFinal","MAPE_final","DIC","WAIC"),
                   function(cc) mean(folds_df[[cc]], na.rm = TRUE))
  resumo <- data.frame(nivel = nivel, regiao = regiao, alvo = alvo_nome, n_folds = nrow(folds_df),
                       t(medias), check.names = FALSE)
  folds_df$nivel <- nivel; folds_df$regiao <- regiao; folds_df$alvo <- alvo_nome
  pred_df <- do.call(rbind, preds); pred_df$nivel <- nivel; pred_df$regiao <- regiao; pred_df$alvo <- alvo_nome
  list(resumo = resumo, folds = folds_df, pred = pred_df)
}
fig_validacao_acum <- function(pred_df, titulo, subtitulo, arquivo, pasta_val) {
  p <- ggplot(pred_df, aes(x = semana)) +
    geom_ribbon(aes(ymin = li, ymax = ls), fill = "#BDBDBD", alpha = 0.6) +
    geom_line(aes(y = pred_acum), color = "#1565C0", linewidth = 0.6) +
    geom_point(aes(y = obs_acum), color = "#212121", size = 1) +
    facet_wrap(~ ano, scales = "free_y") +
    labs(title = titulo, subtitle = subtitulo, x = "Epidemiological week", y = "Cumulative cases") + tema_banda(10)
  ggsave(file.path(pasta_val, arquivo), p, width = config$largura / config$dpi,
         height = config$altura / config$dpi, units = "in", dpi = config$dpi, device = "png", bg = "white")
}
fig_resumo_metricas <- function(resumo_df, titulo, arquivo, pasta_val) {
  mets <- intersect(c("MAE","RMSE","Cobertura","WIS","WIS_rel","ErroTotalFinal","MAPE_final","DIC","WAIC"), names(resumo_df))
  if (!length(mets) || is.null(resumo_df) || !nrow(resumo_df)) return(invisible(NULL))
  d <- resumo_df; d$inst <- paste0(d$nivel, ":", d$regiao, "\n", d$alvo)
  long <- tidyr::pivot_longer(d, cols = dplyr::all_of(mets), names_to = "metrica", values_to = "valor")
  p <- ggplot(long, aes(x = inst, y = valor, fill = inst)) +
    geom_col(show.legend = FALSE) + facet_wrap(~ metrica, scales = "free_y") +
    labs(title = titulo, subtitle = "Mean across folds (lower is better, except Coverage %)", x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(color = "grey40", hjust = 0.5), strip.text = element_text(face = "bold"))
  ggsave(file.path(pasta_val, arquivo), p, width = config$largura / config$dpi,
         height = config$altura / config$dpi + 1, units = "in", dpi = config$dpi, device = "png", bg = "white")
}
executar_validacao_independente <- function(forward, prefixo, titulo_extra) {
  message("\n===== CUMULATIVE VALIDATION", titulo_extra, " (INLA model) =====\n")
  pasta_val <- file.path(config$pasta_saida, "Validacao")
  if (!dir.exists(pasta_val)) dir.create(pasta_val, recursive = TRUE)
  anos_val <- sort(setdiff(unique(c(config$anos_treino, config$ano_atual)), c(2020, 2021, 2022)))
  min_treino <- config$validacao$block_min_treino %||% 3
  instancias <- list()
  if (isTRUE(config$validacao$incluir_estado)) instancias <- c(instancias, list(list(nivel = "Estado", regiao = "SP")))
  for (r in config$validacao$drs)        instancias <- c(instancias, list(list(nivel = "DRS",        regiao = r)))
  for (r in config$validacao$cir)        instancias <- c(instancias, list(list(nivel = "CIR",        regiao = r)))
  for (r in config$validacao$municipios) instancias <- c(instancias, list(list(nivel = "Municipio", regiao = r)))
  alvos_val <- if (identical(config$validacao$alvos, "Todos")) names(config$alvos) else config$validacao$alvos
  resumos <- list(); folds_all <- list()
  for (inst in instancias) for (alvo_nome in alvos_val) {
    rot <- paste0(inst$nivel, ":", inst$regiao, " | ", alvo_nome); message(">>> Validating ", rot)
    res <- tryCatch(validar_instancia_acum(inst$nivel, inst$regiao, alvo_nome,
                                           config$modelos_inla[[inst$nivel]], anos_val, config$validacao$n_amostras,
                                           config$validacao$min_casos_treino, config$validacao$nivel_cobertura,
                                           forward = forward, min_treino = min_treino),
                    error = function(e) { message("    [VALIDATION ERROR] ", e$message); NULL })
    if (is.null(res)) next
    resumos[[rot]] <- res$resumo; folds_all[[rot]] <- res$folds
    rot_reg <- if (inst$nivel == "Estado") "SP_State" else limpar_nome(inst$regiao); rot_alv <- limpar_nome(alvo_nome)
    fig_validacao_acum(res$pred, paste0("Validation (cumulative", titulo_extra, ") — ", alvo_nome),
                       "Blue line = predicted; band = interval; points = observed",
                       paste0(prefixo, tolower(inst$nivel), "_", rot_reg, "_", rot_alv, ".png"), pasta_val)
    message("    -> figure saved.")
  }
  if (length(resumos)) {
    resumo_df <- do.call(rbind, resumos)
    write.csv(resumo_df, file.path(pasta_val, paste0(prefixo, "metricas.csv")), row.names = FALSE, fileEncoding = "UTF-8")
    fig_resumo_metricas(resumo_df, paste0("Summary of metrics (cumulative", titulo_extra, ")"),
                        paste0(prefixo, "resumo_metricas.png"), pasta_val)
    if (length(folds_all))
      write.csv(do.call(rbind, folds_all), file.path(pasta_val, paste0(prefixo, "metricas_por_fold.csv")),
                row.names = FALSE, fileEncoding = "UTF-8")
    message("\n=== Validation saved in: ", normalizePath(pasta_val, mustWork = FALSE), " ===\n")
  } else message("\n[Validation] No instance produced results.\n")
}
if (isTRUE(config$validacao$ativar))
  executar_validacao_independente(forward = FALSE, prefixo = "validacao_loyo_", titulo_extra = ", LOYO")
if (isTRUE(config$validacao$block_ativar))
  executar_validacao_independente(forward = TRUE, prefixo = "validacao_block_", titulo_extra = ", blocks/forward")

# =============================================================================
#  11) SPATIAL VALIDATION (Besag) ON CUMULATIVE — LOYO and/or blocks
# =============================================================================
ajustar_fold_besag_acum <- function(treino, teste, G, modelo, n_amostras, seed = NULL) {
  anos_tr <- sort(unique(treino$ano)); treino$ano_idx <- match(treino$ano, anos_tr); idx_novo <- length(anos_tr) + 1L
  teste_fit <- data.frame(casos = NA_integer_, semana = teste$semana, ano_idx = idx_novo, regiao_idx = teste$regiao_idx)
  dados <- rbind(treino[, c("casos", "semana", "ano_idx", "regiao_idx")], teste_fit)
  familia <- modelo$familia %||% "poisson"; nh <- nh_do_modelo(modelo)
  if (nh > 0) dados <- adicionar_harmonicos(dados, nh)
  formula_besag <- construir_formula_besag(modelo, G)
  fit <- tryCatch(inla(formula_besag, family = familia, data = dados,
                       control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                       control.predictor = list(link = 1, compute = TRUE), verbose = FALSE),
                  error = function(e) { message("      [ERROR besag fold] ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  n_tr <- nrow(treino)
  if (!is.null(seed)) set.seed(seed)
  amostras <- inla.posterior.sample(n_amostras, fit)
  pred_rows <- grep("^Predictor", rownames(amostras[[1]]$latent)); idx_zip <- idx_zip_hyper(amostras, familia)
  list(amostras = amostras, pred_rows = pred_rows, idx_zip = idx_zip, n_tr = n_tr,
       dic = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
       waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_)
}
validar_nivel_besag_acum <- function(nivel, alvo_nome, anos_val, n_amostras, min_casos, nivel_cob, forward, min_treino) {
  gb <- preparar_grafo_besag(nivel); regioes <- gb$regioes; G <- gb$grafo
  if (length(regioes) < 2) { message("  [Besag val] < 2 regions; skipping."); return(NULL) }
  modelo <- modelo_besag_efetivo(nivel)
  series_lst <- lapply(regioes, function(r) {
    cont <- agregar_semanal(base_menin, nivel, r, alvo_nome)
    s <- serie_validacao(cont, anos_val, config$ano_atual); if (!is.null(s)) s$regiao <- r; s
  })
  series <- do.call(rbind, series_lst)
  if (is.null(series) || nrow(series) == 0) return(NULL)
  series$regiao_idx <- match(series$regiao, regioes); anos_disp <- sort(unique(series$ano))
  if (forward) {
    if (length(anos_disp) <= min_treino) return(NULL)
    ks <- (min_treino + 1):length(anos_disp); folds_def <- lapply(ks, function(k) list(yr = anos_disp[k], tr = anos_disp[seq_len(k - 1)]))
  } else folds_def <- lapply(anos_disp, function(yr) list(yr = yr, tr = setdiff(anos_disp, yr)))
  folds_m <- list(); preds <- list()
  for (fd in folds_def) {
    yr <- fd$yr; treino <- series[series$ano %in% fd$tr, ]; teste <- series[series$ano == yr, ]
    if (nrow(teste) == 0 || sum(treino$casos) < min_casos || length(unique(treino$ano)) < 3) next
    message("    [Besag val] ", nivel, " | ", alvo_nome, " | left out year: ", yr, " (", length(regioes), " regions)...")
    fj <- ajustar_fold_besag_acum(treino, teste, G, modelo, n_amostras,
                                  seed = instance_seed(nivel, "ALL", alvo_nome, yr, "fold_besag"))
    if (is.null(fj)) next
    for (r in regioes) {
      idx_r <- which(teste$regiao == r); if (length(idx_r) == 0) next
      teste_r <- teste[idx_r, c("ano", "semana", "casos")]; teste_r <- teste_r[order(teste_r$semana), ]
      linhas <- fj$pred_rows[fj$n_tr + idx_r[order(teste$semana[idx_r])]]
      acum_mat <- simular_acumulado(fj$amostras, linhas, fj$idx_zip)
      if (is.null(dim(acum_mat))) acum_mat <- matrix(acum_mat, nrow = length(idx_r))
      fold_r <- list(acum_mat = acum_mat, dic = fj$dic, waic = fj$waic)
      mp <- metricas_acum(fold_r, teste_r, nivel_cob); mp$metricas$ano <- yr; mp$metricas$regiao <- r; mp$pred$regiao <- r
      folds_m[[paste0(r, "__", yr)]] <- mp$metricas; preds[[paste0(r, "__", yr)]] <- mp$pred
    }
  }
  if (length(folds_m) == 0) return(NULL)
  folds_df <- do.call(rbind, folds_m)
  resumo <- do.call(rbind, lapply(split(folds_df, folds_df$regiao), function(d) {
    medias <- sapply(c("MAE","RMSE","Cobertura","WIS","WIS_rel","ErroTotalFinal","MAPE_final","DIC","WAIC"),
                     function(cc) mean(d[[cc]], na.rm = TRUE))
    data.frame(nivel = nivel, regiao = d$regiao[1], alvo = alvo_nome, n_folds = nrow(d), t(medias), check.names = FALSE)
  }))
  folds_df$nivel <- nivel; folds_df$alvo <- alvo_nome
  pred_df <- do.call(rbind, preds); pred_df$nivel <- nivel; pred_df$alvo <- alvo_nome
  list(resumo = resumo, folds = folds_df, pred = pred_df)
}
executar_validacao_besag <- function(forward, prefixo, titulo_extra) {
  message("\n===== SPATIAL VALIDATION (Besag) ON CUMULATIVE", titulo_extra, " =====\n")
  pasta_val <- file.path(config$pasta_saida, "Validacao")
  if (!dir.exists(pasta_val)) dir.create(pasta_val, recursive = TRUE)
  anos_val <- sort(setdiff(unique(c(config$anos_treino, config$ano_atual)), c(2020, 2021, 2022)))
  min_treino <- config$validacao$block_min_treino %||% 3
  alvos_val <- if (identical(config$validacao$alvos, "Todos")) names(config$alvos) else config$validacao$alvos
  resumos_b <- list(); folds_b <- list()
  for (nivel in config$validacao$besag_niveis) {
    if (nivel == "Estado") { message("  [Besag val] 'Estado' without neighborhood; skipping."); next }
    for (alvo_nome in alvos_val) {
      message(">>> [Besag val] ", nivel, " | ", alvo_nome)
      res <- tryCatch(validar_nivel_besag_acum(nivel, alvo_nome, anos_val, config$validacao$n_amostras,
                                               config$validacao$min_casos_treino, config$validacao$nivel_cobertura,
                                               forward = forward, min_treino = min_treino),
                      error = function(e) { message("    [ERROR Besag val] ", e$message); NULL })
      if (is.null(res)) next
      chave <- paste0(nivel, "|", alvo_nome); resumos_b[[chave]] <- res$resumo; folds_b[[chave]] <- res$folds
      filtro <- config$validacao$besag_regioes[[nivel]]; todas_reg <- unique(res$pred$regiao)
      regioes_rel <- if (!is.null(filtro) && length(filtro))
        todas_reg[padronizar_mun(todas_reg) %in% padronizar_mun(filtro)] else todas_reg
      if ((is.null(filtro) || !length(filtro)) && length(regioes_rel) > 40) {
        message("    [Besag val] ", length(regioes_rel), " regions: skipping figures (see CSV)."); regioes_rel <- character(0)
      }
      for (r in regioes_rel) {
        pr <- res$pred[res$pred$regiao == r, ]; if (nrow(pr) == 0) next
        fig_validacao_acum(pr, paste0("Besag Validation (cumulative", titulo_extra, ") ", ": ", r, " — ", alvo_nome),
                           "Joint spatial fit | line=predicted, band=interval, points=observed",
                           paste0(prefixo, tolower(nivel), "_", limpar_nome(r), "_", limpar_nome(alvo_nome), ".png"), pasta_val)
      }
      message("    -> ", nivel, " / ", alvo_nome, " completed.")
    }
  }
  if (length(resumos_b)) {
    resumo_besag <- do.call(rbind, resumos_b)
    write.csv(resumo_besag, file.path(pasta_val, paste0(prefixo, "metricas.csv")), row.names = FALSE, fileEncoding = "UTF-8")
    fig_resumo_metricas(resumo_besag, paste0("Summary of Besag metrics (cumulative", titulo_extra, ")"),
                        paste0(prefixo, "resumo_metricas.png"), pasta_val)
    if (length(folds_b))
      write.csv(do.call(rbind, folds_b), file.path(pasta_val, paste0(prefixo, "metricas_por_fold.csv")),
                row.names = FALSE, fileEncoding = "UTF-8")
    message("\n=== Besag validation saved in: ", normalizePath(pasta_val, mustWork = FALSE), " ===\n")
  } else message("\n[Besag Validation] No results.\n")
}
if (isTRUE(config$validacao$besag_ativar))
  executar_validacao_besag(forward = FALSE, prefixo = "validacao_besag_loyo_", titulo_extra = ", LOYO")
if (isTRUE(config$validacao$block_besag_ativar))
  executar_validacao_besag(forward = TRUE, prefixo = "validacao_besag_block_", titulo_extra = ", blocks/forward")

# =============================================================================
#  11.5) BLOCK VALIDATION (forward-chaining) OF STATE RECONSTRUCTED BY LEVEL
# =============================================================================
#  For each fold: trains ONLY with PREVIOUS years, predicts the test year in
#  all regions of the granular level, SUMS sample by sample -> predicted STATE, and compares
#  with observed STATE SUMMING the same level (to be apples-to-apples; records
#  without DRS/CIR/municipality are left out on BOTH sides).
#  Runs in both requested versions: without Besag (independent) and with Besag.
# =============================================================================
validar_estado_reconstruido_block <- function(nivel_grao, alvo_nome, usar_besag) {
  tag <- if (isTRUE(usar_besag)) "besag" else "indep"
  regioes <- regioes_do_grao(nivel_grao, usar_besag)
  if (!length(regioes)) return(NULL)
  
  # only COMPLETE years (52 weeks) — current year is incomplete and would distort
  anos_disp  <- sort(unique(config$anos_treino))
  min_treino <- config$reconstrucao$val_min_treino %||% 3
  if (length(anos_disp) <= min_treino) {
    message("    [Val recon] insufficient years for forward-chaining."); return(NULL)
  }
  n_am <- config$reconstrucao$val_n_amostras      %||% 500
  cob  <- config$reconstrucao$val_nivel_cobertura %||% 0.90
  minr <- config$reconstrucao$min_casos_regiao    %||% 0
  
  folds_m <- list(); preds <- list()
  for (k in (min_treino + 1):length(anos_disp)) {
    yr      <- anos_disp[k]
    anos_tr <- anos_disp[seq_len(k - 1)]
    message("    [Val recon ", tag, "] ", nivel_grao, " | ", alvo_nome,
            " | test year: ", yr, " (training: ", paste(range(anos_tr), collapse = "-"), ")")
    res <- tryCatch(
      if (isTRUE(usar_besag))
        soma_grao_besag(nivel_grao, alvo_nome, anos_tr, 1:52, n_am)
      else
        soma_grao_indep(nivel_grao, alvo_nome, regioes, anos_tr, 1:52, n_am, minr),
      error = function(e) { message("      [ERROR recon fold] ", e$message); NULL })
    if (is.null(res)) next
    
    obs   <- obs_acum_estado_por_grao(nivel_grao, alvo_nome, regioes, yr, 52)
    teste <- data.frame(ano = yr, semana = obs$semana, casos = obs$casos)
    mp <- metricas_acum(list(acum_mat = res$acum_mat, dic = res$dic, waic = res$waic), teste, cob)
    mp$metricas$ano <- yr
    folds_m[[as.character(yr)]] <- mp$metricas
    preds[[as.character(yr)]]   <- mp$pred
  }
  if (!length(folds_m)) return(NULL)
  
  folds_df <- do.call(rbind, folds_m)
  mets <- c("MAE","RMSE","Cobertura","WIS","WIS_rel","ErroTotalFinal","MAPE_final","DIC","WAIC")
  medias <- sapply(mets, function(cc) mean(folds_df[[cc]], na.rm = TRUE))
  rot_nivel <- paste0("State<-", nivel_grao, " (", tag, ")")
  resumo <- data.frame(nivel = rot_nivel, regiao = "SP", alvo = alvo_nome,
                       n_folds = nrow(folds_df), t(medias), check.names = FALSE)
  folds_df$nivel <- rot_nivel; folds_df$regiao <- "SP"; folds_df$alvo <- alvo_nome
  pred_df <- do.call(rbind, preds)
  pred_df$nivel <- rot_nivel; pred_df$regiao <- "SP"; pred_df$alvo <- alvo_nome
  list(resumo = resumo, folds = folds_df, pred = pred_df)
}

executar_validacao_reconstrucao <- function() {
  message("\n===== BLOCK VALIDATION — STATE RECONSTRUCTED BY ",
          config$reconstrucao$nivel_grao, " =====\n")
  pasta_val <- file.path(config$pasta_saida, "Validacao")
  if (!dir.exists(pasta_val)) dir.create(pasta_val, recursive = TRUE)
  nivel_grao <- config$reconstrucao$nivel_grao
  alvos_rec  <- if (identical(config$reconstrucao$alvos, "Todos")) names(config$alvos)
  else config$reconstrucao$alvos
  
  resumos <- list(); folds_all <- list()
  for (alvo_nome in alvos_rec) for (ub in config$reconstrucao$besag) {
    tag <- if (isTRUE(ub)) "besag" else "indep"
    res <- tryCatch(validar_estado_reconstruido_block(nivel_grao, alvo_nome, ub),
                    error = function(e) { message("    [ERROR val recon] ", e$message); NULL })
    if (is.null(res)) next
    chave <- paste0(tag, "|", alvo_nome)
    resumos[[chave]] <- res$resumo; folds_all[[chave]] <- res$folds
    fig_validacao_acum(
      res$pred,
      paste0("Block validation — State reconstructed by ", nivel_grao, " (", tag, ") — ", alvo_nome),
      "Training only with previous years | line=predicted, band=interval, points=observed",
      paste0("validacao_block_estado_recon_", tolower(nivel_grao), "_", tag, "_",
             limpar_nome(alvo_nome), ".png"), pasta_val)
    message("    -> figure saved (", tag, " / ", alvo_nome, ").")
  }
  
  if (length(resumos)) {
    resumo_df <- do.call(rbind, resumos)
    write.csv(resumo_df, file.path(pasta_val, "validacao_block_estado_recon_metricas.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    fig_resumo_metricas(resumo_df,
                        paste0("Summary — State reconstructed by ", nivel_grao, " (block)"),
                        "validacao_block_estado_recon_resumo_metricas.png", pasta_val)
    write.csv(do.call(rbind, folds_all),
              file.path(pasta_val, "validacao_block_estado_recon_metricas_por_fold.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    message("\n=== Reconstruction validation saved in: ",
            normalizePath(pasta_val, mustWork = FALSE), " ===\n")
  } else message("\n[Val recon] No results.\n")
}

if (isTRUE(config$reconstrucao$ativar) && isTRUE(config$reconstrucao$validacao_block))
  executar_validacao_reconstrucao()

# =============================================================================
#  12) CONSOLIDATED TXT TABLE
# =============================================================================
{
  pasta_val <- file.path(config$pasta_saida, "Validacao")
  fontes <- list("LOYO" = "validacao_loyo_metricas.csv", "Blocks" = "validacao_block_metricas.csv",
                 "LOYO+Besag" = "validacao_besag_loyo_metricas.csv", "Blocks+Besag" = "validacao_besag_block_metricas.csv",
                 "Blocks+StateRecon" = "validacao_block_estado_recon_metricas.csv")
  partes <- list()
  for (nm in names(fontes)) {
    arq <- file.path(pasta_val, fontes[[nm]])
    if (file.exists(arq)) {
      df <- tryCatch(read.csv(arq, stringsAsFactors = FALSE, fileEncoding = "UTF-8"), error = function(e) NULL)
      if (!is.null(df) && nrow(df)) partes[[nm]] <- cbind(validacao = nm, df)
    }
  }
  if (length(partes)) {
    todas_cols <- unique(unlist(lapply(partes, names)))
    partes <- lapply(partes, function(d) { for (c0 in setdiff(todas_cols, names(d))) d[[c0]] <- NA; d[, todas_cols, drop = FALSE] })
    tab <- do.call(rbind, partes); num <- vapply(tab, is.numeric, logical(1)); tab[num] <- lapply(tab[num], round, 4)
    arq_txt <- file.path(pasta_val, "validation_ALL_metrics.txt")
    con <- file(arq_txt, open = "w", encoding = "UTF-8")
    writeLines(c("==========================================================================",
                 " HISTORICAL MENINGITIS BAND (CUMULATIVE) — VALIDATION METRICS",
                 paste0(" Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
                 " (Validation uses the INLA model as a predictive benchmark.)",
                 "    MAE/RMSE, Coverage, WIS/WIS_rel, TotalFinalError, MAPE_final, DIC/WAIC",
                 "==========================================================================", ""), con)
    utils::write.table(format(tab, justify = "left", na.encode = TRUE), con, sep = "  ", quote = FALSE, row.names = FALSE)
    close(con)
    message("\n=== Consolidated TXT table: ", normalizePath(arq_txt, mustWork = FALSE), " ===\n")
  } else message("\n[TXT Table] No metrics CSV to consolidate.\n")
}

# =============================================================================
#  13) COMBINED FIGURE — side-by-side panels sharing a single legend
# =============================================================================
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))

gerar_painel_alvo <- function(nivel, regiao, alvo_nome, letra) {
  df_todos  <- agregar_semanal(base_menin, nivel, regiao, alvo_nome)
  df_treino <- df_todos %>% filter(ano %in% config$anos_treino) %>% completar_semanas(config$anos_treino)
  
  sem_hoje <- semana_final_do_ano(config$ano_atual) + 1L
  obs_atual <- df_todos %>% filter(ano == config$ano_atual, semana < sem_hoje) %>% select(semana, casos)
  if (sem_hoje > 1) {
    obs_atual <- data.frame(semana = 1:(sem_hoje - 1)) %>%
      left_join(obs_atual, by = "semana") %>%
      mutate(casos = tidyr::replace_na(casos, 0L)) %>% arrange(semana) %>% mutate(acum = cumsum(casos))
  } else obs_atual <- data.frame(semana = integer(0), casos = integer(0), acum = numeric(0))
  
  metodo <- config$metodo_banda %||% "canal_geometrico"
  if (identical(metodo, "canal_geometrico")) {
    canal <- banda_canal_geometrico(df_treino, config$anos_treino,
                                    conf = config$canal$nivel_confianca %||% 0.95,
                                    tipo = config$canal$tipo_intervalo %||% "ic_media",
                                    offset = config$canal$offset_log %||% 1,
                                    min_casos = config$min_casos_treino)
  } else {
    canal <- rodar_inla_acumulado(df_treino, config$n_amostras, config$min_casos_treino,
                                  config$modelos_inla[[nivel]],
                                  seed = instance_seed(nivel, regiao, alvo_nome, config$ano_atual, "banda"))
  }
  if (is.null(canal)) return(NULL)
  
  lab_obs <- paste0("Cumulative observed ", config$ano_atual)
  plt <- grafico_canal_zonas(canal, if (nrow(obs_atual) > 0) obs_atual else NULL,
                             paste0(letra, ") ", alvo_nome), " ", lab_obs, lab_centro = "Median")
  plt + theme(plot.subtitle = element_blank())
}

gerar_figura_combinada <- function(nivel = "Estado", regiao = "SP",
                                   alvos_combo = c("Bacterial Meningitis", "Viral Meningitis",
                                                   "Pneumococcus", "Meningococcus",
                                                   "Meningococcus B", "Meningococcus C")) {
  letras  <- letters[seq_along(alvos_combo)]
  paineis <- Map(function(al, lt) gerar_painel_alvo(nivel, regiao, al, lt), alvos_combo, letras)
  paineis <- Filter(Negate(is.null), paineis)
  if (!length(paineis)) { message("  [Combined figure] no panel was generated."); return(invisible(NULL)) }
  
  combinada <- patchwork::wrap_plots(paineis, ncol = 2) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  pasta <- file.path(config$pasta_saida, "Figura_Combinada")
  if (!dir.exists(pasta)) dir.create(pasta, recursive = TRUE)
  arq <- file.path(pasta, paste0("figura_combinada_", tolower(nivel), "_", regiao, ".png"))
  ggsave(arq, combinada, width = 2400/150, height = 1800/150, units = "in", dpi = 150, bg = "white")
  message("  -> Combined figure saved to: ", arq)
  invisible(combinada)
}

gerar_figura_combinada()
