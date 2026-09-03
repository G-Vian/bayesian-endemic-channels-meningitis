# Bayesian Epidemic Bands for Meningitis Surveillance — State of São Paulo, Brazil

R code for building **probabilistic epidemic bands** for meningitis
surveillance using Bayesian structural time-series models fitted with
[R-INLA](https://www.r-inla.org/), and for estimating the **counterfactual effect
of vaccination coverage** on the expected disease burden.

Developed at the Centro de Vigilância Epidemiológica (CVE) of the São Paulo State
Health Department, using notification data from **SINAN** (Sistema de Informação
de Agravos de Notificação, the Brazilian National Notifiable Disease System).

---

## What problem this solves

An endemic channel answers a practical surveillance question: *given what past
years looked like, is the current year running low, normal, or unusually high?*

The classic approach (the geometric endemic channel, Bortman 1999) summarises
past years with a log-scale mean and confidence interval. It is simple and widely
used, but it describes the training years rather than genuinely predicting a new
one, and it has no natural way to incorporate covariates.

This repository implements a Bayesian alternative that keeps the familiar
four-zone output while fixing both limitations:

- the band is a genuine **predictive** interval for an unseen year, because the
  target year enters the model with a *new* year index whose effect is drawn from
  its prior rather than estimated from data;
- covariates (here, **vaccination coverage**) can be added, enabling a
  with-vaccine vs without-vaccine comparison;
- uncertainty is propagated correctly when trajectories are summed across
  regions or etiologies, because sums are taken **sample by sample** from the
  posterior rather than on quantiles.

The classic geometric channel is retained in the code as a benchmark, reproducing
the `endemic.channels` package exactly, so the two methods can be compared
side by side.

---

## Repository contents

| File | Purpose |
|---|---|
| `historical_band_meningitis.R` | Builds the four-zone historical band per region and etiology; optional spatial (Besag) model, state reconstruction from sub-regions, contingency analysis, and predictive validation. |
| `vaccine_bands_meningitis.R` | Fits paired models (with vs without vaccine-coverage covariates) and compares the resulting expected burden; includes a diagnostic log for uncertainty. |

---

## The model

Weekly case counts are modelled as Poisson process. The linear predictor contains an
intercept, a **seasonal term** across the 52 epidemiological weeks (a first-order
random walk by default, or Fourier harmonics), and an **i.i.d. year effect**:

```
log(mu[week, year]) = alpha + f_seasonal(week) + u[year]      u[year] ~ iid
```

To predict a target year, a block of rows with a **new** year index is appended,
so `u[target]` is drawn from its prior rather than fitted. `N` draws are taken
from the posterior; weekly counts are simulated from each draw and accumulated
within that draw:

```
C[t, m] = sum over i = 1..t of Y[i, m]
```

This yields a full predictive matrix of **52 weeks × N samples**. The four epidemic band zones
are quantiles of that matrix (2.5% / 50% / 97.5% by default), matching the
interval level of the classic endemic channel.

In the vaccine script, the second model adds standardised coverage covariates
(primary series, and booster where applicable) to that same linear predictor.

### Validation

Predictive performance is assessed with leave-one-year-out and forward-chaining
(block) validation, reporting MAE, RMSE, empirical interval coverage, the
**Weighted Interval Score (WIS)** and its relative version, final cumulative
error, MAPE, and INLA's DIC/WAIC. Results are written to CSV and to a formatted
TXT summary.

---

## Input data

### ⚠️ Curated data in, easily adaptable to raw SINAN

These scripts read an **already-curated line-list** (`.RDS`) — a dataset in which
raw SINAN records have already been cleaned and decoded into readable columns.

**We produced that curated file ourselves, starting from raw SINAN data.** The
pipeline is therefore straightforward to adapt to raw SINAN extracts: you only
need a preprocessing step that decodes the raw SINAN fields into the column names
listed below. Nothing downstream needs to change, because every routine in these
scripts consumes those column names — not the raw SINAN numeric codes.

Expected columns, one row per notified case:

| Column | Meaning |
|---|---|
| `anoepi`, `sepi` | Epidemiological year and week (integers). If absent, derived from `data_sintoma`. |
| `data_sintoma` | Symptom onset date (fallback source for year/week). |
| `classi_fin` | Final case classification; `"Confirmado"` marks confirmed cases. |
| `grupo_etiologico` | Broad etiological group (e.g. `"Asseptica"`, `"Bacteriana"`). |
| `etiologia` | Agent (e.g. `"Meningococo"`, `"Pneumococo"`, `"Hemofilo"`). |
| `sorogrupo` | Meningococcal serogroup (e.g. `"Sorogrupo B"`, `"Sorogrupo C"`). |
| `MUNI_RESID` | Municipality of residence (name). |
| `DRS_RESID` | Health region of residence (name); optional. |

### Vaccine coverage files (vaccine script only)

One CSV per vaccine and dose. First column = municipality (optionally prefixed
with its IBGE code); remaining columns = one per **year**, named with four digits
(`"2015"`, `"2016"`, …). Both `;` and `,` separators and decimal commas are
handled automatically.

---

## Requirements

R (≥ 4.0) plus:

```r
install.packages(c("dplyr", "tidyr", "stringr", "lubridate",
                   "ggplot2", "sf", "stringi", "patchwork"))
# optional, only for the spatial (Besag) model:
install.packages("spdep")
```

R-INLA is not on CRAN and installs from its own repository (the scripts do this
automatically on first run):

```r
install.packages("INLA",
  repos = c(CRAN = "https://cloud.r-project.org",
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)
```

---

## Usage

1. **Set your paths.** Every file path in the configuration block is a
   placeholder reading `"insira aqui seu diretório"`. Replace each one with a
   real path before running.

2. **Choose training and target years** in the configuration block:

   ```r
   anos_treino = c(2013, 2014, 2015, 2016, 2017, 2018, 2023, 2024, 2025)
   ano_atual   = 2019   # the year whose observed trajectory is overlaid
   ```

   Pandemic years are excluded from training in our analysis, since 2020–2022
   case counts do not reflect the endemic pattern. This is an epidemiological
   judgement, not a technical requirement — adjust it for your context.

3. **Source the whole file.** Figures, CSVs and TXT summaries are written to the
   configured output folder.

The configuration blocks are annotated with ✍️ (safe to edit) and 🚫 (leave alone
unless you know exactly what you are changing).

---

## Reproducibility

`inla.posterior.sample()` draws random samples, so the same model run twice
without a fixed seed returns slightly different numbers. Several routines here
re-fit the *same* model (for example, the no-vaccine baseline is computed once
for the figures and again inside the validation module). Both scripts therefore
derive a **deterministic seed from the instance identity** — level, region,
target, year and a role tag — so the same instance always produces the same
posterior sample regardless of where in the script it is computed. Re-running a
script reproduces its figures and metrics exactly.

---

## Known limitations

- **Week 53.** Years with 53 epidemiological weeks have their week-53 cases
  dropped, keeping every year at exactly 52 points so that years can be stacked
  in one matrix. This follows the `endemic.channels` convention, but it means
  annual totals are marginally underestimated in those years.

- **Extrapolation in the vaccine models.** When a vaccine's historical coverage
  is very stable, its training standard deviation is small, and a moderate
  real-world drop in the target year becomes an extreme z-score outside the
  range the model was trained on. The resulting interval can widen dramatically.
  This is a genuine statistical limitation, not a coding error; the diagnostic
  log flags it explicitly so that affected bands can be interpreted with caution.

- **Reporting delay.** Surveillance data are revised for months after
  notification. The contingency analysis exists to visualise this: it overlays
  how a consolidated year *looked* at several points during that year, using
  partial database snapshots.

---

## Citation and contact

If you use or adapt this code, please cite the associated publication (details to
be added on acceptance).

Developed by the technical team of the Divisão de Doenças de Transmissão
Respiratória, Centro de Vigilância Epidemiológica (CVE/CCD/SES-SP).
