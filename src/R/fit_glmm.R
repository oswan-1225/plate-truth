# Fit one binomial mixed model with glmmTMB and write its outputs (Phase 3).
#
# Usage: Rscript fit_glmm.R spec.json
# spec fields:
#   input   CSV with y, w, train (1/0), row_id, the grouping columns and the fixed-effect columns
#   outdir  folder for the outputs
#   fixed   right-hand side of the fixed-effect formula (columns of input)
#   random  grouping columns, each given a random intercept
#   known   grouping columns whose level is known for a held-out row (e.g. phylum): their
#           estimated effects are added to the held-out predictions and written out
#   grid    optional CSV of covariate settings: its design matrix and n_draws draws of the
#           fixed effects are written, for the curve
#   n_draws, seed
# Outputs (outdir): test_eta.csv (row_id, eta: fixed effects plus known random effects),
#   coef.csv, vcov.csv, ranef_<group>.csv for each known group, grid_X.csv and draws.csv (if
#   grid), and fit.json last, which marks the fit complete.
# Called by plate_truth.rfit; all data preparation and metrics are done in Python.

suppressPackageStartupMessages({
  library(glmmTMB); library(splines); library(jsonlite); library(MASS)
})
nobars <- if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars else lme4::nobars

spec <- fromJSON(commandArgs(trailingOnly = TRUE)[1])
dir.create(spec$outdir, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(spec$input)
for (g in spec$random) d[[g]] <- factor(d[[g]])
train <- d[d$train == 1, ]
test <- d[d$train == 0, ]

form <- as.formula(paste("y ~", spec$fixed, "+",
                         paste0("(1 | ", spec$random, ")", collapse = " + ")))
warns <- character(0)
t0 <- Sys.time()
m <- withCallingHandlers(
  tryCatch(glmmTMB(form, family = binomial, data = train, weights = w),
           error = function(e) e),
  warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

info <- list(n = nrow(train), events = sum(train$y), n_test = nrow(test), seconds = secs,
             warnings = unique(warns))
if (inherits(m, "error")) {
  info$status <- "error"; info$error <- conditionMessage(m)
} else {
  vc <- VarCorr(m)$cond
  info$status <- "ok"
  info$converged <- m$fit$convergence == 0
  info$pdHess <- isTRUE(m$sdr$pdHess)
  info$sd <- as.list(sapply(spec$random, function(g) attr(vc[[g]], "stddev")[[1]]))
  info$logLik <- as.numeric(logLik(m))
  b <- fixef(m)$cond
  b <- b[!is.na(b)]                      # rank-deficient columns are dropped by glmmTMB
  V <- vcov(m)$cond[names(b), names(b)]
  write.csv(data.frame(term = names(b), estimate = b, se = sqrt(diag(V))),
            file.path(spec$outdir, "coef.csv"), row.names = FALSE)
  write.csv(V, file.path(spec$outdir, "vcov.csv"))
  re <- ranef(m)$cond
  for (g in spec$known) {
    write.csv(data.frame(level = rownames(re[[g]]), effect = re[[g]][, 1]),
              file.path(spec$outdir, paste0("ranef_", g, ".csv")), row.names = FALSE)
  }
  rhs <- reformulate(attr(terms(nobars(form)), "term.labels"))
  if (nrow(test) > 0) {
    X <- model.matrix(rhs, test)[, names(b), drop = FALSE]
    eta <- as.vector(X %*% b)
    for (g in spec$known) {      # levels unseen in training get 0 (the mean effect)
      u <- re[[g]][as.character(test[[g]]), 1]
      eta <- eta + ifelse(is.na(u), 0, u)
    }
    write.csv(data.frame(row_id = test$row_id, eta = eta),
              file.path(spec$outdir, "test_eta.csv"), row.names = FALSE)
  }
  if (!is.null(spec$grid)) {
    grid <- read.csv(spec$grid)
    X <- model.matrix(rhs, grid)[, names(b), drop = FALSE]
    write.csv(X, file.path(spec$outdir, "grid_X.csv"), row.names = FALSE)
    set.seed(spec$seed)
    B <- if (info$pdHess) mvrnorm(spec$n_draws, b, V) else matrix(b, 1, dimnames = list(NULL, names(b)))
    write.csv(B, file.path(spec$outdir, "draws.csv"), row.names = FALSE)
  }
}
write(toJSON(info, auto_unbox = TRUE, digits = NA, null = "null"),
      file.path(spec$outdir, "fit.json"))
