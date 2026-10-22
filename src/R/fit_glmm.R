# Fit one binomial mixed model with glmmTMB and write its outputs (Phase 3).
#
# Usage: Rscript fit_glmm.R spec.json
# spec fields:
#   input   CSV with y, w, train (1/0), row_id, the grouping columns and the fixed-effect columns
#   outdir  folder for the outputs
#   fixed   right-hand side of the fixed-effect formula (columns of input)
#   random  grouping columns, each given a random intercept
#   grid    optional CSV of covariate settings: fixed-effect predictions plus n_draws draws
#   n_draws, seed
# Outputs (outdir): test_eta.csv (row_id, eta; fixed effects only), coef.csv, vcov.csv,
#   grid_eta.csv and grid_draws.csv (if grid), and fit.json last, which marks the fit complete.
# Called by plate_truth.rfit; all data preparation and metrics are done in Python.

suppressPackageStartupMessages({
  library(glmmTMB); library(splines); library(jsonlite); library(MASS)
})

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
  rhs <- reformulate(attr(terms(lme4::nobars(form)), "term.labels"))
  if (nrow(test) > 0) {
    X <- model.matrix(rhs, test)[, names(b), drop = FALSE]
    write.csv(data.frame(row_id = test$row_id, eta = as.vector(X %*% b)),
              file.path(spec$outdir, "test_eta.csv"), row.names = FALSE)
  }
  if (!is.null(spec$grid)) {
    grid <- read.csv(spec$grid)
    X <- model.matrix(rhs, grid)[, names(b), drop = FALSE]
    write.csv(data.frame(grid, eta = as.vector(X %*% b)),
              file.path(spec$outdir, "grid_eta.csv"), row.names = FALSE)
    set.seed(spec$seed)
    B <- if (info$pdHess) mvrnorm(spec$n_draws, b, V) else matrix(b, 1)
    write.csv(X %*% t(B), file.path(spec$outdir, "grid_draws.csv"), row.names = FALSE)
  }
}
write(toJSON(info, auto_unbox = TRUE, digits = NA, null = "null"),
      file.path(spec$outdir, "fit.json"))
