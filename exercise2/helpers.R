# Helper functions shared by group1.qmd and group2.qmd (loaded in the pages, not shown to students)
COLS <- c(chow = "#8A939C", HFD = "#E8851C")
PAL <- c("#566E98", "#E8851C", "#54A24B", "#B279A2", "#9D755D", "#C9A227")

load_data <- function(file) {
  d <- read.csv(file, colClasses = "character", na.strings = character(0))
  d$diet <- factor(d$diet, c("chow", "HFD"))
  d
}

# numeric column `value`; handling of "<LOD" and "NA" entries
prep <- function(d, lod = "drop", na = "drop") {
  v <- d$insulin_ng_ml
  d$value <- suppressWarnings(as.numeric(v))
  isl <- v == "<LOD"; isn <- v == "NA"
  gm <- tapply(d$value, d$diet, mean, na.rm = TRUE)
  d$value[isl] <- switch(lod, "drop" = NA, "replace by LOD/2" = LOD / 2,
                         "replace by LOD" = LOD, "replace by 0" = 0)
  d$value[isn] <- switch(na, "drop" = NA, "replace by group mean" = gm[as.character(d$diet[isn])])
  d
}

on_scale <- function(x, scale) if (scale == "log2") log2(x) else x
ylab_of <- function(scale) if (scale == "log2") "log2 insulin (ng/ml)" else "insulin (ng/ml)"
show <- function(x, w = 125) print(x, row.names = FALSE, width = w)
spl <- function(x) { x <- as.character(unlist(x)); if (length(x) == 0 || identical(x, "")) character(0) else unlist(strsplit(x, ",")) }

GRP <- list("all mice pooled" = character(0), "diet" = "diet", "plate (diets pooled)" = "plate",
            "plate and diet" = c("plate", "diet"), "cage (diets pooled)" = "cage", "cage and diet" = c("cage", "diet"))
make_g <- function(d, vars) {
  if (length(vars) == 0) factor(rep("all mice", nrow(d)))
  else interaction(d[vars], drop = TRUE, sep = " / ", lex.order = TRUE)
}

counts <- function(d, by) {
  if (by == "diet x cage") table(diet = d$diet, cage = d$cage)
  else if (by == "diet x plate") table(diet = d$diet, plate = d$plate)
  else if (by == "cage x plate") table(cage = d$cage, plate = d$plate)
  else if (by == "diet x cage x plate") table(diet = d$diet, cage = d$cage, plate = d$plate)
  else table(d[[by]], dnn = by)
}

# ---- descriptive numbers per group
STATS <- c("n", "missing", "mean", "geometric mean", "SD", "variance", "median", "IQR", "min", "max", "skewness")
stats_tab <- function(d, by, scale, which) {
  which <- intersect(STATS, spl(which))
  g <- make_g(d, GRP[[by]])
  r <- do.call(rbind, lapply(levels(g), function(l) {
    v <- d$value[g == l]; miss <- sum(is.na(v)); v <- v[!is.na(v)]; y <- on_scale(v, scale)
    c(n = length(y), missing = miss, mean = mean(y), "geometric mean" = exp(mean(log(v))), SD = sd(y),
      variance = var(y), median = median(y), IQR = IQR(y), min = min(y), max = max(y),
      skewness = mean((y - mean(y))^3) / mean((y - mean(y))^2)^1.5)
  }))
  data.frame(group = levels(g), round(r[, which, drop = FALSE], 2), check.names = FALSE)
}

# ---- plots: spec is "<kind>, by <grouping>" or "<kind>, <grouping>"
draw_plot <- function(d, spec, scale = "original", label = FALSE) {
  kind <- sub(",.*$", "", spec); grp <- sub("^[^,]*, (by )?", "", spec)
  d <- d[!is.na(d$value), ]; d$y <- on_scale(d$value, scale)
  bad <- !is.finite(d$y)
  if (any(bad)) { cat("Not shown (log2 of 0):", sum(bad), "value(s)\n"); d <- d[!bad, ] }
  vars <- GRP[[grp]]; g <- make_g(d, vars); K <- nlevels(g); y <- d$y; yl <- ylab_of(scale)
  gc <- if (identical(vars, "diet")) COLS else setNames(PAL[seq_len(K)], levels(g))
  pc <- COLS[as.character(d$diet)]
  jit <- function() { set.seed(1); as.integer(g) + runif(nrow(d), -.15, .15) }
  panels <- function() if (K > 1) par(mfrow = c(ceiling(K / 2), 2))
  op <- par(mar = c(4, 4, 2, 1)); on.exit(par(op))
  if (kind == "Dot plot") {
    xj <- jit()
    plot(xj, y, xaxt = "n", xlim = c(.5, K + .5), xlab = "", ylab = yl, pch = 19, col = pc)
    axis(1, 1:K, levels(g), cex.axis = .9)
    if (label) text(xj, y, d$mouse_id, pos = 4, cex = .7)
    if (K == 2 && identical(vars, "diet")) NULL
    else if (!identical(vars, "diet")) legend("topleft", names(COLS), pch = 19, col = COLS, bty = "n", horiz = TRUE)
  } else if (kind == "Mean and 95% CI") {
    mu <- tapply(y, g, mean); n <- tapply(y, g, length); se <- tapply(y, g, sd) / sqrt(n); h <- qt(.975, n - 1) * se
    xj <- jit()
    plot(xj, y, xaxt = "n", xlim = c(.5, K + .5), xlab = "", ylab = yl, pch = 19, col = adjustcolor(pc, .5))
    axis(1, 1:K, levels(g), cex.axis = .9)
    segments(1:K, mu - h, 1:K, mu + h, lwd = 2); points(1:K, mu, pch = 18, cex = 2)
  } else if (kind == "Boxplot") {
    boxplot(y ~ g, ylab = yl, xlab = "", outline = !label, cex.axis = .9)
    if (label) { xj <- jit(); points(xj, y, pch = 19, col = pc); text(xj, y, d$mouse_id, pos = 4, cex = .7) }
  } else if (kind == "Histogram") {
    panels()
    for (l in levels(g)) hist(y[g == l], xlim = range(y), xlab = yl, main = l, col = gc[l])
  } else if (kind == "ECDF") {
    plot(NA, xlim = range(y), ylim = c(0, 1), xlab = yl, ylab = "empirical distribution function")
    for (l in levels(g)) lines(ecdf(y[g == l]), col = gc[l], verticals = TRUE)
    legend("bottomright", levels(g), col = gc, lty = 1, bty = "n")
  } else {  # QQ plot
    panels()
    for (l in levels(g)) { qqnorm(y[g == l], main = l, pch = 19, col = gc[l]); qqline(y[g == l]) }
  }
}


# ======================================================= Exercise 2
spl <- function(x) { x <- as.character(unlist(x)); if (length(x) == 0 || identical(x, "")) character(0) else unlist(strsplit(x, ",")) }

load2 <- function(g) {
  d <- read.csv(sprintf("ex2_group%d.csv", g), colClasses = c(plate = "character", cage = "character"))
  d$diet <- factor(d$diet, c("chow", "HFD"))
  m <- unique(d[if (g == 1) d$week == 12 else TRUE, c("mouse_id", "cage", "diet", "start_weight_g", "final_weight_g")])
  names(m)[4:5] <- c("sw", "fw"); w <- NULL
  if (g == 1) {
    e <- d[d$week == 12, ]; b <- d[d$week == 0, ]
    m$plate <- e$plate[match(m$mouse_id, e$mouse_id)]
    m$end <- e$insulin_ng_ml[match(m$mouse_id, e$mouse_id)]; m$base <- b$insulin_ng_ml[match(m$mouse_id, b$mouse_id)]
  } else {
    m$plate <- d$plate[match(m$mouse_id, d$mouse_id)]
    m$end <- tapply(d$insulin_ng_ml, d$mouse_id, mean)[m$mouse_id]
    m$well1 <- d$insulin_ng_ml[d$well == 1][match(m$mouse_id, d$mouse_id[d$well == 1])]
    m$well2 <- d$insulin_ng_ml[d$well == 2][match(m$mouse_id, d$mouse_id[d$well == 2])]
    w <- d; names(w)[names(w) == "start_weight_g"] <- "sw"; names(w)[names(w) == "final_weight_g"] <- "fw"; w$value <- w$insulin_ng_ml
  }
  m <- m[order(m$mouse_id), ]; m$value <- m$end; rownames(m) <- NULL
  list(d = d, m = m, w = w)
}

counts <- function(d, by) {   # rows and distinct mice per group
  v <- strsplit(by, " x ")[[1]]
  cat("Number of rows\n"); print(table(d[v], useNA = "ifany"))
  mi <- d[!duplicated(d[c("mouse_id", v)]), ]
  na <- !complete.cases(mi[v]); mi <- mi[!na, ]
  cat("\nNumber of mice\n"); print(table(mi[v]))
  if (any(na)) cat("\n", sum(na), " mouse-by-group combinations with NA in the grouping variable (week 0 rows, no plate recorded) are not counted as mice.\n", sep = "")
}

# ---- fitted-model helpers
coef_vec <- function(fit, term = "dietHFD") {
  s <- summary(fit)$coefficients[term, ]; ci <- confint(fit)[term, ]
  c(estimate = s[[1]], SE = s[[2]], t = s[[3]], df = fit$df.residual, lo = ci[[1]], hi = ci[[2]], p = s[[4]])
}
fmt_p <- function(p) unname(ifelse(p < 1e-4, "< 0.0001", formatC(p, format = "f", digits = 4)))
STATS2 <- c("estimate", "SE", "df", "95% CI", "p-value", "fold change (2^estimate)")
fit_row <- function(label, fit, unit, stats) {
  v <- coef_vec(fit); log2u <- unit == "log2"
  all <- c(estimate = sprintf("%.3f", v["estimate"]), SE = sprintf("%.3f", v["SE"]), df = sprintf("%.0f", v["df"]),
           "95% CI" = sprintf("%.3f to %.3f", v["lo"], v["hi"]), "p-value" = fmt_p(v["p"]),
           "fold change (2^estimate)" = if (log2u) sprintf("%.2f (%.2f to %.2f)", 2^v["estimate"], 2^v["lo"], 2^v["hi"]) else "n/a")
  data.frame(analysis = label, unit = unit, "n (rows)" = nobs(fit), t(all[intersect(STATS2, stats)]), check.names = FALSE)
}
rows_df <- function(l) { r <- do.call(rbind, l); rownames(r) <- NULL; r }

# Group 1: paired design
AN1 <- c("unpaired: end values", "paired: end minus baseline (log2 ratio)", "paired: end minus baseline (difference in ng/ml)",
         "baseline as covariate: log2 end ~ diet + log2 baseline")
fit1 <- function(m, a) switch(a,
  "unpaired: end values" = lm(log2(end) ~ diet, m),
  "paired: end minus baseline (log2 ratio)" = lm(I(log2(end) - log2(base)) ~ diet, m),
  "paired: end minus baseline (difference in ng/ml)" = lm(I(end - base) ~ diet, m),
  "baseline as covariate: log2 end ~ diet + log2 baseline" = lm(log2(end) ~ diet + I(log2(base)), m))
unit1 <- function(a) ifelse(a == "paired: end minus baseline (difference in ng/ml)", "ng/ml", "log2")
analyses1 <- function(m, which, stats) {
  which <- intersect(AN1, spl(which)); stats <- spl(stats)
  rows_df(lapply(which, function(a) fit_row(a, fit1(m, a), unit1(a), stats)))
}
# Group 2: nested design
AN2 <- c("all wells as independent observations", "mouse means", "one well per mouse (well 1)", "one well per mouse (well 2)")
fit2 <- function(m, w, a) switch(a,
  "all wells as independent observations" = lm(log2(value) ~ diet, w),
  "mouse means" = lm(log2(end) ~ diet, m),
  "one well per mouse (well 1)" = lm(log2(well1) ~ diet, m),
  "one well per mouse (well 2)" = lm(log2(well2) ~ diet, m))
analyses2 <- function(m, w, which, stats) {
  which <- intersect(AN2, spl(which)); stats <- spl(stats)
  rows_df(lapply(which, function(a) fit_row(a, fit2(m, w, a), "log2", stats)))
}
est_plot <- function(fits, labels, units, as) {   # estimates with 95% CI; fits: list of lm
  v <- sapply(fits, coef_vec); fc <- as == "fold change"
  est <- v["estimate", ]; lo <- v["lo", ]; hi <- v["hi", ]
  keep <- if (fc) units == "log2" else rep(TRUE, length(est))
  if (fc) { est <- 2^est; lo <- 2^lo; hi <- 2^hi }
  est <- est[keep]; lo <- lo[keep]; hi <- hi[keep]; labels <- labels[keep]; K <- length(est)
  if (K == 0) return(invisible())
  op <- par(mar = c(4, 17, 1, 1)); on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, if (fc) 1 else 0)), ylim = c(.5, K + .5), yaxt = "n", ylab = "",
       xlab = if (fc) "fold change HFD / chow" else "difference HFD - chow (log2, or ng/ml where stated)")
  abline(v = if (fc) 1 else 0, lty = 3)
  segments(lo, K:1, hi, K:1, lwd = 2, col = COLS["HFD"]); points(est, K:1, pch = 19, cex = 1.4)
  axis(2, K:1, labels, las = 1, cex.axis = .72)
}
tab_est <- function(fits, labels, units, fc = FALSE) {
  v <- sapply(fits, coef_vec)
  if (!fc) return(data.frame(analysis = labels, unit = units, estimate = round(v["estimate", ], 3), CI_low = round(v["lo", ], 3), CI_high = round(v["hi", ], 3), row.names = NULL))
  k <- units == "log2"
  data.frame(analysis = labels[k], fold_change = round(2^v["estimate", k], 3), CI_low = round(2^v["lo", k], 3), CI_high = round(2^v["hi", k], 3), row.names = NULL)
}
sel_fits1 <- function(m, which) { which <- intersect(AN1, spl(which)); list(fits = lapply(which, function(a) fit1(m, a)), labels = which, units = sapply(which, unit1)) }
sel_fits2 <- function(m, w, which) { which <- intersect(AN2, spl(which)); list(fits = lapply(which, function(a) fit2(m, w, a)), labels = which, units = rep("log2", length(which))) }
est_item <- function(s, as, plot) {
  if (length(s$fits) == 0) { cat("Select at least one analysis.\n"); return(invisible()) }
  if (plot != "(none)") est_plot(s$fits, s$labels, s$units, as)
  if (as == "fold change" && all(s$units != "log2")) { cat("No analysis on the log2 scale selected; analyses on the ng/ml scale have no fold change.\n"); return(invisible()) }
  show(tab_est(s$fits, s$labels, s$units, as == "fold change"))
  if (as == "fold change" && any(s$units != "log2")) cat("\nAnalyses on the ng/ml scale have no fold change and are not shown.\n")
}

# Group 1: relation between baseline and end value
CORS <- c("correlation of all mice", "correlation within chow", "correlation within HFD", "SD of baseline", "SD of end value", "SD of the change (end - baseline)")
rel_item <- function(m, plot, nums, scale) {
  nums <- intersect(CORS, spl(nums)); b <- on_scale(m$base, scale); e <- on_scale(m$end, scale); ch <- e - b
  if (plot == "(none)" && length(nums) == 0) { cat("Select a plot and/or numbers.\n"); return(invisible()) }
  pc <- COLS[as.character(m$diet)]
  if (plot == "Scatter, baseline against end value") {
    op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op))
    plot(b, e, pch = 19, col = pc, xlab = paste("baseline,", ylab_of(scale)), ylab = paste("end,", ylab_of(scale)))
    legend("topleft", names(COLS), pch = 19, col = COLS, bty = "n")
  } else if (plot == "Lines, baseline to end per mouse") {
    op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op))
    plot(NA, xlim = c(.8, 2.2), ylim = range(c(b, e)), xaxt = "n", xlab = "", ylab = ylab_of(scale)); axis(1, 1:2, c("baseline", "end"))
    segments(1, b, 2, e, col = pc, lwd = 1.5); points(rep(1, 20), b, pch = 19, col = pc); points(rep(2, 20), e, pch = 19, col = pc)
  } else if (plot == "Dot plot, change (end - baseline), by diet") {
    op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op)); set.seed(1)
    xj <- as.integer(m$diet) + runif(20, -.15, .15)
    plot(xj, ch, xaxt = "n", xlim = c(.5, 2.5), pch = 19, col = pc, xlab = "", ylab = paste("change,", ylab_of(scale))); axis(1, 1:2, levels(m$diet)); abline(h = 0, lty = 3)
  }
  if (length(nums)) {
    v <- c("correlation of all mice" = cor(b, e), "correlation within chow" = cor(b[m$diet == "chow"], e[m$diet == "chow"]),
           "correlation within HFD" = cor(b[m$diet == "HFD"], e[m$diet == "HFD"]), "SD of baseline" = sd(b), "SD of end value" = sd(e),
           "SD of the change (end - baseline)" = sd(ch))
    show(data.frame(measure = nums, value = round(unname(v[nums]), 3)))
  }
}
# Group 2: wells and mice
NUMS2 <- c("n wells per diet", "n mice per diet", "wells per mouse", "SD between wells of a mouse", "SD between mouse means within diet", "correlation between the two wells")
unit_item <- function(m, w, plot, nums, scale) {
  nums <- intersect(NUMS2, spl(nums))
  if (plot == "(none)" && length(nums) == 0) { cat("Select a plot and/or numbers.\n"); return(invisible()) }
  if (plot != "(none)") {
    kind <- sub(",.*$", "", plot); what <- sub("^[^,]*, ", "", plot)
    if (what == "wells joined by mouse") {
      op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op)); y1 <- on_scale(m$well1, scale); y2 <- on_scale(m$well2, scale); x <- 1:20
      plot(NA, xlim = c(.5, 20.5), ylim = range(c(y1, y2)), xaxt = "n", xlab = "mouse", ylab = ylab_of(scale))
      axis(1, x, sub("M", "", m$mouse_id), cex.axis = .7); segments(x, y1, x, y2, col = COLS[as.character(m$diet)], lwd = 2)
      points(x, y1, pch = 19, col = COLS[as.character(m$diet)]); points(x, y2, pch = 1, col = COLS[as.character(m$diet)])
      legend("topleft", names(COLS), pch = 19, col = COLS, bty = "n", horiz = TRUE)
    } else draw_plot(if (what == "all wells by diet") w else m, paste0(kind, ", by diet"), scale)
  }
  if (length(nums)) {
    y1 <- on_scale(m$well1, scale); y2 <- on_scale(m$well2, scale); ym <- on_scale(m$end, scale)
    v <- list("n wells per diet" = table(w$diet), "n mice per diet" = table(m$diet), "wells per mouse" = table(as.vector(table(w$mouse_id)), dnn = "wells per mouse (the count is the number of mice)"),
              "SD between wells of a mouse" = sd(y1 - y2) / sqrt(2),
              "SD between mouse means within diet" = sqrt(mean(tapply(ym, m$diet, var))),
              "correlation between the two wells" = cor(y1, y2))
    for (k in nums) { cat(k, "\n"); print(round(v[[k]], 3)); cat("\n") }
  }
}

# ---- Exercise 3: two-group model
METHODS <- c("linear model (equal variances)", "Welch t-test", "Wilcoxon rank-sum test")
TGS <- c("estimate", "SE", "t statistic", "df", "95% CI", "p-value", "fold change with 95% CI (log2 scale)")
two_group <- function(m, scale, method, stats, plot) {
  stats <- intersect(TGS, spl(stats))
  if (plot == "(none)" && length(stats) == 0) { cat("Select a plot and/or numbers.\n"); return(invisible()) }
  if (plot != "(none)") draw_plot(m, plot, scale)
  if (length(stats) == 0) return(invisible())
  y <- on_scale(m$end, scale); x <- m$diet
  if (method == "Wilcoxon rank-sum test") {
    r <- c("p-value" = fmt_p(suppressWarnings(wilcox.test(y ~ x)$p.value)))
    show(data.frame(method = method, t(r[intersect(stats, names(r))]), check.names = FALSE))
    cat("\nThe rank-sum test gives a p-value only.\n"); return(invisible())
  }
  if (method == "Welch t-test") {
    t <- t.test(y ~ x); est <- unname(diff(t$estimate)); se <- unname(t$stderr); df <- unname(t$parameter); lo <- -t$conf.int[2]; hi <- -t$conf.int[1]
    ts <- -unname(t$statistic); p <- t$p.value
  } else {
    v <- coef_vec(lm(y ~ diet, data.frame(y = y, diet = x))); est <- v[["estimate"]]; se <- v[["SE"]]; df <- v[["df"]]; lo <- v[["lo"]]; hi <- v[["hi"]]; ts <- v[["t"]]; p <- v[["p"]]
  }
  all <- c(estimate = sprintf("%.3f", est), SE = sprintf("%.3f", se), "t statistic" = sprintf("%.2f", ts), df = sprintf("%.1f", df),
           "95% CI" = sprintf("%.3f to %.3f", lo, hi), "p-value" = fmt_p(p),
           "fold change with 95% CI (log2 scale)" = if (scale == "log2") sprintf("%.2f (%.2f to %.2f)", 2^est, 2^lo, 2^hi) else "n/a")
  show(data.frame(method = method, scale = scale, t(all[stats]), check.names = FALSE))
}

# ---- Exercise 4: assumptions
ATESTS <- c("Shapiro-Wilk test of the residuals", "Shapiro-Wilk test of the raw values (diets pooled)", "SD per diet", "skewness of the residuals")
assump <- function(m, scale, plot, tests) {
  tests <- intersect(ATESTS, spl(tests))
  if (plot == "(none)" && length(tests) == 0) { cat("Select a plot and/or tests.\n"); return(invisible()) }
  y <- on_scale(m$end, scale); fit <- lm(y ~ diet, m); r <- resid(fit)
  if (plot == "QQ plot of the residuals") { qqnorm(r, pch = 19, main = ""); qqline(r) }
  else if (plot == "Histogram of the residuals") hist(r, main = "", xlab = "residual")
  else if (plot == "Residuals against fitted values") { plot(fitted(fit), r, pch = 19, xlab = "fitted", ylab = "residual"); abline(h = 0, lty = 3) }
  else if (plot == "QQ plot of the raw values (diets pooled)") { qqnorm(y, pch = 19, main = ""); qqline(y) }
  for (k in tests) {
    if (k == ATESTS[1]) { s <- shapiro.test(r); cat(sprintf("Shapiro-Wilk, residuals: W = %.3f, p = %s\n", s$statistic, fmt_p(s$p.value))) }
    if (k == ATESTS[2]) { s <- shapiro.test(y); cat(sprintf("Shapiro-Wilk, raw values pooled: W = %.3f, p = %s\n", s$statistic, fmt_p(s$p.value))) }
    if (k == ATESTS[3]) { cat("SD per diet: "); cat(paste(names(tapply(y, m$diet, sd)), round(tapply(y, m$diet, sd), 3), collapse = ", "), "\n") }
    if (k == ATESTS[4]) cat(sprintf("Skewness of the residuals: %.3f\n", mean(r^3) / mean(r^2)^1.5))
  }
}

# ---- Exercise 5: model builder
COVS <- c("start weight", "final weight", "plate", "cage")
TERM <- c("start weight" = "sw", "final weight" = "fw", "plate" = "plate", "cage" = "cage", "baseline" = "base_s", "fasting glucose (week 0)" = "gl")
model_df <- function(m, w, obs, scale) {
  if (grepl("wells", obs)) { df <- w; df$y <- on_scale(df$value, scale) } else { df <- m; df$y <- on_scale(df$end, scale) }
  df$plate <- factor(df$plate); df$cage <- factor(df$cage)
  if ("base" %in% names(df)) df$base_s <- on_scale(df$base, scale)
  df
}
fit_spec <- function(m, w, spec) {
  p <- strsplit(spec, ";")[[1]]; covs <- if (nchar(p[1])) strsplit(p[1], "+", fixed = TRUE)[[1]] else character(0)
  df <- model_df(m, w, p[3], p[2])
  list(fit = lm(reformulate(c("diet", TERM[covs]), "y"), df), covs = covs, scale = p[2], obs = p[3])
}
pretty_term <- function(x) { x <- sub("dietHFD", "diet (HFD vs chow)", x); x <- sub("^sw$", "start weight (per g)", x); x <- sub("^fw$", "final weight (per g)", x)
  x <- sub("plate2", "plate 2 vs 1", x); x <- sub("cageB", "cage B vs A", x); x <- sub("^gl$", "fasting glucose (per mg/dl)", x); sub("base_s", "baseline (per unit)", x) }
model_report <- function(m, w, spec, show_what) {
  show_what <- spl(show_what)
  if (length(show_what) == 0) { cat("Select what to display.\n"); return(invisible()) }
  s <- fit_spec(m, w, spec); fit <- s$fit
  cat("Model: ", if (s$scale == "log2") "log2 insulin" else "insulin", " ~ ", paste(c("diet", s$covs), collapse = " + "), "   (", nobs(fit), " rows)\n\n", sep = "")
  if ("coefficient table" %in% show_what) {
    cf <- summary(fit)$coefficients; ci <- confint(fit)
    show(data.frame(term = pretty_term(rownames(cf)), estimate = round(cf[, 1], 3), SE = round(cf[, 2], 3),
                    CI_low = round(ci[, 1], 3), CI_high = round(ci[, 2], 3), p_value = fmt_p(cf[, 4])))
    cat("\n")
  }
  if ("fit statistics" %in% show_what) {
    sm <- summary(fit)
    show(data.frame(R2 = round(sm$r.squared, 3), adj_R2 = round(sm$adj.r.squared, 3), AIC = round(AIC(fit), 1), BIC = round(BIC(fit), 1),
                    residual_SD = round(sm$sigma, 3), "n (rows)" = nobs(fit), check.names = FALSE))
  }
  if ("QQ plot of the residuals" %in% show_what) { qqnorm(resid(fit), pch = 19, main = ""); qqline(resid(fit)) }
  if ("residuals against fitted values" %in% show_what) { plot(fitted(fit), resid(fit), pch = 19, xlab = "fitted", ylab = "residual"); abline(h = 0, lty = 3) }
}
CMPCOLS <- c("covariates", "observations", "scale", "diet estimate", "SE", "95% CI", "p-value", "R2", "adjusted R2", "AIC", "BIC", "n")
compare_models <- function(m, w, hist, cols) {
  cols <- intersect(CMPCOLS, spl(cols)); specs <- if (nchar(hist)) strsplit(hist, "|", fixed = TRUE)[[1]] else character(0)
  if (length(specs) == 0) { cat("No model in the table. Add models in item 5.1.\n"); return(invisible()) }
  if (length(cols) == 0) { cat("Select at least one column.\n"); return(invisible()) }
  rows <- lapply(seq_along(specs), function(i) {
    s <- fit_spec(m, w, specs[i]); v <- coef_vec(s$fit); sm <- summary(s$fit)
    c(model = i, covariates = if (length(s$covs)) paste(s$covs, collapse = " + ") else "diet only", observations = s$obs, scale = s$scale,
      "diet estimate" = sprintf("%.3f", v["estimate"]), SE = sprintf("%.3f", v["SE"]), "95% CI" = sprintf("%.3f to %.3f", v["lo"], v["hi"]),
      "p-value" = fmt_p(v["p"]), R2 = sprintf("%.3f", sm$r.squared), "adjusted R2" = sprintf("%.3f", sm$adj.r.squared),
      AIC = sprintf("%.1f", AIC(s$fit)), BIC = sprintf("%.1f", BIC(s$fit)), n = nobs(s$fit))
  })
  r <- as.data.frame(do.call(rbind, rows), check.names = FALSE); o <- r[, c("model", cols), drop = FALSE]; names(o)[names(o) == "n"] <- "n (rows)"; show(o, 250)
  if (length(unique(r$observations)) > 1 || length(unique(r$scale)) > 1)
    cat("\nAIC and BIC are comparable only between models fitted to the same observations and the same scale.\n")
  if (length(unique(r$scale)) > 1) cat("Diet estimate, SE and CI are in log2 units on the log2 scale and in ng/ml on the original scale", if (!"scale" %in% cols) " (select the column 'scale' to see which is which)", ".\n", sep = "")
}

# ---- Exercise 6: outlook (precomputed in simulate_ex2.R with lme4/lmerTest)
OUTCOLS <- c("estimate", "SE", "df", "CI_low", "CI_high", "p_value", "SD_mouse", "SD_residual")
outlook <- function(g, methods, cols, allowed = OUTCOLS) {
  o <- read.csv(sprintf("outlook_group%d.csv", g), check.names = FALSE)
  # method names may contain commas (group 3): a selection joined by "|" is split there, a single exact name is kept whole
  methods <- intersect(o$method, if (methods %in% o$method) methods else if (grepl("|", methods, fixed = TRUE)) strsplit(methods, "|", fixed = TRUE)[[1]] else spl(methods)); cols <- intersect(allowed, spl(cols))
  if (length(methods) == 0 || length(cols) == 0) { cat("Select at least one method and one column.\n"); return(invisible()) }
  show(o[match(methods, o$method), c("method", cols), drop = FALSE])
}

# ======================================================= Exercise 2, group 3: covariates and model choice
load3 <- function() {
  d <- read.csv("ex2_group3.csv", colClasses = c(plate = "character", cage = "character"))
  d$diet <- factor(d$diet, c("chow", "HFD"))
  m <- d[, c("mouse_id", "cage", "diet", "plate", "start_weight_g", "glucose_week0_mg_dl", "final_weight_g", "insulin_ng_ml")]
  names(m)[5:8] <- c("sw", "gl", "fw", "end"); m$value <- m$end
  list(d = d, m = m, w = NULL)
}
CV <- c("start weight" = "sw", "fasting glucose (week 0)" = "gl", "final weight" = "fw")
CV_UNIT <- c(sw = "g", gl = "mg/dl", fw = "g")
cov_item <- function(m, plot, nums) {
  nums <- intersect(c("mean and SD per diet: start weight", "mean and SD per diet: fasting glucose (week 0)", "mean and SD per diet: final weight",
                      "difference HFD - chow with 95% CI: start weight", "difference HFD - chow with 95% CI: fasting glucose (week 0)",
                      "difference HFD - chow with 95% CI: final weight"), spl(nums))
  if (plot == "(none)" && length(nums) == 0) { cat("Select a plot and/or numbers.\n"); return(invisible()) }
  if (plot != "(none)") {
    v <- CV[[sub("^[^,]*, ", "", plot)]]; kind <- sub(",.*$", "", plot); lab <- names(CV)[CV == v]
    op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op)); y <- m[[v]]; set.seed(1); xj <- as.integer(m$diet) + runif(nrow(m), -.15, .15)
    if (kind == "Dot plot") { plot(xj, y, xaxt = "n", xlim = c(.5, 2.5), xlab = "", ylab = paste0(lab, " (", CV_UNIT[v], ")"), pch = 19, col = COLS[as.character(m$diet)]); axis(1, 1:2, levels(m$diet)) }
    else boxplot(y ~ m$diet, xlab = "", ylab = paste0(lab, " (", CV_UNIT[v], ")"))
  }
  if (length(nums)) {
    rows <- lapply(nums, function(k) {
      v <- CV[[sub("^[^:]*: ", "", k)]]; y <- m[[v]]; g <- m$diet
      if (startsWith(k, "mean and SD")) data.frame(measure = k, chow = sprintf("%.2f (SD %.2f)", mean(y[g == "chow"]), sd(y[g == "chow"])), HFD = sprintf("%.2f (SD %.2f)", mean(y[g == "HFD"]), sd(y[g == "HFD"])), check.names = FALSE)
      else { t <- t.test(y ~ g); data.frame(measure = k, chow = "", HFD = sprintf("%.2f (%.2f to %.2f)", diff(t$estimate), -t$conf.int[2], -t$conf.int[1])) }
    })
    show(rows_df(rows))
    cat("\nMouse counts: ", paste(names(table(m$diet)), table(m$diet), collapse = ", "), "(n mice per diet).\n")
  }
}
RELN <- c("correlation with start weight, all mice", "correlation with fasting glucose (week 0), all mice", "correlation with final weight, all mice",
          "correlation with start weight, within chow", "correlation with fasting glucose (week 0), within chow", "correlation with final weight, within chow",
          "correlation with start weight, within HFD", "correlation with fasting glucose (week 0), within HFD", "correlation with final weight, within HFD",
          "mean insulin per plate", "mean insulin per cage")
rel3_item <- function(m, plot, nums, scale) {
  nums <- intersect(RELN, spl(nums)); y <- on_scale(m$end, scale); yl <- ylab_of(scale)
  if (plot == "(none)" && length(nums) == 0) { cat("Select a plot and/or numbers.\n"); return(invisible()) }
  if (plot != "(none)") {
    op <- par(mar = c(4, 4, 1, 1)); on.exit(par(op))
    if (startsWith(plot, "Scatter")) {
      v <- CV[[sub("^Scatter, insulin against ", "", plot)]]
      plot(m[[v]], y, pch = 19, col = COLS[as.character(m$diet)], xlab = paste0(names(CV)[CV == v], " (", CV_UNIT[v], ")"), ylab = yl)
      legend("topleft", names(COLS), pch = 19, col = COLS, bty = "n")
    } else {
      f <- factor(m[[if (grepl("plate", plot)) "plate" else "cage"]]); set.seed(1); xj <- as.integer(f) + runif(nrow(m), -.15, .15)
      plot(xj, y, xaxt = "n", xlim = c(.5, nlevels(f) + .5), xlab = if (grepl("plate", plot)) "plate" else "cage", ylab = yl, pch = 19, col = COLS[as.character(m$diet)])
      axis(1, seq_len(nlevels(f)), levels(f)); legend("topleft", names(COLS), pch = 19, col = COLS, bty = "n")
    }
  }
  if (length(nums)) {
    val <- sapply(nums, function(k) {
      if (k == "mean insulin per plate") return(paste(paste0("plate ", levels(factor(m$plate))), round(tapply(y, m$plate, mean), 3), collapse = ", "))
      if (k == "mean insulin per cage") return(paste(paste0("cage ", levels(factor(m$cage))), round(tapply(y, m$cage, mean), 3), collapse = ", "))
      v <- CV[[sub(",.*$", "", sub("correlation with ", "", k))]]; i <- if (grepl("within chow", k)) m$diet == "chow" else if (grepl("within HFD", k)) m$diet == "HFD" else rep(TRUE, nrow(m))
      sprintf("%.3f (n mice = %d)", cor(m[[v]][i], y[i]), sum(i))
    })
    show(data.frame(measure = nums, value = unname(val)))
  }
}
# model of item 2.3: coefficient table and the diet effect with and without the chosen covariates
cov_model <- function(m, covs, show_what) {
  show_what <- spl(show_what); covs <- spl(covs)
  if (length(show_what) == 0) { cat("Select what to display.\n"); return(invisible()) }
  spec <- paste(c(paste(covs, collapse = "+"), "log2", "mice (20 rows)"), collapse = ";")
  if ("coefficient table" %in% show_what) model_report(m, NULL, spec, "coefficient table")
  if ("diet effect with and without the chosen covariates" %in% show_what) {
    f0 <- fit_spec(m, NULL, ";log2;mice (20 rows)")$fit; f1 <- fit_spec(m, NULL, spec)$fit
    r <- rbind(coef_vec(f0), coef_vec(f1))
    cat("\n"); show(data.frame(model = c("diet only", paste(c("diet", covs), collapse = " + ")), estimate = round(r[, "estimate"], 3), SE = round(r[, "SE"], 3),
                                CI_low = round(r[, "lo"], 3), CI_high = round(r[, "hi"], 3), "n (rows)" = nobs(f0)[c(1, 1)], check.names = FALSE))
  }
}
