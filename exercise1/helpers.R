# Helper functions shared by group1.qmd, group2.qmd and group3.qmd (loaded in the pages, not shown to students)
LOD <- 0.3
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
UNIT <- "ng/ml"   # unit of the measurement (the tutorial page sets "pg/mg")
ylab_of <- function(scale) if (scale == "log2") "log2 insulin (ng/ml)" else "insulin (ng/ml)"
show <- function(x) print(x, row.names = FALSE)
spl <- function(x) { x <- as.character(unlist(x)); if (length(x) == 0 || identical(x, "")) character(0) else unlist(strsplit(x, ",")) }

GRP <- list("all mice pooled" = character(0), "diet" = "diet", "plate (diets pooled)" = "plate",
            "plate and diet" = c("plate", "diet"), "cage (diets pooled)" = "cage", "cage and diet" = c("cage", "diet"))
make_g <- function(d, vars) {
  if (length(vars) == 0) factor(rep("all mice", nrow(d)))
  else { l <- d[vars]; for (v in intersect(vars, c("plate", "cage"))) l[[v]] <- paste(v, l[[v]])   # labels name the variable: "plate 1 / chow"
         interaction(l, drop = TRUE, sep = " / ", lex.order = TRUE) }
}

counts <- function(d, by) {   # one row per mouse
  cat("Number of mice\n")
  print(if (by == "diet x cage") table(diet = d$diet, cage = d$cage)
  else if (by == "diet x plate") table(diet = d$diet, plate = d$plate)
  else if (by == "cage x plate") table(cage = d$cage, plate = d$plate)
  else if (by == "diet x cage x plate") table(diet = d$diet, cage = d$cage, plate = d$plate)
  else table(d[[by]], dnn = by))
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
  out <- data.frame(group = levels(g), round(r[, which, drop = FALSE], 2), check.names = FALSE)
  if (scale == "log2") names(out)[names(out) == "geometric mean"] <- paste0("geometric mean (", UNIT, ")")   # back-transformed, not on the log2 scale
  names(out)[names(out) %in% c("n", "missing")] <- paste(names(out)[names(out) %in% c("n", "missing")], "(mice)")
  out
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

# ---- one item: chosen plot and chosen numbers
item <- function(d, plot, stats, by, scale, label = FALSE) {
  stats <- spl(stats)
  if (plot == "(none)" && length(stats) == 0) { cat("Select a plot and/or numbers.\n"); return(invisible()) }
  if (plot != "(none)") draw_plot(d, plot, scale, label)
  if (length(stats)) { cat(if (scale == "log2") "Numbers on the log2 scale (geometric mean in original units)\n" else "Numbers on the original scale\n"); show(stats_tab(d, by, scale, stats)) }
}

EFFECTS <- c("difference of means", "ratio of means", "ratio of medians", "log2 fold change", "2^(log2 fold change)")
effect <- function(d, which) {
  d <- d[!is.na(d$value), ]
  m <- tapply(d$value, d$diet, mean); md <- tapply(d$value, d$diet, median); l <- tapply(log2(d$value), d$diet, mean)
  a <- c("difference of means" = m[["HFD"]] - m[["chow"]], "ratio of means" = m[["HFD"]] / m[["chow"]],
         "ratio of medians" = md[["HFD"]] / md[["chow"]], "log2 fold change" = l[["HFD"]] - l[["chow"]],
         "2^(log2 fold change)" = 2^(l[["HFD"]] - l[["chow"]]))
  which <- intersect(EFFECTS, spl(which))
  data.frame(measure = which, value = round(unname(a[which]), 3))
}

# ---- Group 1: missing values
# numeric column `value`; handling of "<LOD" and "NA" entries
LOD_H <- c("drop", "replace by LOD/2", "replace by LOD", "replace by 0")
prep <- function(d, lod = "drop", na = "drop") {
  v <- d$insulin_ng_ml
  d$value <- suppressWarnings(as.numeric(v))
  isl <- v == "<LOD"; isn <- v == "NA"
  gm <- tapply(d$value, d$diet, mean, na.rm = TRUE)
  d$value[isl] <- switch(lod, "drop" = NA, "replace by LOD/2" = LOD / 2, "replace by LOD" = LOD, "replace by 0" = 0)
  d$value[isn] <- switch(na, "drop" = NA, "replace by group mean" = gm[as.character(d$diet[isn])])
  d
}
missing_view <- function(d, lab, which) {
  nn <- d[!grepl("^[0-9.]+$", d$insulin_ng_ml), ]; nn$coding <- nn$insulin_ng_ml
  if (which == "All raw values") show(d[, c("mouse_id", "cage", "diet", "plate", "insulin_ng_ml")])
  else if (which == "Non-numeric entries only") show(nn[, c("mouse_id", "diet", "plate", "coding")])
  else if (which == "Non-numeric entries with lab record") show(merge(nn[, c("mouse_id", "diet", "plate", "coding")], lab, by = "mouse_id"))
  else if (which == "Count of non-numeric entries by diet") print(table(diet = nn$diet, coding = nn$coding))
  else print(table(plate = nn$plate, diet = nn$diet, coding = nn$coding))
}
MECH <- c(MCAR = "MCAR: missingness is unrelated to the insulin value and to the other variables. Dropping the values reduces n and leaves the group means unbiased.",
          MAR = "MAR: missingness depends on observed variables (diet, plate, cage) but not on the insulin value itself. Dropping can bias estimates unless these variables are accounted for.",
          MNAR = "MNAR: missingness depends on the unobserved value itself. Dropping biases the estimates; the missingness mechanism has to be modelled.")
ml_censored <- function(v, k) {   # log-normal, k values censored below LOD
  y <- log(v)
  nll <- function(p) -(sum(dnorm(y, p[1], exp(p[2]), log = TRUE)) + k * pnorm(log(LOD), p[1], exp(p[2]), log.p = TRUE))
  o <- optim(c(mean(y), log(sd(y))), nll); c(mu = o$par[1], sigma = exp(o$par[2]))
}
compare_handling <- function(d, na) {
  r <- do.call(rbind, lapply(LOD_H, function(h) {
    dd <- prep(d, h, na); m <- tapply(dd$value, dd$diet, mean, na.rm = TRUE)
    data.frame(handling = h, mean_chow = m[["chow"]], mean_HFD = m[["HFD"]], log2FC = log2fc(dd))
  }))
  dd <- prep(d, "drop", "drop")
  p <- sapply(c("chow", "HFD"), function(l) {
    k <- sum(d$insulin_ng_ml == "<LOD" & d$diet == l); v <- dd$value[dd$diet == l & !is.na(dd$value)]
    e <- ml_censored(v, k); c(exp(e[1] + e[2]^2 / 2), e[1] / log(2))
  })
  r <- rbind(r, data.frame(handling = "maximum likelihood, censored log-normal (NA dropped)", mean_chow = p[1, 1], mean_HFD = p[1, 2],
                           log2FC = p[2, 2] - p[2, 1]))
  r[, -1] <- round(r[, -1], 3); r
}
log2fc <- function(d) {
  y <- log2(d$value); ok <- !is.na(y)
  m <- tapply(y[ok], d$diet[ok], mean)
  unname(m["HFD"] - m["chow"])
}

# ---- Group 2: outlier rules, evaluated within each diet group
RULES <- c("boxplot rule", "z-score", "robust z-score")
rule_stat <- function(v, rule) {
  if (rule == "boxplot rule") {   # distance beyond the nearer quartile in IQR units; flagged if |s| > 1.5
    q <- quantile(v, c(.25, .75)); i <- q[2] - q[1]
    s <- ifelse(v > median(v), (v - q[2]) / i, (v - q[1]) / i); f <- abs(s) > 1.5
  } else if (rule == "z-score") { s <- (v - mean(v)) / sd(v); f <- abs(s) > 3
  } else { s <- 0.6745 * (v - median(v)) / mad(v, constant = 1); f <- abs(s) > 3.5 }
  list(s = as.numeric(s), f = as.logical(f))
}
rule_table <- function(d, rules, scales) {
  d <- d[!is.na(d$value), ]; d <- d[order(d$mouse_id), ]
  out <- d[, c("mouse_id", "diet", "plate", "insulin_ng_ml")]
  for (sc in intersect(c("original", "log2"), spl(scales))) for (ru in intersect(RULES, spl(rules))) {
    cell <- character(nrow(d))
    for (l in levels(d$diet)) {
      j <- d$diet == l; r <- rule_stat(on_scale(d$value[j], sc), ru)
      cell[j] <- sprintf("%.2f%s", r$s, ifelse(r$f, " *", ""))
    }
    out[[paste(sub(" rule| score", "", ru), sc)]] <- cell
  }
  out
}
flagged_any <- function(d) {
  d <- d[!is.na(d$value), ]; ids <- character(0)
  for (ru in RULES) for (sc in c("original", "log2")) for (l in levels(d$diet)) {
    j <- d$diet == l; r <- rule_stat(on_scale(d$value[j], sc), ru); ids <- union(ids, d$mouse_id[j][r$f])
  }
  sort(ids)
}
clean_data <- function(d, ids, vals) {
  d$value <- suppressWarnings(as.numeric(d$insulin_ng_ml))
  ids <- as.character(unlist(ids)); vals <- as.character(unlist(vals))
  for (k in seq_along(ids)) {
    j <- d$mouse_id == ids[k]
    if (vals[k] == "exclude") d$value[j] <- NA
    if (vals[k] == "correct by factor 0.1") d$value[j] <- d$value[j] * 0.1
  }
  d
}

# ---- Group 3: skewness and scale (added; groups 1 and 2 and the tutorial do not use these)
SPREADS <- c("SD", "IQR", "variance")
spread_tab <- function(d, measures, scales) {   # spread per diet and ratio HFD/chow, per scale
  d <- d[!is.na(d$value), ]
  f <- list(SD = sd, IQR = IQR, variance = var)
  do.call(rbind, lapply(intersect(c("original", "log2"), spl(scales)), function(sc)
    do.call(rbind, lapply(intersect(SPREADS, spl(measures)), function(m) {
      v <- sapply(levels(d$diet), function(l) f[[m]](on_scale(d$value[d$diet == l], sc)))
      data.frame(scale = sc, measure = m, chow = round(v[["chow"]], 2), HFD = round(v[["HFD"]], 2),
                 "ratio HFD/chow" = round(v[["HFD"]] / v[["chow"]], 2), check.names = FALSE)
    }))))
}
