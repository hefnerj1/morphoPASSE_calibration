library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(randomForest)
library(grid)
library(patchwork)

dat_clean <- readRDS("data/derived/dat_clean.rds")
dir.create("figures/generated", showWarnings = FALSE, recursive = TRUE)

set.seed(1234)

non_trait_cols <- c(
  "dataset", "region", "sheet",
  "ind_id", "Ind_ID", "ID",
  "Sex",
  "CasePrediction", "postprob",
  "Accuracy", "Kappa", "Sensitivity", "Specificity", "PPV", "NPV", "Balanced",
  "Collection", "Institution", "Dataset", "Notes"
)

coerce_traits_1to5 <- function(df, candidates, min_valid_prop = 0.95) {
  for (cn in candidates) {
    x <- df[[cn]]
    if (is.numeric(x) || is.integer(x)) next
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      suppressWarnings(x_num <- as.numeric(x))
      ok <- !is.na(x_num)
      if (sum(ok) == 0) next
      valid <- x_num[ok] %in% 1:5
      prop_valid <- mean(valid)
      if (prop_valid >= min_valid_prop) df[[cn]] <- x_num
    }
  }
  df
}

get_trait_cols <- function(df) {
  candidates <- setdiff(names(df), non_trait_cols)
  numeric_candidates <- candidates[vapply(df[candidates], \(x) is.numeric(x) || is.integer(x), logical(1))]
  numeric_candidates <- numeric_candidates[vapply(df[numeric_candidates], \(x) !all(is.na(x)), logical(1))]
  numeric_candidates <- numeric_candidates[vapply(df[numeric_candidates], \(x) {
    x <- x[!is.na(x)]
    length(x) > 0 && all(x %in% 1:5)
  }, logical(1))]
  numeric_candidates
}

perm_ceiling_rfm <- function(df, trait_cols, n_perm = 500, seed = 1234) {
  set.seed(seed)
  X <- dplyr::select(df, all_of(trait_cols))
  n <- nrow(X)
  acc <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    y <- factor(sample(rep(c("F", "M"), length.out = n), size = n, replace = FALSE))
    fit <- randomForest(x = X, y = y)
    pred <- fit$predicted
    acc[i] <- mean(pred == y)
  }
  list(
    acc = acc,
    summary = tibble(
      perm_mean = mean(acc),
      perm_p95 = as.numeric(quantile(acc, 0.95)),
      perm_p99 = as.numeric(quantile(acc, 0.99))
    )
  )
}

summarize_dataset_region <- function(df, dataset_name, region_name, n_perm = 500) {
  trait_cols <- get_trait_cols(df)
  conf <- df %>%
    summarise(
      n = n(),
      post_mean = mean(postprob, na.rm = TRUE),
      post_median = median(postprob, na.rm = TRUE)
    )
  perm <- perm_ceiling_rfm(df, trait_cols, n_perm = n_perm, seed = 1234)
  acc <- df %>%
    summarise(accuracy = if (dataset_name == "known") mean(CasePrediction == Sex, na.rm = TRUE) else NA_real_)
  tibble(
    dataset = dataset_name,
    region = region_name,
    n = conf$n,
    n_traits = length(trait_cols),
    post_mean = conf$post_mean,
    post_median = conf$post_median,
    perm_mean = perm$summary$perm_mean,
    perm_p95 = perm$summary$perm_p95,
    perm_p99 = perm$summary$perm_p99,
    accuracy = acc$accuracy
  )
}

trait_signal_test <- function(df_region, trait_cols, n_perm = 500, seed = 1234) {
  sim_long <- df_region %>%
    dplyr::select(all_of(trait_cols)) %>%
    pivot_longer(cols = everything(), names_to = "trait", values_to = "value") %>%
    mutate(value = as.integer(value))

  p_trait_dist <- ggplot(sim_long, aes(x = factor(value))) +
    geom_bar() +
    facet_wrap(~trait, ncol = 4) +
    labs(
      x = "Trait score (1–5)",
      y = "Count",
      title = paste0("Simulated trait distributions (", unique(df_region$region), ")")
    )

  set.seed(seed)
  X <- dplyr::select(df_region, all_of(trait_cols))
  n <- nrow(X)
  acc <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    y <- factor(sample(rep(c("F", "M"), length.out = n), size = n, replace = FALSE))
    fit <- randomForest(x = X, y = y)
    pred <- fit$predicted
    acc[i] <- mean(pred == y)
  }

  perm_summary <- tibble(
    mean_acc = mean(acc),
    p95_acc = as.numeric(quantile(acc, 0.95)),
    p99_acc = as.numeric(quantile(acc, 0.99))
  )

  conf_sim <- df_region %>%
    summarise(
      n = n(),
      mean_postprob = mean(postprob, na.rm = TRUE),
      median_postprob = median(postprob, na.rm = TRUE),
      pct_ge_0_90 = mean(postprob >= 0.90, na.rm = TRUE),
      pct_ge_0_95 = mean(postprob >= 0.95, na.rm = TRUE),
      min_postprob = min(postprob, na.rm = TRUE),
      max_postprob = max(postprob, na.rm = TRUE)
    )

  list(
    trait_cols = trait_cols,
    trait_distribution_plot = p_trait_dist,
    perm_summary = perm_summary,
    conf_sim = conf_sim,
    perm_acc = acc
  )
}

sim <- dat_clean %>%
  dplyr::filter(dataset == "simulated")

candidate_cols <- setdiff(names(sim), non_trait_cols)
sim <- coerce_traits_1to5(sim, candidate_cols, min_valid_prop = 0.95)

trait_cols_by_region <- sim %>%
  group_by(region) %>%
  group_map(~ tibble(region = unique(.x$region), trait_cols = list(get_trait_cols(.x))), .keep = TRUE) %>%
  bind_rows()

print(trait_cols_by_region)

if (all(lengths(trait_cols_by_region$trait_cols) == 0)) {
  cat("\nNo trait columns detected after coercion.\nHere are candidate columns and their classes (first region only):\n")
  df0 <- sim %>% dplyr::filter(region == unique(sim$region)[1])
  cand0 <- setdiff(names(df0), non_trait_cols)
  diag_tbl <- tibble(
    col = cand0,
    class = vapply(df0[cand0], \(x) paste(class(x), collapse = "/"), character(1)),
    example = vapply(df0[cand0], \(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return("all NA")
      paste(head(unique(x), 5), collapse = ", ")
    }, character(1))
  )
  print(diag_tbl)
}

results <- list()

for (i in seq_len(nrow(trait_cols_by_region))) {
  reg <- trait_cols_by_region$region[i]
  tcols <- trait_cols_by_region$trait_cols[[i]]
  df_reg <- sim %>% dplyr::filter(region == reg)
  if (length(tcols) == 0) {
    cat("\n--- Region:", reg, "---\nNo trait columns detected for this region.\n")
    next
  }
  results[[reg]] <- trait_signal_test(df_reg, tcols, n_perm = 500, seed = 1234)
  print(results[[reg]]$trait_distribution_plot)
  cat("\n--- Region:", reg, "---\n")
  cat("Trait columns detected:\n")
  print(tcols)
  print(results[[reg]]$perm_summary)
  print(results[[reg]]$conf_sim)
}

if (length(results) > 0) {
  perm_plot_df <- bind_rows(lapply(names(results), function(reg) {
    tibble(region = reg, acc = results[[reg]]$perm_acc)
  }))
  perm_plot <- ggplot(perm_plot_df, aes(x = acc)) +
    geom_density() +
    facet_wrap(~region, ncol = 1) +
    labs(
      x = "Permutation RFM accuracy",
      y = "Density",
      title = "Permutation test: separability of simulated traits under random sex labels"
    )
  print(perm_plot)
}

df <- dat_clean %>%
  dplyr::filter(dataset == "simulated", region == "cranial") %>%
  dplyr::select(CasePrediction, G, NC, SOL, SOR, MPL, MPR) %>%
  na.omit()

X <- dplyr::select(df, -CasePrediction)
y <- factor(df$CasePrediction)

set.seed(1234)
rf_fit <- randomForest(x = X, y = y)
cv_acc <- mean(rf_fit$predicted == y)

cv_acc

df <- dat_clean %>%
  dplyr::filter(dataset == "simulated", region == "cranial") %>%
  dplyr::select(CasePrediction, postprob, G, NC, SOL, SOR, MPL, MPR) %>%
  na.omit()

X <- dplyr::select(df, G, NC, SOL, SOR, MPL, MPR)
y <- factor(df$CasePrediction)

set.seed(1234)
fit <- randomForest(x = X, y = y)
p_rf <- fit$votes
pred_class <- fit$predicted
pred_idx <- match(pred_class, colnames(p_rf))
rf_conf <- p_rf[cbind(seq_len(nrow(p_rf)), pred_idx)]

cor(df$postprob, rf_conf, method = "spearman")

n_perm <- 500

sim_fig <- sim %>%
  mutate(region = factor(region, levels = c("cranial", "pelvic"), labels = c("Cranial", "Pelvic")))

perm_acc_df <- sim_fig %>%
  group_by(region) %>%
  group_split() %>%
  map_dfr(\(df) {
    reg <- unique(df$region)
    tcols <- get_trait_cols(df)
    perm <- perm_ceiling_rfm(df, tcols, n_perm = n_perm, seed = 1234)
    tibble(region = reg, acc = perm$acc)
  })

perm_summ_df <- perm_acc_df %>%
  group_by(region) %>%
  summarise(
    perm_mean = mean(acc),
    perm_p95 = as.numeric(quantile(acc, 0.95)),
    perm_p99 = as.numeric(quantile(acc, 0.99)),
    .groups = "drop"
  )

pA <- ggplot(perm_acc_df, aes(x = acc)) +
  geom_density(linewidth = 0.9) +
  geom_vline(data = perm_summ_df, aes(xintercept = perm_mean), linewidth = 0.9) +
  geom_vline(data = perm_summ_df, aes(xintercept = perm_p95), linetype = "dashed", linewidth = 0.8) +
  geom_vline(data = perm_summ_df, aes(xintercept = perm_p99), linetype = "dotted", linewidth = 0.8) +
  geom_text(data = perm_summ_df, aes(x = perm_mean, y = Inf, label = "Mean"), vjust = 1.5, size = 3) +
  geom_text(data = perm_summ_df, aes(x = perm_p95, y = Inf, label = "95%"), vjust = 1.5, size = 3) +
  geom_text(data = perm_summ_df, aes(x = perm_p99, y = Inf, label = "99%"), vjust = 1.5, size = 3) +
  facet_wrap(~region, ncol = 2) +
  coord_cartesian(xlim = c(0.35, 0.95)) +
  labs(
    title = "A. Trait-signal ceiling (simulated data)",
    x = "Permutation RFM accuracy",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

post_df <- dat_clean %>%
  mutate(
    dataset = factor(dataset, levels = c("simulated", "known"), labels = c("Simulated", "Known-sex")),
    region = factor(region, levels = c("cranial", "pelvic"), labels = c("Cranial", "Pelvic"))
  )

pB <- ggplot(post_df, aes(x = postprob)) +
  geom_histogram(bins = 15) +
  geom_vline(xintercept = 0.95, linetype = "dashed", linewidth = 0.6) +
  facet_grid(dataset ~ region) +
  coord_cartesian(xlim = c(0.45, 1.00)) +
  labs(
    title = "B. MorphoPASSE posterior probabilities",
    x = "Posterior probability for predicted class",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

summ_tbl <- dat_clean %>%
  dplyr::filter(dataset %in% c("simulated", "known")) %>%
  group_by(dataset, region) %>%
  group_split() %>%
  map_dfr(\(df) {
    ds <- unique(df$dataset)
    rg <- unique(df$region)
    summarize_dataset_region(df, ds, rg, n_perm = n_perm)
  }) %>%
  mutate(
    dataset = factor(dataset, levels = c("simulated", "known"), labels = c("Simulated", "Known-sex")),
    region = factor(region, levels = c("cranial", "pelvic"), labels = c("Cranial", "Pelvic"))
  )

points_df <- summ_tbl %>%
  dplyr::select(dataset, region, post_mean, post_median, perm_mean, accuracy) %>%
  pivot_longer(
    cols = c(perm_mean, post_mean, post_median, accuracy),
    names_to = "stat",
    values_to = "value"
  ) %>%
  dplyr::filter(!(stat == "accuracy" & is.na(value))) %>%
  mutate(
    stat = recode(
      stat,
      perm_mean = "Permutation mean",
      post_mean = "MorphoPASSE mean postprob",
      post_median = "MorphoPASSE median postprob",
      accuracy = "Observed accuracy (known only)"
    ),
    label = sprintf("%.2f", value)
  )

summ_tbl2 <- summ_tbl %>%
  mutate(region = factor(region, levels = c("Cranial", "Pelvic")), x_pos = as.numeric(region))

points_df2 <- points_df %>%
  mutate(region = factor(region, levels = c("Cranial", "Pelvic")), x_pos = as.numeric(region)) %>%
  mutate(pt_size = ifelse(stat == "Permutation mean", 2.6, 3.4))

pC <- ggplot(summ_tbl2, aes(x = x_pos)) +
  geom_hline(yintercept = 0.5, linewidth = 0.9, linetype = "dashed", color = "firebrick") +
  annotate(
    "curve",
    x = 1.35, y = 0.60, xend = 1.02, yend = 0.50,
    curvature = -0.35,
    arrow = arrow(length = unit(0.16, "inches")),
    linewidth = 0.7,
    color = "firebrick"
  ) +
  annotate(
    "text",
    x = 1.38, y = 0.62,
    label = "Chance = 0.50",
    hjust = 0,
    size = 3.4,
    color = "firebrick"
  ) +
  geom_linerange(aes(ymin = perm_p95, ymax = perm_p99), linewidth = 10, alpha = 0.18) +
  geom_point(aes(y = perm_p95), shape = 95, size = 10, alpha = 0.55) +
  geom_text(
    data = summ_tbl2 %>%
      group_by(dataset) %>%
      summarise(x_pos = 1.05, y = mean(c(mean(perm_p95), mean(perm_p99))), .groups = "drop"),
    aes(x = x_pos, y = y, label = "Permutation benchmark\n(95–99%)"),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.3
  ) +
  geom_point(
    data = points_df2,
    aes(x = x_pos, y = value, shape = stat, size = pt_size),
    stroke = 1.1,
    inherit.aes = FALSE
  ) +
  scale_size_identity() +
  geom_text(
    data = points_df2 %>% dplyr::filter(stat != "Permutation mean"),
    aes(x = x_pos, y = value, label = label),
    nudge_x = 0.18,
    vjust = -0.5,
    size = 3.3,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  facet_wrap(~dataset) +
  coord_cartesian(ylim = c(0.45, 1.02), clip = "off") +
  scale_x_continuous(breaks = c(1, 2), labels = c("Cranial", "Pelvic")) +
  scale_shape_manual(values = c(
    "Permutation mean" = 16,
    "MorphoPASSE mean postprob" = 17,
    "MorphoPASSE median postprob" = 15,
    "Observed accuracy (known only)" = 1
  )) +
  labs(
    title = "C. Permutation benchmark vs posterior confidence",
    x = NULL,
    y = "Value",
    shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

multi <- (pA / pB / pC) + plot_layout(heights = c(1.0, 1.15, 1.1))

print(multi)

ggsave(
  filename = "figures/generated/Figure1_MultiPanel_A-C.png",
  plot = multi,
  width = 11,
  height = 13,
  dpi = 300
)
