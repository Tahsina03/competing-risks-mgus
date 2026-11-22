# Competing risks analysis of progression to plasma cell malignancy (PCM)
# in patients with MGUS, with death before progression as a competing event.
#
# Data: survival::mgus2 (Mayo Clinic, n = 1384)
# Methods: Aalen-Johansen cumulative incidence, Gray's test,
#          cause-specific Cox models, Fine-Gray subdistribution hazard models,
#          and a comparison with the naive 1 - Kaplan-Meier estimate.
#
# Run from the repository root:  Rscript R/analysis.R

suppressPackageStartupMessages({
  library(survival)
  library(cmprsk)
  library(ggplot2)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

# ---------------------------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------------------------
d <- survival::mgus2

# Time to first event: PCM if it happened, otherwise death or censoring
d$etime  <- with(d, ifelse(pstat == 1, ptime, futime))
d$status <- with(d, ifelse(pstat == 1, 1, 2 * death))   # 0 = censored, 1 = PCM, 2 = death
d$event  <- factor(d$status, levels = 0:2, labels = c("censor", "PCM", "death"))
d$etime_yr <- d$etime / 12

# Covariates: age (per 10 years), sex, M-spike size
d$age10 <- d$age / 10
d <- d[complete.cases(d[, c("etime", "event", "age10", "sex", "mspike")]), ]

summary_tab <- data.frame(
  n          = nrow(d),
  PCM        = sum(d$event == "PCM"),
  death      = sum(d$event == "death"),
  censored   = sum(d$event == "censor"),
  median_age = median(d$age),
  pct_male   = round(100 * mean(d$sex == "M"), 1)
)
write.csv(summary_tab, "output/tables/sample_summary.csv", row.names = FALSE)
print(summary_tab)

# ---------------------------------------------------------------------------
# 2. Cumulative incidence (Aalen-Johansen) vs naive 1 - KM
# ---------------------------------------------------------------------------
aj <- survfit(Surv(etime_yr, event) ~ 1, data = d)
km <- survfit(Surv(etime_yr, event == "PCM") ~ 1, data = d)

grid <- seq(0, 30, by = 0.25)
s_aj <- summary(aj, times = grid, extend = TRUE)
s_km <- summary(km, times = grid, extend = TRUE)
pcm_col <- which(aj$states == "PCM")
death_col <- which(aj$states == "death")

cif_df <- rbind(
  data.frame(time = grid, risk = s_aj$pstate[, pcm_col],   curve = "PCM (Aalen-Johansen)"),
  data.frame(time = grid, risk = s_aj$pstate[, death_col], curve = "Death before PCM (Aalen-Johansen)"),
  data.frame(time = grid, risk = 1 - s_km$surv,            curve = "PCM (naive 1 - KM)")
)

p1 <- ggplot(cif_df, aes(time, risk, linetype = curve, colour = curve)) +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(values = c("PCM (Aalen-Johansen)" = "#1b6ca8",
                                 "Death before PCM (Aalen-Johansen)" = "grey45",
                                 "PCM (naive 1 - KM)" = "#c0392b")) +
  scale_linetype_manual(values = c("PCM (Aalen-Johansen)" = "solid",
                                   "Death before PCM (Aalen-Johansen)" = "solid",
                                   "PCM (naive 1 - KM)" = "dashed")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Years since MGUS diagnosis", y = "Cumulative incidence",
       colour = NULL, linetype = NULL,
       title = "Ignoring the competing risk of death overstates PCM risk") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave("output/figures/cif_vs_km.png", p1, width = 7, height = 5.5, dpi = 150, bg = "white")

at <- c(5, 10, 20)
risk_tab <- data.frame(
  years       = at,
  PCM_AJ      = round(100 * summary(aj, times = at)$pstate[, pcm_col], 1),
  PCM_naiveKM = round(100 * (1 - summary(km, times = at)$surv), 1),
  Death_AJ    = round(100 * summary(aj, times = at)$pstate[, death_col], 1)
)
write.csv(risk_tab, "output/tables/cumulative_incidence.csv", row.names = FALSE)
print(risk_tab)

# ---------------------------------------------------------------------------
# 3. CIF by sex and Gray's test
# ---------------------------------------------------------------------------
gray <- cmprsk::cuminc(ftime = d$etime_yr, fstatus = d$status, group = d$sex, cencode = 0)
gray_tab <- data.frame(event = c("PCM", "death"),
                       statistic = round(gray$Tests[, "stat"], 2),
                       p_value = signif(gray$Tests[, "pv"], 3))
write.csv(gray_tab, "output/tables/gray_test_by_sex.csv", row.names = FALSE)
print(gray_tab)

aj_sex <- survfit(Surv(etime_yr, event) ~ sex, data = d)
s_sex <- summary(aj_sex, times = grid, extend = TRUE)
sex_df <- data.frame(
  time  = s_sex$time,
  sex   = sub("sex=", "", s_sex$strata),
  PCM   = s_sex$pstate[, pcm_col],
  death = s_sex$pstate[, death_col]
)
sex_long <- rbind(
  data.frame(sex_df[, c("time", "sex")], risk = sex_df$PCM,   event = "PCM"),
  data.frame(sex_df[, c("time", "sex")], risk = sex_df$death, event = "Death before PCM")
)
sex_long$sex <- factor(sex_long$sex, levels = c("F", "M"), labels = c("Female", "Male"))
p2 <- ggplot(sex_long, aes(time, risk, colour = sex)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~ event, scales = "free_y") +
  scale_colour_manual(values = c("#8e44ad", "#16a085")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Years since MGUS diagnosis", y = "Cumulative incidence", colour = NULL,
       title = "Cumulative incidence by sex") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/cif_by_sex.png", p2, width = 8, height = 4.5, dpi = 150, bg = "white")

# ---------------------------------------------------------------------------
# 4. Regression: cause-specific hazards vs Fine-Gray subdistribution hazards
# ---------------------------------------------------------------------------
cs_pcm   <- coxph(Surv(etime_yr, event == "PCM")   ~ age10 + sex + mspike, data = d)
cs_death <- coxph(Surv(etime_yr, event == "death") ~ age10 + sex + mspike, data = d)

fg_pcm_data   <- finegray(Surv(etime_yr, event) ~ ., data = d, etype = "PCM")
fg_death_data <- finegray(Surv(etime_yr, event) ~ ., data = d, etype = "death")
fg_pcm   <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age10 + sex + mspike,
                  weights = fgwt, data = fg_pcm_data)
fg_death <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age10 + sex + mspike,
                  weights = fgwt, data = fg_death_data)

tidy_hr <- function(fit, model, outcome) {
  s  <- summary(fit)$coefficients
  ci <- exp(confint(fit))
  data.frame(outcome = outcome, model = model, term = rownames(s),
             HR = round(exp(s[, "coef"]), 2),
             lower = round(ci[, 1], 2), upper = round(ci[, 2], 2),
             p = signif(s[, ncol(s)], 3), row.names = NULL)
}
hr_tab <- rbind(
  tidy_hr(cs_pcm,   "Cause-specific HR",  "PCM"),
  tidy_hr(fg_pcm,   "Subdistribution HR", "PCM"),
  tidy_hr(cs_death, "Cause-specific HR",  "Death"),
  tidy_hr(fg_death, "Subdistribution HR", "Death")
)
hr_tab$term <- factor(hr_tab$term, levels = c("age10", "sexM", "mspike"),
                      labels = c("Age (per 10 years)", "Male sex", "M-spike (g/dL)"))
write.csv(hr_tab, "output/tables/hazard_ratios.csv", row.names = FALSE)
print(hr_tab)

p3 <- ggplot(hr_tab, aes(HR, term, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ outcome) +
  scale_x_log10() +
  scale_colour_manual(values = c("#1b6ca8", "#e67e22")) +
  labs(x = "Hazard ratio (log scale, 95% CI)", y = NULL, colour = NULL,
       title = "Cause-specific vs subdistribution hazard ratios") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/hazard_ratios.png", p3, width = 8, height = 4.5, dpi = 150, bg = "white")

# Proportional hazards check for the cause-specific PCM model
ph <- cox.zph(cs_pcm)
write.csv(data.frame(term = rownames(ph$table), round(ph$table, 3)),
          "output/tables/ph_test_cs_pcm.csv", row.names = FALSE)
print(ph)

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat("\nDone. Results written to output/\n")
