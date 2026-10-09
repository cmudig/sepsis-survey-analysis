### This script contains the prespecified analyses from osf.io/h5urm. Run
### 01_preprocessing.ipynb in Python before executing this script.

install_if_missing <- function(pkgs){
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if(length(to_install)) install.packages(to_install)
}
install_if_missing(c("AICcmodavg", "ordinal", "emmeans", "broom.mixed", "DHARMa", "dplyr", "ggplot2", "nnet", "car", "effectsize", "multcomp"))

library(ordinal)     # clmm
library(emmeans)     # predicted probabilities and contrasts
library(broom.mixed) # tidy() for mixed models
library(DHARMa)      # diagnostics
library(dplyr)
library(ggplot2)
library(nnet)
library(scales)
library(car)
library(effectsize)
library(multcomp)
library(AICcmodavg)
library(lme4)

# Change this to the date from the preprocessing script
analysis_date <- "261008"

decision_df <- read.csv(paste0("data/cleaned_decision_level_", analysis_date, ".csv"), stringsAsFactors = TRUE)

respondent_df <- subset(read.csv(paste0("data/cleaned_respondent_level_", analysis_date, ".csv"), stringsAsFactors = TRUE), 
                        select = -AIInterface)
decision_df <- merge(decision_df, respondent_df, by = 'ResponseId')

decision_df$Case <- relevel(factor(decision_df$Case), ref = "Esther") # arbitrary

decision_df$QualityScore <- ordered(decision_df$QualityScore, levels = sort(unique(as.integer(decision_df$QualityScore))))
decision_df$Confidence <- ordered(decision_df$Confidence, levels = c("Not at all confident", "Slightly confident", "Somewhat confident", "Fairly confident", "Extremely confident"))
decision_df$HighConfidence <- decision_df$Confidence == "Fairly confident" |  decision_df$Confidence == "Extremely confident"
decision_df$HighPerceivedQuality <- decision_df$AIQuality == "Satisfied" |  decision_df$AIQuality == "Very satisfied"
decision_df$ExperiencedFlag <- decision_df$Experience == 5
decision_df$ExpertFlag <- decision_df$SepsisExpertise == 5
decision_df$ClinicalTime <- relevel(factor(decision_df$clinical_time), ref = "Full")
decision_df$PracticeSettingCommunity <- factor(decision_df$practice_setting == "Community")
decision_df$Specialty <- factor(decision_df$specialty)

# Condition: factor with control "0" as reference, ensure levels sorted 0..5
decision_df$AIInterface <- factor(decision_df$AIInterface, levels = c("No AI", "Recommended Treatment", "Options to Consider", "Options to Avoid"))
decision_df$AIInterfaceBinary <- factor(decision_df$AIInterface != "No AI")

# Quick sanity checks
cat("Rows:", nrow(decision_df), "\n")
cat("Condition levels:", paste(levels(decision_df$AIInterface), collapse = ", "), "\n")
table(decision_df$QualityScore, useNA = "ifany")


### 1. Decision Quality Analysis

model_cond_only <- clmm(QualityScore ~ AIInterface + (1 | Case) + (1 | ResponseId),
                        data = decision_df,
                        link = "logit",
                        Hess = TRUE)

model_cond_binary <- clmm(QualityScore ~ AIInterfaceBinary + (1 | Case) + (1 | ResponseId),
                          data = decision_df,
                          link = "logit",
                          Hess = TRUE)

model_bg <- clmm(QualityScore ~ SelfReportedUnaffected + (1 | Case) + (1 | ResponseId),
                 data = decision_df,
                 link = "logit",
                 Hess = TRUE)

model_demographics <- clmm(QualityScore ~ ExperiencedFlag + ExpertFlag + PracticeSettingCommunity + ClinicalTime + Specialty + LongTime + (1 | Case) + (1 | ResponseId),
                           data = decision_df,
                           link = "logit",
                           Hess = TRUE)

model_cond_influences <- clmm(QualityScore ~ AIInterface * SelfReportedUnaffected + (1 | Case) + (1 | ResponseId),
                              data = decision_df,
                              link = "logit",
                              Hess = TRUE)

model_cond_time <- clmm(QualityScore ~ AIInterface * LongTime + (1 | Case) + (1 | ResponseId),
                        data = decision_df,
                        link = "logit",
                        Hess = TRUE)

model_simple <- clmm(QualityScore ~ (1 | Case) + (1 | ResponseId),
                     data = decision_df,
                     link = "logit",
                     Hess = TRUE)

model_list <- list(Demographics = model_demographics,
                   ConditionOnly = model_cond_only,
                   ConditionBinary = model_cond_binary,
                   ConditionInfluences = model_cond_influences,
                   ConditionTime = model_cond_time,
                   BehaviorGroup = model_bg,
                   NoFixedEffects = model_simple)

aictab(cand.set = model_list)

print("Comparison between simple model and condition-only model:")
print(anova(model_cond_only, model_simple))

print("Odds ratios:")
# --- 4) Extract fixed-effect estimates and compute odds ratios + 95% CIs ---
s <- coef(summary(model_cond_only))
# s rows include thresholds (alpha) and fixed effects. We'll keep only fixed effects by excluding threshold names:
threshold_names <- names(model_cond_only$alpha)
fixed_rows <- rownames(s)[!(rownames(s) %in% threshold_names)]
fixed_table <- s[fixed_rows, , drop = FALSE]

est <- fixed_table[, "Estimate"]
se  <- fixed_table[, "Std. Error"]
OR  <- exp(est)
z <- qnorm(0.975)
ci_low <- exp(est - z * se)
ci_high <- exp(est + z * se)

fixed_results <- data.frame(
  term = rownames(fixed_table),
  estimate = est,
  SE = se,
  OR = OR,
  CI_low = ci_low,
  CI_high = ci_high,
  row.names = NULL,
  stringsAsFactors = FALSE
)
print(fixed_results, digits = 4)
write.csv(fixed_results, "result_data/prespecified_decision_quality_ors.csv")

print("Test that all three weight coefficients are equal:")
coef_names <- c("AIInterfaceRecommended Treatment", "AIInterfaceOptions to Consider", "AIInterfaceOptions to Avoid")

# extract coefficient vector and covariance matrix for those terms
beta_all <- coef(model_cond_only)
vcov_all <- vcov(model_cond_only)

L <- rbind(
  c(1, -1, 0),
  c(1,  0, -1)
)
colnames(L) <- coef_names

beta <- beta_all[coef_names]
V <- vcov_all[coef_names, coef_names, drop = FALSE]

# basic sanity checks
if (any(is.na(beta)) || any(is.na(V))) stop("One or more names in coef_names not found in model. Check names(coef(model_cond_only)).")
if (nrow(V) != length(beta)) stop("Covariance submatrix has wrong dimensions.")

Lb <- as.numeric(L %*% beta)
LVL <- L %*% V %*% t(L)

W <- t(Lb) %*% solve(LVL) %*% Lb
df <- nrow(L)
p_value <- 1 - pchisq(W, df)

cat("Wald chi-square =", format(W, digits = 6), "on", df, "df; p =", format(p_value, digits = 6), "\n")

## Visualize the condition-only model

emm_latent <- emmeans(model_cond_only, ~ AIInterface, mode = "latent", bias.adjust = TRUE)

# pairwise contrasts on latent scale
ctr <- contrast(emm_latent, method = "pairwise", adjust = "tukey")
print(ctr)

# turn to data.frame and compute ORs + 95% CIs by delta method (exp on estimate ± 1.96*SE)
ctr_df <- as.data.frame(ctr)
ctr_df$OR <- exp(ctr_df$estimate)
ctr_df$lower.OR <- exp(ctr_df$estimate - 1.96 * ctr_df$SE)
ctr_df$upper.OR <- exp(ctr_df$estimate + 1.96 * ctr_df$SE)

# tidy output
ctr_df <- ctr_df[, c("contrast", "estimate", "SE", "OR", "lower.OR", "upper.OR", "z.ratio", "p.value")]
print(ctr_df)

emm_df <- as.data.frame(emm_latent)
write.csv(emm_df, "result_data/prespecified_decision_quality_emmeans.csv")

write.csv(ctr_df, "result_data/prespecified_decision_quality_pairs.csv")

write.csv(as.data.frame(model_cond_only$Theta), 
          "result_data/prespecified_decision_quality_thresholds.csv")

### 2. Satisfaction

respondent_df <- read.csv(paste0("data/cleaned_respondent_level_", analysis_date, ".csv"), stringsAsFactors = TRUE)

# Fit one-way ANOVA
responses <- c("Total.Satisfaction", "Productivity", "Effectiveness", "EaseOfUse", "HighOutputQuality", "WouldUse")  # put your variables here

p_values <- numeric()

for (resp in responses) {
  cat("====================================\n")
  cat("Response variable:", resp, "\n")
  
  formula <- as.formula(paste(resp, "~ AIInterface"))
  aov_model <- aov(formula, data = respondent_df)
  
  print(summary(aov_model))
  p_values <- append(p_values, summary(aov_model)[[1]][["Pr(>F)"]][1])

  # Numerical assumption checks
  cat("\n=== Shapiro-Wilk test for normality of residuals ===\n")
  shapiro_res <- shapiro.test(residuals(aov_model))
  print(shapiro_res)
  
  cat("\n=== Levene's test for homogeneity of variances (center = median) ===\n")
  levene_res <- car::leveneTest(formula, data = respondent_df, center = median)
  print(levene_res)
  
  # Effect size (eta-squared, partial = FALSE for classical one-way)
  cat("\n=== Effect size: eta-squared ===\n")
  eta2_res <- effectsize::eta_squared(aov_model, partial = FALSE)
  print(eta2_res)
  
  # Estimated marginal means and pairwise contrasts (Holm adjustment)
  emm <- emmeans(aov_model, ~ AIInterface)
  
  write.csv(as.data.frame(emm), paste("result_data/prespecified_satisfaction_", resp, "_emmeans.csv", sep = ""))
  write.csv(as.data.frame(pairs(emm, adjust = "holm")), paste("result_data/prespecified_satisfaction_", resp, "_emmeans_pairs.csv", sep = ""))
}

print(p.adjust(p_values[1:length(p_values)], method = "holm"))

respondent_df$ExperiencedFlag <- respondent_df$Experience == 5
respondent_df$ExpertFlag <- respondent_df$SepsisExpertise == 5

respondent_df$ClinicalTime <- relevel(factor(respondent_df$clinical_time), ref = "Full")
respondent_df$PracticeSettingCommunity <- factor(respondent_df$practice_setting == "Community")
respondent_df$Specialty <- factor(respondent_df$specialty)

model_full <- lm(Total.Satisfaction ~ AIInterface + ExpertFlag + ExperiencedFlag + SelfReportedUnaffected + PracticeSettingCommunity + ClinicalTime + Specialty,
                 data = respondent_df)

model_simple <- lm(Total.Satisfaction ~ AIInterface,
                   data = respondent_df)

print(anova(model_full, model_simple, test = "LRT"))

print(summary(model_full))
