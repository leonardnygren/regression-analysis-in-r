library(readxl)
library(tidyverse); theme_set(theme_bw() + theme(text = element_text(size = 18)))
library(car)
library(rstatix)
library(GGally)
library(patchwork)
library(plyr)
library(dplyr)
library(tidyr)


rm(list=ls()) # i.e. clear all

carotene <- read_xlsx("data/carotene.xlsx")

#turning categorical variables into factor variables with suitable references
mutate(carotene,
       sex = factor(sex,
                    levels = c(1, 2),
                    labels = c("Male", "Female"))) -> carotene
carotene$sex = relevel(carotene$sex, ref = "Female")

mutate(carotene,
       vituse = factor(vituse,
                       levels = c(1, 2, 3),
                       labels = c("Regularly", "Irregularly", "No"))) -> carotene
carotene$vituse = relevel(carotene$vituse, ref = "No")

mutate(carotene,
       smokstat = factor(smokstat,
              levels = c(1,2,3),
              labels = c("Never","Former","Current Smoker"))) -> carotene
carotene$smokstat = relevel(carotene$smokstat, ref = "Never")
glimpse(carotene)



#Section 1 
#1(a)
a <-225.5

#1(b)
carotene |> mutate(
  lowplasma_01 = as.numeric(betaplasma < a),
  lowplasma_hl = factor(lowplasma_01,
                        levels = c(0, 1),
                        labels = c("high", "low"))) -> carotene

carotene$vituse
glimpse(carotene)

#1(c)
carotene |> 
  count(vituse, lowplasma_hl) |> 
  pivot_wider(id_cols = vituse,
              names_from = lowplasma_hl,
              values_from = n)
names(carotene)
carotene

#calculating probabilities and odds
table3c <- tibble::tribble(
  ~vituse,       ~high, ~low,
  "No",           13,    98,
  "Regularly",    46,    76,
  "Irregularly",  21,    61
)

table3c <- table3c |>
  mutate(
    total = high + low,
    prob_low = low / total,
    odds_low = low / high
  )

# Get reference odds
ref_odds <- table3c |> filter(vituse == "No") |> pull(odds_low)

# Add odds ratio
table3c <- table3c |>
  mutate(
    odds_ratio = odds_low / ref_odds
  )

table3c

#1(d)
#create model
model_1d <- glm(lowplasma_hl ~ vituse, 
                family = "binomial", 
                data = carotene
                )

sum_model_1d = summary(model_1d)
sum_model_1d

exp(sum_model_1d$coefficients)

#conf-int for beta-params and e^beta-params
conf1d = confint(model_1d)
conf1d

conf1dexp = exp(conf1d)
conf1dexp

#1(d) pt 2: linear predictor, etc

new_data <- data.frame(vituse = levels(carotene$vituse))

# logit linear predictor and standard error
pred <- predict(model_1d, newdata = new_data, type = "link", se.fit = TRUE)

# Compute confidence intervals on the logit (linear predictor) scale
logit_fit <- pred$fit
logit_se <- pred$se.fit
logit_lower <- logit_fit - 1.96 * logit_se
logit_upper <- logit_fit + 1.96 * logit_se

# Transform logit to probability
prob_fit <- plogis(logit_fit)
prob_lower <- plogis(logit_lower)
prob_upper <- plogis(logit_upper)

# Compute odds
odds_fit <- exp(logit_fit)
odds_lower <- exp(logit_lower)
odds_upper <- exp(logit_upper)

# Combine everything
results <- data.frame(
  vituse = new_data$vituse,
  logit = logit_fit,
  logit_lower = logit_lower,
  logit_upper = logit_upper,
  odds = odds_fit,
  odds_lower = odds_lower,
  odds_upper = odds_upper,
  probability = prob_fit,
  prob_lower = prob_lower,
  prob_upper = prob_upper
)

print(results)



#1(e)
#global likelihood test with H_0: beta_1=beta_2=0
D_diff <- sum_model_1d$null.deviance - sum_model_1d$deviance
df_diff <- sum_model_1d$df.null - sum_model_1d$df.residual
chi2_alpha <- qchisq(p = 1 - 0.05, df = df_diff)
Pvalue <- pchisq(q = D_diff, df = df_diff, lower.tail = FALSE)

cbind(D_diff, df_diff, chi2_alpha, Pvalue)



#Section 3
#3(a)
model_temp = glm(lowplasma_hl ~ vituse +  bmi + age + calories + fat + 
                 cholesterol + fiber + alcohol + betadiet + smokstat + 
                 sex, family = "binomial", data=carotene)

vif(model_temp)

model_3a = glm(lowplasma_hl ~ vituse +  bmi + age + fat + 
                 cholesterol + fiber + alcohol + betadiet + smokstat + 
                 sex, family = "binomial", data=carotene)
vif(model_3a) 

#3(b) 
#fit model w/ all covariates & their two-way interactions 
#Full model
full_model = glm(lowplasma_hl ~ (vituse +  bmi + age + fat + 
                 cholesterol + fiber + alcohol + betadiet + smokstat + 
                 sex)^2, family = "binomial", data=carotene)

#Null model
null_model = glm(lowplasma_hl ~1, family = "binomial", data = carotene)


#Backward elimination
#note: "k = log(nrow(carotene))" changes it from AIC to BIC,even if it says "AIC"
#the best model (according to the BIC criterion) is stored in backward_model
backward_model <- step(full_model, direction = "backward", 
                       scope = list(lower = formula(null_model), 
                                    upper = formula(full_model)), 
                       k = log(nrow(carotene)))
summary(backward_model)

#Forward selection
forward_model <- step(null_model, direction = "forward", scope = list(lower = formula(null_model), upper = formula(full_model)), k = log(nrow(carotene)))
summary(forward_model)

#Stepwise regression
stepwise_from_backward <- step(backward_model, direction = "both", scope = list(lower = formula(null_model), upper = formula(full_model)), k = log(nrow(carotene)))

stepwise_from_forward <- step(forward_model, direction = "both", scope = list(lower = formula(null_model), upper = formula(full_model)), k = log(nrow(carotene)))

summary(stepwise_from_backward)
summary(stepwise_from_forward) #contains categorical variable vituse


#test if we can reduce the stepwise forward model to only use two categories

#binary vituse version
carotene$vituse_reduced <- factor(carotene$vituse, levels = c("No", "Irregularly", "Regularly"))
levels(carotene$vituse_reduced)[levels(carotene$vituse_reduced) != "No"] <- "User"

#reduced model (w binary vituse)
model_reduced_vituse <- glm(lowplasma_hl ~ betadiet + bmi + age + vituse_reduced + betadiet:bmi, 
                            family = binomial, data = carotene)

summary(model_reduced_vituse)

#compare models with likelihood ratio test 
anova(stepwise_from_forward, model_reduced_vituse, test = "Chisq")
#reduced model is fine (fail to reject H_0)


#create table 3(b)
df_backward <- data.frame(variable = names(coef(backward_model)),
                          b_backward = coef(backward_model), row.names = NULL)

df_forward <- data.frame(variable = names(coef(forward_model)),
                         b_forward = coef(forward_model), row.names = NULL)

df_step_f <- data.frame(variable = names(coef(stepwise_from_forward)),
                        b_step_f = coef(stepwise_from_forward), row.names = NULL)

df_step_b <- data.frame(variable = names(coef(stepwise_from_backward)),
                        b_step_b = coef(stepwise_from_backward), row.names = NULL)

df_reduced <- data.frame(variable = names(coef(model_reduced_vituse)),
                         b_reduced = coef(model_reduced_vituse), row.names = NULL)

table_3b <- df_backward %>%
  full_join(df_forward, by = "variable") %>%
  full_join(df_step_f, by = "variable") %>%
  full_join(df_step_b, by = "variable") %>%
  full_join(df_reduced, by = "variable")

print(table_3b)

#note that forward and stepwise fwd are the same



#3(c)
aic <- AIC(backward_model, forward_model, stepwise_from_backward, model_reduced_vituse)

bic <- BIC(backward_model, forward_model, stepwise_from_backward, model_reduced_vituse)

collect.AICetc <- data.frame(aic, bic)
collect.AICetc |> mutate(df.1 = NULL) -> collect.AICetc

lnL0 <- logLik(null_model)[1]

collect.AICetc |> mutate(
  loglik =  c(logLik(backward_model)[1],
              logLik(forward_model)[1],
              logLik(stepwise_from_backward)[1],
              logLik(model_reduced_vituse)[1])) -> collect.AICetc

collect.AICetc |> mutate(
  R2McF = 1 - loglik/lnL0,
  R2McF.adj = 1 - (loglik - (df - 1)/2)/lnL0) -> collect.AICetc
collect.AICetc


#3(d)
#stepwise_from_backward
#model_reduced_vituse

#linear predictors
pred1 <- predict(stepwise_from_backward)
pred2 <- predict(model_reduced_vituse)

#standardized deviance for stepwise from backward model and reduced vituse model
stepwise_from_backward_infl <- influence(stepwise_from_backward)
model_reduced_vituse_infl <- influence(model_reduced_vituse)

v1 = stepwise_from_backward_infl$hat
v2 = model_reduced_vituse_infl$hat

devresid1 = stepwise_from_backward_infl$dev.res
std_devresid1 = devresid1/sqrt(1-v1)

devresid2 = model_reduced_vituse_infl$dev.res
std_devresid2 = devresid2/sqrt(1-v2)

#plots
# Set up side-by-side plotting layout
par(mfrow = c(2, 2))

# === Stepwise from Backward Model ===

# 1. Residuals vs Linear Predictor
plot(pred1, std_devresid1,
     main = "Std Deviance Residuals vs Linear Predictor\n(Stepwise Model)",
     xlab = "Linear Predictor",
     ylab = "Std. Deviance Residuals",
     pch = 20, col = "blue")
abline(h = 0, col = "red", lty = 2)
abline(h = c(-2, 2), col = "gray", lty = 3)

# 2. QQ-Plot
qqnorm(std_devresid1,
       main = "QQ-Plot of Std. Deviance Residuals\n(Stepwise Model)",
       pch = 20, col = "blue")
qqline(std_devresid1, col = "red")

# === Reduced Vituse Model ===

# 3. Residuals vs Linear Predictor
plot(pred2, std_devresid2,
     main = "Std Deviance Residuals vs Linear Predictor\n(Reduced Vituse Model)",
     xlab = "Linear Predictor",
     ylab = "Std. Deviance Residuals",
     pch = 20, col = "darkgreen")
abline(h = 0, col = "red", lty = 2)
abline(h = c(-2, 2), col = "gray", lty = 3)

# 4. QQ-Plot
qqnorm(std_devresid2,
       main = "QQ-Plot of Std. Deviance Residuals\n(Reduced Vituse Model)",
       pch = 20, col = "darkgreen")
qqline(std_devresid2, col = "red")

