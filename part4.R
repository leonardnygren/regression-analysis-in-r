library(pROC)
library(caret)
library(tidyverse)
library(readxl)
library(broom)
library(dplyr)
carotene <- read_excel("/data/carotene.xlsx")


carotene$sex <- factor(carotene$sex, levels = c("1", "2"), labels = c("Male", "Female"))

carotene$vituse <- factor(carotene$vituse, levels = c("1", "2", "3"),
                          labels = c("Regularly", "Irregularly", "No"))

carotene$smokstat <- factor(carotene$smokstat, levels = c("1", "2", "3"),
                            labels = c("Never", "Former", "Current Smoker"))

carotene$vituse <- relevel(carotene$vituse, ref = "No")

carotene <- carotene |>
  mutate(vituse_red= factor (if_else(vituse %in% c("Regularly", "Irregularly"),
                                         "User", "No"),
                                 levels = c("No", "User")))
carotene$vituse_red <- relevel(carotene$vituse_red, ref = "No")
a = 0.42*536.88
carotene |> mutate(
  lowplasma_01 = as.numeric(betaplasma < a),
  lowplasma_hl = factor(lowplasma_01,
  levels = c(0, 1),
  labels = c("high", "low"))) -> carotene

# Stepwise from Backward model selected in Section 3
stepwise_from_backward <- glm(
  lowplasma_01 ~ bmi + age + fat + cholesterol + fiber + betadiet +
    fat:betadiet + cholesterol:betadiet + fiber:betadiet +
    bmi:betadiet,
  data = carotene,
  family = binomial
)

AIC(stepwise_from_backward)
BIC(stepwise_from_backward)
coef(stepwise_from_backward)

summary(stepwise_from_backward)


summary(backward_model)

reduce_vituse_model <- glm(lowplasma_01 ~ bmi + age +  betadiet +  vituse_red +
                        betadiet:bmi,
                      data = carotene,
                      family = binomial)

summary(reduce_vituse_model)

pred_phat <- carotene |>
  mutate(
    p_b = predict(stepwise_from_backward, type = "response"),
    p_v = predict(reduce_vituse_model, type = "response")
  )

pred_phat <- pred_phat |> mutate(
  actual_lowplasma = factor(lowplasma_01,
                            levels = c(0, 1),
                            labels = c("High", "Low"))
)

pred_phat |> mutate(
  yhat_b = factor(p_b > 0.5,
                     levels = c(FALSE, TRUE),
                     labels = c( "High","Low")),
    yhat_v = factor(p_v > 0.5,
                     levels = c(FALSE, TRUE),
                     labels = c( "High","Low"))
  
  ) -> pred_phat




confusion_matrix_b <- confusionMatrix(
  data = pred_phat$yhat_b, 
  reference = pred_phat$actual_lowplasma,
  positive = "Low")

confusion_matrix_v <- confusionMatrix(
  data = pred_phat$yhat_v, 
  reference = pred_phat$actual_lowplasma,
  positive = "Low")

confusion_matrix_b
confusion_matrix_v

library(pROC)
roc_b <- roc(actual_lowplasma ~ p_b, data = pred_phat)
roc_v <- roc(actual_lowplasma ~ p_v, data = pred_phat)

#coords(roc_oslo) |> head()
png("../Figures/4/ROC.png", width=8, height=6, units="in", res=600)
plot(roc_b, col = "blue", lwd = 2, main = "ROC Curves for Stepwise Bwd vs Reduced Vituse")
plot(roc_v, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c("Stepwise Bwd", "Reduced Vituse Model"),
       col = c("blue", "red"), lwd = 2)
dev.off()
#ggsave("../Figures/4/ROC.png", width = 8, height = 6, dpi = 600)
auc_b <-auc(roc_b)
auc_v <-auc(roc_v)
ci_b <-ci(roc_b)
ci_v <-ci(roc_v)
data.frame(
  Model = c("Stepwise Bwd"," Reduced Vituse"),
  AUC = c(ci_b[2],ci_v[2]),
  LowerCI = c(ci_b[1],ci_v[1]),
  UpperCI = c(ci_b[3],ci_v[3])
)

roc.test(roc_b, roc_v)

topleft_b <- coords(roc_b, "best", best.method = "closest.topleft")
topleft_b
topleft_v <- coords(roc_v, "best", best.method = "closest.topleft")
topleft_v


pred_phat |> mutate(
  yhat_b_optimal = factor(p_b > as.numeric(topleft_b["threshold"]),
                     levels = c(FALSE, TRUE),
                     labels = c( "High","Low")),
    yhat_v_optimal = factor(p_v > as.numeric(topleft_v["threshold"]),
                     levels = c(FALSE, TRUE),
                     labels = c( "High","Low"))
  
  ) -> pred_phat


confusion_matrix_b_optimal <- confusionMatrix(
  data = pred_phat$yhat_b_optimal, 
  reference = pred_phat$actual_lowplasma,
  positive = "Low")

confusion_matrix_v_optimal <- confusionMatrix(
  data = pred_phat$yhat_v_optimal, 
  reference = pred_phat$actual_lowplasma,
  positive = "Low")

confusion_matrix_b_optimal
confusion_matrix_v_optimal
exp(cbind(OddsRatio = coef(reduce_vituse_model), 
           confint(reduce_vituse_model)))


