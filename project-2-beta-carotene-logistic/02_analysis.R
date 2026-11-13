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

a <- 0.42 * 536.88

carotene <- carotene |>
  mutate(
    lowplasma_01 = as.numeric(betaplasma < a),
    lowplasma_hl = factor(
      lowplasma_01,
      levels = c(0, 1),
      labels = c("high", "low")
    )
  )

ggplot(carotene, aes(bmi, lowplasma_01)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  xlab("BMI") +
  ylab("Low plasma β-carotene (0 = No, 1 = Yes)") +
  labs(title = "Low plasma β-carotene vs BMI")




model_bmi_glm <- glm(lowplasma_01 ~ bmi, 
                       family = "binomial", data = carotene)
tidy(model_bmi_glm, exponentiate = FALSE, conf.int = TRUE)
tidy(model_bmi_glm, exponentiate = TRUE, conf.int = TRUE)
beta <- coef(model_bmi_glm)["bmi"]
CI <- confint(model_bmi_glm)["bmi", ]
OR_plus1 <- exp(beta)           
OR_minus1 <- exp(-beta)         
OR_minus10 <- exp(-10 *beta)  
OR_results <- data.frame(
  Change_in_BMI = c("+1 unit", "-1 unit", "-10 units"),
  OR = c(exp(beta), exp(-beta), exp(-10 * beta)),
  Lower95CI = c(exp(CI[1]), exp(-CI[2]), exp(-10 * CI[2])),
  Upper95CI = c(exp(CI[2]), exp(-CI[1]), exp(-10 * CI[1]))
)
OR_results
carotene |> mutate(phat = predict(model_bmi_glm, type = "response")) -> carotene_pred



ggplot(carotene_pred, aes(bmi, lowplasma_01)) +
  geom_point() +
  geom_smooth(se = FALSE, linetype = "dashed",span=0.8) +
  geom_line(aes(y = phat), color = "red", size = 1) +
  xlab("bmi") +
  ylab("High carotene") +
  labs(title = "High carotene (=1) or Not high carotene (=0) vs bmi",
       caption = "red = fitted line, blue dashed = moving average")



carotene_pred <- cbind(
  carotene_pred,
  logit = predict(model_bmi_glm, se.fit = TRUE))
glimpse(carotene_pred)
#carotene_pred |> mutate(logit.residual.scale = NULL) -> carotene_pred

lambda <- qnorm(1 - 0.05/2)
carotene_pred |> mutate(
  logit.lwr = logit.fit - lambda*logit.se.fit,
  logit.upr = logit.fit + lambda*logit.se.fit) -> carotene_pred
glimpse(carotene_pred)

carotene_pred |> mutate(
  odds.lwr = exp(logit.lwr),
  odds.upr = exp(logit.upr)) -> carotene_pred
glimpse(carotene_pred)

carotene_pred |> mutate(
  p.lwr = odds.lwr/(1 + odds.lwr),
  p.upr = odds.upr/(1 + odds.upr)) -> carotene_pred
glimpse(carotene_pred)

ggplot(carotene_pred, aes(bmi, lowplasma_01)) +
  geom_point() +
  geom_line(aes(y = phat), color = "red", size = 1) +
  geom_ribbon(aes(ymin = p.lwr, ymax = p.upr), alpha = 0.2) +
  geom_smooth(se = FALSE, linetype = "dashed",span=2) +
  xlab("BMI") +
  ylab("High carotene") +
  labs(title = "High carotene (=0) or Low carotene (=1) vs bmi",
       caption = "red = fitted line, with 95% confidence interval, blue dashed = moving average")
ggsave("../Figures/2/mean_pred_log.png", width = 8, height = 6, dpi = 600)

model_bmi_sum <- summary(model_bmi_glm)
model_bmi_sum$coefficients

null_model <- glm(lowplasma_01 ~ 1, family = binomial, data = carotene)

anova(null_model,model_bmi_glm, test = "LRT")


model_bmi_infl <- influence(model_bmi_glm)
glimpse(model_bmi_infl)

carotene_pred <- cbind(carotene,
                   xbeta = predict(model_bmi_glm),
                   v = model_bmi_infl$hat)
glimpse(carotene_pred)

model_bmi_lin <- lm(lowplasma_01 ~ bmi, data = carotene)

pplus1_carotene <- length(model_bmi_glm$coefficients)
n <- nobs(model_bmi_glm)
pplus1_lin <- length(coef(model_bmi_lin))
pplus1_glm <- length(coef(model_bmi_glm))
n <- nobs(model_bmi_lin)  # same n for both models


pplus1 <- length(model_bmi_lin$coefficients)
n <- nobs(model_bmi_lin)

carotene_pred_lin <- mutate(carotene,
                    yhat = predict(model_bmi_lin),
                    r = rstudent(model_bmi_lin),
                    v = hatvalues(model_bmi_lin),
                    D = cooks.distance(model_bmi_lin))

carotene_leverage <- carotene |> 
  mutate(
    leverage_lin = hatvalues(model_bmi_lin),
    leverage_glm = hatvalues(model_bmi_glm),
    xbeta_glm = predict(model_bmi_glm, type = "link")
  )


carotene_leverage_long <- carotene_leverage |> 
  pivot_longer(cols = c(leverage_lin, leverage_glm), 
               names_to = "model", values_to = "leverage") |> 
  mutate(model = recode(model,
                        leverage_lin = "Linear regression",
                        leverage_glm = "Logistic regression"))

ggplot(carotene_pred, aes(x = xbeta, y = v)) +
  geom_point() +
  geom_hline(yintercept = c(2*pplus1_carotene/n)) +
  facet_wrap(~ lowplasma_01)

ggplot(carotene_leverage_long, aes(x = bmi, y = leverage, color = model)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 2*pplus1_lin/n, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 2*pplus1_glm/n, color = "red", linetype = "dashed") +
  labs(x = "BMI",
       y = "Leverage",
       title = "Comparison of Leverage: Linear vs Logistic Regression",
       caption = "Dashed line represent leverage threshold",
       color = "Model Type") +
  theme_minimal()
ggsave("../Figures/2/leverage.png", width = 8, height = 6, dpi = 600)
carotene_pred |> mutate(devresid = model_bmi_infl$dev.res,
                    stddevresid = devresid/sqrt(1 - v)) -> carotene_pred
glimpse(carotene_pred)
ggplot(carotene_pred, aes(sample = stddevresid)) +
  geom_qq() + geom_qq_line() +
  labs(title = "QQ-plot Standardized deviance residuals")
ggsave("../Figures/2/standardized_deviance.png", width = 8, height = 6, dpi = 600)
ggplot(carotene_pred, aes(x = bmi, 
                      y = stddevresid, 
                      color = as.factor(lowplasma_hl))) +
  geom_point() +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"),
             linewidth = 1) +
  labs(title = "Standardised deviance residuals vs linear predictor",
       x = "BMI", y = "devstd",
       color = "Y")
ggsave("../Figures/2/standardized_deviance_bmi.png", width = 8, height = 6, dpi = 600)
carotene_pred <- carotene_pred |> mutate(
  large_resid = abs(stddevresid) > 2
)
carotene_pred <- mutate(carotene_pred, 
                    Dcook = cooks.distance(model_bmi_glm))


ggplot(carotene_pred, aes(x = bmi, y = Dcook, color = factor(lowplasma_01))) +
  geom_point() +
  geom_point(data = filter(carotene_pred, abs(stddevresid) > 2), aes(color = "|devstd| > 2"),
             size = 3) +
  geom_hline(yintercept = 4/n, linewidth = 1) +
   #facet_wrap(~ lowplasma_01, labeller = label_both) +
  labs(title = "Cook's distance vs linear predictor by temp diff",
       x = "BMI", 
       color = "Highlight") +
  theme(legend.position = "top")
ggplot(carotene_pred, aes(x = bmi, y = Dcook, color = factor(lowplasma_01))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_point(data = filter(carotene_pred, large_resid), 
             shape = 21, size = 4, stroke = 1.2, color = "black") +
  geom_hline(yintercept = 4 / nrow(carotene_pred), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("blue", "orange"), 
                     labels = c("High β-carotene", "Low β-carotene")) +
  labs(
    title = "Cook's Distance vs BMI",
    x = "BMI",
    y = "Cook's Distance",
    color = "Plasma β-carotene level",
    caption = "Red dashed line: Cook's distance threshold (4/n)\nBlack circles: Large residual observations"
  ) +
  theme_minimal()
ggsave("../Figures/2/cook.png", width = 8, height = 6, dpi = 600)
