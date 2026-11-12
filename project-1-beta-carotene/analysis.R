library(readxl)
library(tidyverse); theme_set(theme_bw() + theme(text = element_text(size = 18)))
library(car)
library(rstatix)
library(GGally)
library(patchwork)


carotene <- readxl::read_xlsx("data/carotene.xlsx")
head(carotene)


#1 BMI

ggplot(carotene, aes(y = betaplasma, x=bmi)) + geom_point()

#making the Lin and Log models
bmi_lm = lm(data = carotene , betaplasma ~ bmi)       
bmi_log_lm = lm(data = carotene, log(betaplasma) ~ bmi)

#adding the bmi residuals to the main data
carotene <- mutate(carotene, e_bmi = bmi_lm$residuals, e_bmi_log = bmi_log_lm$residuals)
bmi_elims <- max(abs(carotene$e_bmi))*c(-1, 1)
bmi_log_elims <- max(abs(carotene$e_bmi_log))*c(-1, 1)

#1(a)
#Linear model residuals vs bmi
ggplot(data = carotene, aes(x = bmi, y = e_bmi)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  expand_limits(y = bmi_elims) +
  xlab("BMI") +
  ylab("Residual") +
  labs(tag = "A") +
  labs(title = "Linear Model : Residuals vs BMI")

#Log model residuals vs bmi 
ggplot(data = carotene, aes(x = bmi, y = e_bmi_log )) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  expand_limits(y = bmi_log_elims) +
  xlab("BMI") +
  ylab("Residual") +
  labs(tag = "A") +
  labs(title = "Logarithmic Model : Residuals vs x-values")

#The Logarithmic Model is more suitable (check asymetry of residuals of linear model, and severity of some of them)


#1(b)

#table with b estimates and 95% C.I.
cbind(beta =bmi_log_lm$coefficients, confint(bmi_log_lm)) |> round(digits=3)

#plotting

c(min(carotene$bmi), max(carotene$bmi))
bmi_seq = data.frame(bmi = seq(15,55))

bmi_seq |> mutate(
  fit = predict(bmi_log_lm, newdata = bmi_seq),
  conf = predict(bmi_log_lm, newdata = bmi_seq, interval = "confidence"),
  pred = predict(bmi_log_lm, newdata = bmi_seq, interval = "prediction")) -> 
  bmi_log_ints

#plot

#Linear plot in log(y) scale
ggplot(bmi_log_ints, aes(x = bmi)) + 
  geom_point(data = carotene, aes(y = log(betaplasma))) + 
  geom_line(aes(y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = conf[, "lwr"], ymax = conf[, "upr"]), alpha = 0.2) +
  geom_line(aes(y = pred[, "lwr"]),
            color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = pred[, "upr"]),
            color = "red", linetype = "dashed", linewidth = 1) +
  xlab("BMI") +
  ylab("ln(betaplasma (ng/ml) )") +
  labs(title = "Plasma β-carotene levels : log-betaplasma by BMIl",
       caption = "fitted log-lin model, 95% conf. and pred. intervals")

#Log plot with y scale
ggplot(bmi_log_ints, aes(x = bmi)) + 
  geom_point(data = carotene, aes(y = betaplasma)) + 
  geom_line(aes(y = exp(fit)), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = exp(conf[, "lwr"]), ymax = exp(conf[, "upr"])), alpha = 0.2) +
  geom_line(aes(y = exp(pred[, "lwr"])),
            color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = exp(pred[, "upr"])),
            color = "red", linetype = "dashed", linewidth = 1) +
  xlab("BMI") +
  ylab("betaplasma (ng/ml)") +
  labs(title = "Plasma β-carotene levels : betaplasma by BMI",
       caption = "fitted log-lin model, 95% conf. and pred. intervals")


paste0("betaplasma = ",round(exp(bmi_log_lm$coefficients[1]),2),"*exp(",round(bmi_log_lm$coefficients[2],2),"x). Where x is the BMI level (kg/m2)")

#1(c)
#Expected betaplasma changes

#One-unit BMI increase
paste0("The betaplasma levels decay exponentially with respect to BMI, with a decay factor of ",round(exp(bmi_log_lm$coefficients[2]),3), 
       " per unit change in BMI (95% CI:[",round(exp(confint(bmi_log_lm)[2,][1]),3)," - ",round(exp(confint(bmi_log_lm)[2,][2]),3) ,"]), while holding other factors constant, corresponding to an approximate of ",100*(1-round(exp(bmi_log_lm$coefficients[2]),3)),"% decrease per unit change.")

#One-unit BMI decrease
paste0("Holding other variables constant, a one-unit decrease in BMI is associated with a ",round((1/exp(bmi_log_lm$coefficients[2])-1)*100,2),"% increase in betaplasma on average (multiplicative factor :",round(1/exp(bmi_log_lm$coefficients[2]),3) ," ,95% CI:[",round(1/exp(confint(bmi_log_lm)[2,][2]),3)," - ",round(1/exp(confint(bmi_log_lm)[2,][1]),3) ,"]), (the inverse of the ",100*(1-round(exp(bmi_log_lm$coefficients[2]),3)), "% decay).")

#Ten-unit BMI decrease
paste0("Holding other variables constant, a ten-unit decrease in BMI is associated with a ",round((exp(bmi_log_lm$coefficients[2])^(-10)-1)*100,2),"% increase in betaplasma on average (multiplicative factor :",round(exp(bmi_log_lm$coefficients[2])^(-10),2)," ,95% CI:[", round(exp(confint(bmi_log_lm)[2,][2])^(-10), 3)," - ",round(exp(confint(bmi_log_lm)[2,][1])^(-10), 3),"])") 


#1(d)

bmi_sum <- summary(bmi_log_lm)

bmi_sum$coefficients
bmi_sum$coefficients["bmi","Pr(>|t|)"] #p value

paste0("We reject the Null Hypothesis at a 95% confidence level and can claim that BMI does have an effect on Betaplasma levels")

#2 Smoking Habits

#2(a)

#Releveling with "Non-Smokers" as the reference and turning smokstat to "factor" type

carotene$smokstat
carotene$smokstat = factor(carotene$smokstat, levels = c(1,2,3), labels = c("Never","Former","Current Smoker"))
carotene$smokstat = relevel(carotene$smokstat, ref = "Never")

#betaplasma mean and sd per smoking status
carotene %>%
  group_by(smokstat) %>%
  summarise(mean(betaplasma), sd(betaplasma))

#log(betaplasma) mean and sd per smoking status
carotene %>%
  group_by(smokstat) %>%
  summarise(mean(log(betaplasma)), sd(log(betaplasma)))

#
#lin plasma
ggplot(carotene, aes(x = smokstat, y = betaplasma)) +
  geom_point(aes(color = smokstat), size = 3) +
  expand_limits(y = c(0,1500)) +
  ylab("Betaplasma (ng/ml)") +
  xlab("Smoking Status") +
  labs(title = "Betaplasma vs Smoking status")


#boxplot and violin-plot for betaplasma
ggplot(carotene, aes(x = smokstat, y = betaplasma, color = smokstat)) +
  geom_violin() +
  stat_summary(fun = mean, geom="point", size = 2, color = "red") +
  geom_boxplot(width = 0.1) +
  labs(title = "Boxplot and violin-plot for betaplasma vs Smoking status") +
  theme(legend.position = "bottom")

#log plasma
ggplot(carotene, aes(x = smokstat, y = log(betaplasma))) +
  geom_point(aes(color = smokstat), size = 3) +
  expand_limits(y = c(2,7.5)) +
  ylab("Betaplasma (ng/ml)") +
  xlab("Smoking Status") +
  labs(title = "Log-Betaplasma vs Smoking status")

#boxplot and violin-plot for log(betaplasma) per smoking status
ggplot(carotene, aes(x = smokstat, y = log(betaplasma), color = smokstat)) +
  geom_violin() +
  stat_summary(fun = mean, geom="point", size = 2, color = "red") +
  geom_boxplot(width = 0.1) +
  labs(title = "Boxplot and violin-plot for Log-betaplasma vs Smoking status") +
  theme(legend.position = "bottom")

#2(b)

#Creating versions of log-betaplasma model with different smoking status references

#"No Smoker" as reference
carotene$smokstat = relevel(carotene$smokstat, ref = "Never")
Smok_n_lm <- lm(data = carotene, log(betaplasma)~smokstat)

#"Current Smoker" as reference
carotene$smokstat = relevel(carotene$smokstat, ref = "Current Smoker")
Smok_s_lm <- lm(data = carotene, log(betaplasma)~smokstat)

paste0("In the Log-Lin model with non-smokers as the reference category, the dummy variable x1i equals to 1 for Former smokers and 0 otherwise, and is connected to parameter b1(smokstatFormer).The dummy variable x2i for Current Smokers works in a similar way and is connected to b2(smokstatCurrent Smoker. The Non-Smoker category is captured by the intercept term when both dummy variables equal 0")
paste0("Likewise for the log-lin model with Current Smoker as the reference.")

Smok_n_sum <- summary(Smok_n_lm)
Smok_s_sum <- summary(Smok_s_lm)

Smok_n_sum$coefficients[,1:2]
Smok_s_sum$coefficients[,1:2]

paste0("The standard error of b0 in model 1 is 0.059 and corresponds to the standard error of log(Betaplasma) for non-smokers. The 2 dummy variables are compared to the non-smoker status, so their standard errors have a minimum of the standard error of b0. Likewise for model 2, and since the standard error of log(betaplasma) for Non-smokers is smaller than the s.e. for Current smokers, the standard errors in the first model are smaller. ")

Model2b <- Smok_n_lm
Model2b_sum <- Smok_n_sum

#2(c)


#calculating the expected log-plasma b-carotene with 95% C.I. VS Smoking Status

smst = data.frame(smokstat = c("Never", "Former", "Current Smoker"))


predict(Smok_n_lm, smst, interval = "confidence")
predict(Smok_s_lm, smst, interval = "confidence")

exp(predict(Smok_n_lm, smst, interval = "confidence"))
exp(predict(Smok_s_lm, smst, interval = "confidence"))

paste0("The expected values and C.I.s are the same because both models are the same model, with different parameterizations")

#2(d)
anova(Model2b)
Model2b_sum

paste0("We perform an F-test to see whether there are any significant changes in log plasma b-carotene between any of the smokstat categories. H0: b1=b2=0 vs H1: either b1 or b2 =/= 0. The observed value of the F-test is 5.75.If the Null Hypothesis is true F ~ F(2, 312) , p-value = 0.00353. Therefore we can reject the Null Hypothesis at a level of significance a=5% and conclude that smoking status is relevant to the levels of betaplasma")


#from part 2: 
head(carotene)

## PART 3 ##
# 3(a) #
mutate(carotene,
       sex = factor(sex,
                     levels = c(1, 2),
                     labels = c("Male", "Female"))) -> carotene
carotene$sex = relevel(carotene$sex, ref = "Female")
mutate(carotene,
       vituse = factor(vituse,
                    levels = c(1, 2, 3),
                    labels = c("Yes, fairly often", "Yes, not often", "No"))) -> carotene

carotene$vituse = relevel(carotene$vituse, ref = "No")
glimpse(carotene)
head(carotene)

# 3(b) #

#Calculate pairwise correlations between the continuous x-variables bmi, age, calories, fat, choles
# terol, fiber, alcohol, and betadiet, and present, and plot, the pairs where the correlations are
# stronger than 06. 


cortab<-carotene |> select(bmi, age, calories, fat, cholesterol, fiber, alcohol,betadiet) |> 
  #mutate(x1x2 = bmi*age) |> 
  cor_test() |> filter(var1 < var2)

cortab_high<-cortab[abs(cortab$cor) > 0.6, ]
cortab_high
#cortab_low<-cortab[cortab$cor < -0.6, ] #to double check; no negative corr that is significant
#cortab_low

#we have three-way high correlation numbers between calories, fat and cholesterol 

carotene |> select(calories, fat, cholesterol) |>
  #mutate(x1x2 = temp*pressure) |> 
  ggpairs()


# Plots

ggplot(carotene, aes(x = alcohol, y = calories)) + geom_point() +
  labs(x = "alcohol",
       y = "calories")

ggplot(carotene, aes(x = alcohol, y = fat)) + geom_point() +
  labs(x = "alcohol",
       y = "fat")

ggplot(carotene, aes(x = alcohol, y = cholesterol)) + geom_point() +
  labs(x = "alcohol",
       y = "cholersterol")


ggplot(carotene, aes(x = alcohol, y = betaplasma)) + geom_point() +
  labs(x = "alcohol",
       y = "betaplasma")

# comment: there is an extreme value!

# 3(c) #

mod_all<-lm(log(betaplasma)~  bmi + age + calories + fat + cholesterol + 
              fiber + alcohol + betadiet + smokstat + sex + vituse, data=carotene)
mod_all_sum<-summary(mod_all)
mod_all_sum

# VIF/GVIF-values 

vif(mod_all, type = "predictor")
vif(mod_all)


#from lectures we know VIF=1/(1-R^2), where R^2 is prop that can be explained by other x-vars.
#so more than 80% of variability is explained by vars that have VIF > 5. Calories (!), fat most problematic.


# Refit the model without the most problematic x-variable (Model.3(c)). Check new VIF/GVIF-values
head(carotene)
mod_ref<-lm(log(betaplasma)~  bmi + age + fat + cholesterol + 
              fiber + alcohol + betadiet + smokstat + sex + vituse, data=carotene)
mod_ref_sum<-summary(mod_ref)
mod_ref_sum
length(mod_ref$coefficients) # 13 coefficients

#the GVIF vals are much closer to one, none are bigger than 1.5 for GVIF^(1/(2*Df))


# 3(d) #

# Use Model 3(c) and state how you define the different variables and which-parameter they each
# are connected to

# betaplasma~e^sum(beta_i);
# beta estimates
mod_ref_sum #gives beta_i as in log(betaplasma)~sum(beta)
CIs<- cbind(expbeta = mod_ref$coefficients,
               confint(mod_ref)) |> round(digits = 4)
CIs

# e^beta estimates, CIS for e^beta:
expCIs<- cbind(expbeta = exp(mod_ref$coefficients),
      exp(confint(mod_ref))) |> round(digits = 4)

expCIs


# Tests

# (i) H0: beta_{bmi}=0 (t-test)

mod_ref_sum$coefficients["bmi",]
# the p-value is small, therefore can reject H0; bmi does have an influence; 
#statistic, beta/(d(beta)), is t(df)-disttributed when H0 is true; the df given by: 
df<-mod_ref$df.residual 

# (ii) H0: beta_{age}= ...[all but beta_{bmi}] ...= beta_{vituse 'yes..'}=0 (partial F-test)

mod1b = lm(data = carotene, log(betaplasma) ~ bmi)
anova(mod1b, mod_ref)

# p-value is small, there are significant differences between the models (we reject the H0 that all non-bmi betas are 0);
# the statistic, (Q/k)/s^2 is F(k, n-(p+1))-distributed when H0 is true, and:
# Q=SS(Err_reduced)-SS(Err_full);
# n-(p+1) is same as df above=302;
# k is 11, number of betas in 'H0:[all but beta_{bmi}]=0' (beware the categorical vars that give >1 beta)
# *ref: lecture 4 slide 31/34


# (iii) H0: beata_{age}= ...[all but betas_{smokstat}] ...= beta_{vituse 'yes..'}=0 (partial F-test)

# same approach, now we are mindful that smokstat is cat. var. with two betas in the model, so k=10.

mod2b <- lm(data = carotene, log(betaplasma)~smokstat)
anova(mod2b, mod_ref)

# result is largely the same as above in (ii).


# 3(e) #

#outliers, non-constant variance, and non-normality

bp_pred <- mutate(carotene, 
                   yhat = predict(mod_ref),
                   r = rstudent(mod_ref))


highlightcolors <- c("|r*|>3" = "red")
ggplot(bp_pred, aes(x = yhat, y = r)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_point(data = filter(bp_pred, abs(r) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Studentized residuals vs predictor",
       subtitle = "log-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")

#check for constant var: 
ggplot(bp_pred, aes(x = yhat, y = sqrt(abs(r)))) +
  geom_point() +
  geom_hline(yintercept = c(0, sqrt(qnorm(0.75)), sqrt(2))) +
  geom_hline(yintercept = sqrt(3), linetype = 2) +
  geom_point(data = filter(bp_pred, abs(r) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Sqrt absolute studentized residuals vs predictor",
       subtitle = "Log-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")


# there is no noticeable systematic non-linear pattern. res variance could be considered constant. 

# 3 (f) #

# calculate the leverage for Model.3(c) and plot them against the linear predictor, with suitable reference
# lines.
################################
car_pred<-bp_pred
Model3<-mod_ref
Model3_sum<-summary(Model3)

pplus1 <- length(Model3$coefficients)
n <- nobs(Model3)

################################
car_pred <- mutate(carotene,
                    yhat = predict(Model3),
                    r = rstudent(Model3),
                    v = hatvalues(Model3),
                    D = cooks.distance(Model3))
hist(carotene$betaplasma)
unique(sort(carotene$betaplasma))

highlightcolors <- c("|r*|>3" = "red",
                     "betaplasma>1250" = "magenta",
                     "alcohol>200"= "cyan", 
                     "cholesterol>900" = "pink",
                     "all data" = "orange",
                     "excl.length>200" = "blue")

# leverage vs fitted values: 
leverage_plot <- ggplot(car_pred, aes(x = yhat, y = v)) +
  geom_point(size = 2) +
  geom_point(data = filter(car_pred, betaplasma > 1250),
             aes(color = "betaplasma>1250"), size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200),
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900),
             aes(color = "cholesterol>900"), size = 3) +
  geom_hline(yintercept = pplus1/n) +
  geom_hline(yintercept = 2*pplus1/n, color = "red") +
  labs(
    title = "Leverage vs fitted values",
    x = "Fitted log(beta-carotene)",
    y = "Leverage",
    caption = "Average leverage (black); twice average leverage (red)",
    color = "Highlight"
  ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)

print(leverage_plot)

#clearly visible that the alcohol outlier has very high leverage; cholesterol is a bit high too
#see: 
  
ggplot(carotene, aes(x = alcohol, y = betaplasma)) + geom_point() +
  labs(x = "alcohol",
       y = "betaplasma")
ggplot(carotene, aes(x = cholesterol, y = betaplasma)) + geom_point() +
  labs(x = "cholesterol",
       y = "betaplasma")

################################

#Looking in more detail the value with the highest leverage
summary(car_pred)
#View(car_pred[car_pred$v == max(car_pred$v),])

paste0("We can see the observation with the highest leverage corresponds to the person with the highest alcohol consumption per week (203 drinks per week) and the highest consumption of calories per day (27 MJ) ")

#updating highlightcolors
highlightcolors = c(highlightcolors,
                    "v>0.8" = "magenta")

#Plotting alcohol vs betaplasma, focusing on the value with the highest leverage
ggplot(car_pred, aes(y = log(betaplasma), x = alcohol)) +
  geom_point() +
  geom_point(data = filter(car_pred, v>0.8),
             aes(color = "v>0.8"), size = 3) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors) +
  labs(title = "Log-betaplasma vs alcohol.")

#Plotting calories vs betaplasma, focusing on the value with the highest leverage
ggplot(car_pred, aes(y = log(betaplasma), x = calories)) +
  geom_point() +
  geom_point(data = filter(car_pred, v>0.8),
             aes(color = "v>0.8"), size = 3) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)  +
  labs(title = "Log-betaplasma vs calories.")

# The observation with extreme alcohol consumption has very high leverage,
# indicating an unusual combination of predictor values.
# However, high leverage does not necessarily imply high influence.
# We investigate influence using Cook's distance and DFBETAS below.


#3(g)

n = length(carotene$smokstat)
pplus1 = length(coef(Model3))

#Plotting Cook's D , vs Yhat
f1.carotene = pplus1
f2.carotene = Model3$df.residual
cook.limit.carotene = qf(0.5, f1.carotene, f2.carotene)


n = length(carotene$smokstat)
pplus1 = length(coef(Model3))
k1 = length(coef(bmi_log_lm))
k2 = length(coef(Model2b))


highlightcolors <- c("|r*|>3" = "red",
                     "betaplasma>1250" = "magenta",
                     "alcohol>200"= "cyan", 
                     "cholesterol>900" = "pink",
                     "all data" = "orange",
                     "excl.alc>200" = "blue",
                     "excl.chol>900"="green")

ggplot(car_pred, aes(yhat, D)) + 
  geom_point(size = 3) +
  geom_point(data = filter(car_pred, betaplasma > 1250), 
             aes(color = "betaplasma>1250"), size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900), 
             aes(color = "cholesterol>900"), size = 3) +
  geom_hline(yintercept = cook.limit.carotene, color = "red") +
  geom_hline(yintercept = 4/n, linetype = 2, color = "red") +
  xlab("Fitted values") +
  ylab("D_i") +
  labs(title = "Carotene: Cook's D",
       caption = "4/n (dashed), F_0.5, p+1, n-(p+1) (solid)",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors)


paste0("The observation with the largest Cook's distance, is NOT the observation with the highest leverage")
# it is the cholesterol!the alcohol >200 has limited inluence. see: 
ggplot(car_pred,aes (v, D))+
        geom_point(size = 3) +
        geom_point(data = filter(car_pred, betaplasma > 1250), 
                   aes(color = "betaplasma>1250"), size = 3) +
        geom_point(data = filter(car_pred, alcohol > 200), 
                   aes(color = "alcohol>200"), size = 3) +
        geom_point(data = filter(car_pred, cholesterol > 900), 
                   aes(color = "cholesterol>900"), size = 3) +
        xlab("leverage") +
        ylab("D_i") +
        labs(title = "Carotene: leverage vs Cook's D",
             color = "Highlight") +
        scale_color_manual(values = highlightcolors)

# was cholesterol >900 obs the one with high var? 
ggplot(car_pred, aes(x = yhat, y = r)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_point(data = filter(bp_pred, abs(r) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Studentized residuals vs predictor",
       subtitle = "log-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")
# 
ggplot(car_pred, aes(x = yhat, y = r)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_point(data = filter(car_pred, betaplasma > 1250), 
             aes(color = "betaplasma>1250"), size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900), 
             aes(color = "cholesterol>900"), size = 3) +
  labs(title = "Studentized residuals vs predictor",
       subtitle = "log-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")

# yes it was! 
hist(carotene$cholesterol)



head(dfbetas((Model3)))

car_pred <- mutate(
  car_pred,
  df0 = dfbetas(Model3)[, "(Intercept)"],
  df_fat = dfbetas(Model3)[, "fat"],
  df_cholesterol = dfbetas(Model3)[, "cholesterol"],
  df_alcohol = dfbetas(Model3)[, "alcohol"])

## investigating the cholesterol>900 measure 8the one with large cooks D_i) vs.
# alcohol>200 (the one with large leverage) ##
#

f1.carotene = pplus1
f2.carotene = Model3$df.residual
cook.limit.carotene = qf(0.5, f1.carotene, f2.carotene)

highlightcolors <- c("|r*|>3" = "red",
                     "betaplasma>1250" = "magenta",
                     "alcohol>200"= "cyan", 
                     "cholesterol>900" = "pink",
                     "all data" = "orange",
                     "excl.alc>200" = "blue",
                     "excl.chol>900"="green")

highlightshapes <- c("v>0.8" = 22, # square 
                     "Cook's D>0.09"= 24) # triangle ; in lecture ex it is 0.1, but we have 0.09 as largest: 
unique(sort(car_pred$D)) # + see carotene: leverage vs cooks d plot above
filter(car_pred, D > 0.09)
filter(car_pred, v > 0.8)

library(patchwork)
df0<-ggplot(car_pred, aes(x = yhat, y = df0)) +
  geom_point(size = 2) +
  geom_point(data = filter(car_pred, abs(r) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(car_pred, D > 0.09),
             aes(shape = "Cook's D>0.09"),
             size = 3) +
  geom_point(data = filter(car_pred, v > 0.8),
             aes(shape = "v>0.8"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.carotene)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_0(i)") +
  xlab("Fitted values") +
  theme(legend.position = "none") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)


dffat<-ggplot(car_pred, aes(x = yhat, y = df_fat)) +
  geom_point(size = 2) +
  geom_point(data = filter(car_pred, abs(r) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(car_pred, D > 0.09),
             aes(shape = "Cook's D>0.09"),
             size = 3) +
  geom_point(data = filter(car_pred, v > 0.8),
             aes(shape = "v>0.8"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.carotene)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_fat(i)") +
  xlab("Fitted values") +
  theme(legend.position = "none") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

dfchol<-ggplot(car_pred, aes(x = yhat, y = df_cholesterol)) +
  geom_point(size = 2) +
  geom_point(data = filter(car_pred, abs(r) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(car_pred, D > 0.09),
             aes(shape = "Cook's D>0.09"),
             size = 3) +
  geom_point(data = filter(car_pred, v > 0.8),
             aes(shape = "v>0.8"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.carotene)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_chol(i)") +
  xlab("Fitted values") +
  theme(legend.position = "none") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

dfalc<-ggplot(car_pred, aes(x = yhat, y = df_alcohol)) +
  geom_point(size = 2) +
  geom_point(data = filter(car_pred, abs(r) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(car_pred, D > 0.09),
             aes(shape = "Cook's D>0.09"),
             size = 3) +
  geom_point(data = filter(car_pred, v > 0.8),
             aes(shape = "v>0.8"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.carotene)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_alc(i)") +
  xlab("Fitted values") +
  theme(legend.position = "none") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)


dfbetas_plot <- (df0 | dffat) / (dfchol | dfalc) +  plot_annotation(
  title = "Carotene: DFBETAS Impacts",
  subtitle = "full model",
  caption = "y = sqrt(F_0.5) and 2/sqrt(n)"
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


dfbetas_plot

# seems that the observation with the highest leverage is not the one that influences obs most
# instead, it is the one with the highest D_i value. Let us try to investigate with 
# more plots. 

# this is likely explained by the fact that alcohol cov is not significant,
# while cholsterol is on alpha=0.1. 

# also: Note cook's upper limit is not that large (interesting)

#might still be interesting to fit models with and without alc and chol outliers

summary(Model3)

car_excl_alc <- update(Model3, 
                   data = filter(carotene, alcohol<=200))

car_excl_chol <- update(Model3, 
                       data = filter(carotene, cholesterol<=900))
car_pred <- mutate(car_pred,
                    yhat_excl_alc = predict(car_excl_alc,
                                        newdata = carotene),
                   yhat_excl_chol = predict(car_excl_chol,
                                           newdata = carotene))

#summary()
# --
# betaplasma vs alc
#log(weight)=log(betaplasma)
#log(length)=alcohol
plot(car_pred$yhat_excl_alc, car_pred$yhat_excl_chol)

# hinting at why it was the observation with the largest cooks distance that mattered
names(car_pred)
p1<-ggplot(car_pred, aes(alcohol, yhat))+
  geom_point(size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900), 
             aes(color = "cholesterol>900"), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +  #  simple reg line
  
  xlab("alcohol") +
  ylab("yhat") +
  theme(legend.position = "none") +
  scale_color_manual(values = highlightcolors)

p2<-ggplot(car_pred, aes(cholesterol, yhat))+
  geom_point(size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900), 
             aes(color = "cholesterol>900"), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +  #  simple reg line
  theme(legend.position = "none") +
  xlab("cholesterol") +
  ylab("yhat") +
  scale_color_manual(values = highlightcolors)


trend_plot<-(p1| p2) + plot_annotation(
  title = "fitted values vs covariates",
  subtitle = "full model"
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

trend_plot

# plot without the problematic ones 

names(car_pred)
p3<-ggplot(car_pred, aes(alcohol, yhat_excl_alc))+
  geom_point(size = 3) +
   geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900), 
             aes(color = "cholesterol>900"), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +  #  simple reg line
  
  xlab("alcohol") +
  ylab("yhat_excl_alc") +
  theme(legend.position = "none") +
  scale_color_manual(values = highlightcolors)

p4<-ggplot(car_pred, aes(cholesterol, yhat_excl_chol))+
  geom_point(size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
   geom_point(data = filter(car_pred, cholesterol > 900), 
              aes(color = "cholesterol>900"), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +  #  simple reg line
  theme(legend.position = "none") +
  xlab("cholesterol") +
  ylab("yhat_excl_chol") +
  scale_color_manual(values = highlightcolors)


trend_plot2<-(p3| p4) + plot_annotation(
  title = "fitted values vs covariates",
  subtitle = "Models refitted after excluding the respective influential observations"
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

trend_plot2



# part 4a:reduced + cooks #
# 4(a)
# Calculate and plot Cook’s distance against the linear predictor, with suitable reference lines, for
# Model.3(c) on the reduced data set 

car_excl_lm <- update(Model3, 
                      data = filter(carotene, cholesterol<=900))

car_pred <- mutate(car_pred,
                   yhat_excl = predict(car_excl_lm,
                                       newdata = carotene))
#plot(car_pred$yhat_excl_chol, car_pred$yhat_excl)

f1.excl <- length(car_excl_lm$coefficients)
f2.excl <- car_excl_lm$df.residual
cook.limit.excl <- qf(0.5, f1.excl, f2.excl)

n_excl <- nobs(car_excl_lm)
car_excl_pred <- cbind(filter(carotene, cholesterol <= 900),
                       yhat = predict(car_excl_lm),
                       r = rstudent(car_excl_lm),
                       D = cooks.distance(car_excl_lm))

unique(sort(car_excl_pred$cholesterol))
max(car_excl_pred$D)

highlightcolors <- c("|r*|>3" = "red",
                     "betaplasma>1250" = "magenta",
                     "alcohol>200"= "cyan", 
                     "cholesterol>900" = "pink",
                     "all data" = "orange",
                     "excl.alc>200" = "blue",
                     "excl.chol>900"="green",
                     "cholesterol>800" = "brown"
)

pfull<- ggplot(car_pred, aes(yhat, D)) + 
  geom_point(size = 3) +
  geom_point(data = filter(car_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_pred, cholesterol > 900), 
             aes(color = "cholesterol>900"), size = 3) +
  geom_hline(yintercept = cook.limit.carotene, color = "red") +
  geom_hline(yintercept = 4/n, linetype = 2, color = "red") +
  xlab("Fitted values") +
  ylab("D_i") +
  theme(legend.position= "none")+
  scale_color_manual(values = highlightcolors)

pexcl<-ggplot(car_excl_pred, aes(yhat, D)) + 
  geom_point(size = 3) +
  geom_point(data = filter(car_excl_pred, alcohol > 200), 
             aes(color = "alcohol>200"), size = 3) +
  geom_point(data = filter(car_excl_pred, cholesterol > 800), 
             aes(color = "cholesterol>800"), size = 3) +
  geom_hline(yintercept = cook.limit.carotene, color = "red") +
  geom_hline(yintercept = 4/n, linetype = 2, color = "red") +
  xlab("Fitted values") +
  ylab("D_i") +
  theme(legend.position ="none") +
  scale_color_manual(values = highlightcolors)

plt_cooks<-(pfull| pexcl) + plot_annotation(
  title = "fitted values vs cooks D",
  subtitle = "full (left) vs. excl (rhs) model"
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

plt_cooks


# the highest cooks value is no longer the cholesterol outlier (the next in line that is)


#Calculating VIFS in both full and reduced models
#
vifreduced = vif(car_excl_lm, type = 'predictor')
viffull = vif(Model3, type = 'predictor')

vifreduced
viffull
# no interactions 


# 4a:Calculating VIFS in both full and reduced models
vifreduced = vif(car_excl_lm)
viffull = vif(Model3)

#Calculating the diffence of GVIF values between models
vifdiff = data.frame( viffull, GVIFdiff = vifreduced[,1]-viffull[,1], "GVIF^(1/(2*Df)DIFF" = vifreduced[,3]-viffull[,3])

#Ordering the absolute differences of VIF in descending order
vifdiff[order(abs(vifdiff$GVIFdiff), decreasing = T),]

# 4a:correlations #
# replace carotene by filter(carotene, cholesterol<=900)
cortab<-filter(carotene, cholesterol<=900) |> select(bmi, age, calories, fat, cholesterol, fiber, alcohol,betadiet) |> 
  #mutate(x1x2 = bmi*age) |> 
  cor_test() |> filter(var1 < var2)

cortab_high<-cortab[abs(cortab$cor) > 0.6, ]
cortab_high
#cortab_low<-cortab[cortab$cor < -0.6, ] #to double check; no negative corr that is significant
#cortab_low
#same
carotene |> select(calories, fat, cholesterol) |>
  #mutate(x1x2 = temp*pressure) |> 
  ggpairs()
### 
# 4(b)

#excluding the outlier from data

data_excl=filter(carotene, cholesterol<= 900)
# calories + bmi + age + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse
ex_null_lm<-lm(log(betaplasma) ~ 1, data = data_excl)
ex_mod_bmi<-lm(log(betaplasma)~  bmi, data=data_excl)
ex_mod4bmi_fiber<-lm(log(betaplasma)~  bmi+fiber, data=data_excl)
ex_mod4bmi_fiber_calories<-lm(log(betaplasma)~  bmi+fiber+ calories, data=data_excl)
ex_mod4bmi_fiber_calories_vituse<-lm(log(betaplasma)~  bmi+fiber+ calories + vituse, data=data_excl)
BIC(ex_mod4bmi_fiber_calories_vituse)

# null -> + bmi -> +fiber > +calories > +vituse 
summary(ex_mod4bmi_fiber_calories_vituse)
cbind(beta = ex_mod4bmi_fiber_calories_vituse$coefficients, 
      confint(ex_mod4bmi_fiber_calories_vituse)) |>
  round(digits = 5)


#full model, order reflects sequence of trials 
names(carotene)
# calories + bmi + age + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse
null_lm<-lm(log(betaplasma) ~ 1, data = carotene)
modf4bmi<- lm(log(betaplasma)~ bmi, data= carotene)
modf4bmi_fiber<- lm(log(betaplasma)~ bmi+fiber, data= carotene) # age has more inf > fat
modf4bmi_fiber_cholesterol<- lm(log(betaplasma)~ bmi+fiber + cholesterol, data= carotene) # age has more inf > fat
modf4bmi_fiber_cholesterol_vituse<- lm(log(betaplasma)~ bmi+fiber + cholesterol + vituse, data= carotene) # age has more inf > fat
BIC(modf4bmi_fiber_cholesterol_vituse) # best one one can get

summary(modf4bmi_fiber_cholesterol_vituse)
cbind(beta = modf4bmi_fiber_cholesterol_vituse$coefficients, 
      confint(modf4bmi_fiber_cholesterol_vituse)) |>
  round(digits = 5)


# 
Model4full<-modf4bmi_fiber_cholesterol_vituse
Model4full_sum<-summary(Model4full)
Model4excl<-ex_mod4bmi_fiber_calories_vituse
Model4excl_sum<-summary(Model4excl)
# part 5 #
# 5 (a)

# For the reduced data set, determine whether it is necessary to have three categories in the categorical
# x-variable(s) in Model.4(b). 

# 5a: Model4full
summary(Model4full) # we have vituse as the categorical var
#
plyr::count(carotene, "vituse")

carotene <-  mutate(carotene, vit_yn = factor(vituse %in% c("Yes, not often", "Yes, fairly often"),
                         levels = c(FALSE, TRUE),
                         labels = c("No", "Yes")))

plyr::count(carotene, "vit_yn")
#
carotene <-  mutate(carotene, vitlots_yn = factor(vituse %in% c("Yes, not often", "No"),
                                              levels = c(FALSE, TRUE),
                                              labels = c("No", "Yes")))
plyr::count(carotene, "vitlots_yn")

# Tests

Model4full_vityn <- lm(log(betaplasma)~ bmi+fiber + cholesterol + vit_yn, data= carotene) 
Model4full_vityn_sum <- summary(Model4full_vityn)
Model4full_vityn_sum
anova(Model4full_vityn, Model4full )
Model4full_vitlots <- lm(log(betaplasma)~ bmi+fiber + cholesterol + vitlots_yn, data= carotene) 
Model4full_vitlots_sum <- summary(Model4full_vitlots)
Model4full_vitlots_sum
anova(Model4full_vitlots,Model4full)

Model4full_sum$adj.r.squared
Model4full_vityn_sum$adj.r.squared # this is a bit better
Model4full_vitlots_sum$adj.r.squared # this makes no sense

AIC(Model4full)
BIC(Model4full)

AIC(Model4full_vityn)
BIC(Model4full_vityn)# BIC got reduced 

Model5full<-Model4full_vityn

# 95% confidence intervals for new estimates
summary(Model5full)
cbind(beta = Model5full$coefficients, 
      confint(Model5full)) |>
  round(digits = 5)


# 5a:Model4excl

Model4excl # again, vituse
data_excl
plyr::count(data_excl, "vituse")
data_excl <-  mutate(data_excl, vit_yn = factor(vituse %in% c("Yes, not often", "Yes, fairly often"),
                                              levels = c(FALSE, TRUE),
                                              labels = c("No", "Yes")))

plyr::count(data_excl, "vit_yn")
# Perform a suitable test and update the variable(s) and model as necessary.
Model4excl_vityn <-lm(log(betaplasma)~  bmi+fiber+ calories + vit_yn, data=data_excl)
Model4excl_vityn_sum <- summary(Model4excl_vityn)
Model4excl_vityn_sum
anova(Model4excl_vityn, Model4excl)
#
Model4excl_sum$adj.r.squared
Model4full_vityn_sum$adj.r.squared # this is better
AIC(Model4full)
BIC(Model4full)
AIC(Model4excl_vityn)
BIC(Model4excl_vityn)# BIC got reduced 

Model5excl<-Model4excl_vityn
Model5excl_sum<-summary(Model5excl)
# 95% confidence intervals for new estimates
summary(Model5excl)
cbind(beta = Model5excl$coefficients, 
      confint(Model5excl)) |>
  round(digits = 5)

#5(b)

Model1full <- lm(data = carotene, log(betaplasma) ~ bmi)
Model1excl <- lm(data = data_excl, log(betaplasma) ~ bmi)
Model1excl_sum<-summary(Model1excl)
#
Model2full <- lm(data = carotene, log(betaplasma)~smokstat)
Model2excl <- lm(data = data_excl, log(betaplasma)~smokstat)
Model2excl_sum<-summary(Model2excl)
#
Model3excl<-car_excl_lm # for notation consistency
Model3excl_sum<-summary(Model3excl)
# For the reduced data set, present a table with the number of-parameters, the residual standard
# deviation, the R2, adjusted R2, AIC and BIC for the all five models (1(b), 2(b), 3(c), 4(b), and 5(a)).
# State which model you find best, and the reasons for your choice

# Model1excl, Model2excl, Model3excl, Model4excl, Model5excl
collect.R2AIC <- data.frame(
  nr = seq(1, 5),
  modnames=c("1.bmi", "2.smk", "3.all", "4.3flvl", "5.2flvl"),
  num_of_betas=c(length(Model1excl$coefficients),
             length(Model2excl$coefficients),
             length(Model3excl$coefficients),
             length(Model4excl$coefficients),
             length(Model5excl$coefficients)),
  sigma= c(Model1excl_sum$sigma, #residual_stdev
                   Model2excl_sum$sigma,
                   Model3excl_sum$sigma,
                   Model4excl_sum$sigma,
                   Model5excl_sum$sigma),
  R2 = c(Model1excl_sum$r.squared,
       Model2excl_sum$r.squared,
       Model3excl_sum$r.squared,
       Model4excl_sum$r.squared,
       Model5excl_sum$r.squared),
  R2.adj = c(Model1excl_sum$adj.r.squared,
           Model2excl_sum$adj.r.squared,
           Model3excl_sum$adj.r.squared,
           Model4excl_sum$adj.r.squared,
           Model5excl_sum$adj.r.squared)) 

collect.R2AIC
cbind(collect.R2AIC,
      AIC = AIC(Model1excl, Model2excl, Model3excl, Model4excl, Model5excl),
      BIC = BIC(Model1excl, Model2excl, Model3excl, Model4excl, Model5excl)) -> collect.R2AIC
glimpse(collect.R2AIC)
collect.R2AIC |>
  dplyr::rename(df = AIC.df, AIC = AIC.AIC, BIC = BIC.BIC) |>
  mutate(BIC.df = NULL) -> collect.R2AIC

collect.R2AIC

collect.R2AIC |> arrange(BIC)

