######Install and load the necessary packages
install.packages("dplyr")
library(dplyr)

install.packages("tidyr")
library(tidyr)

install.packages("lme4")
library(lme4)

install.packages("VIM")
library(VIM)

install.packages("naniar")
library(naniar)

install.packages("mice")
library(mice)

install.packages("lmerTest")
library(lmerTest)

library("broom")
library("broom.mixed")

install.packages("ggplot2")
library(ggplot2)

install.packages("robustlmm")
library(robustlmm)


#####################################################################################

baseline <- read.csv('baseline.csv')
hscrp <- read.csv('hscrp.csv')

glimpse(hscrp)
non_numeric_values <- hscrp$hsCRP.result..m12[!grepl("^[0-9.]+$", hscrp$hsCRP.result..m12)]
non_numeric_values
hscrp$hsCRP.result..m12 <- gsub("[^0-9.]", "", hscrp$hsCRP.result..m12)
hscrp$hsCRP.result..m12 <- as.numeric(hscrp$hsCRP.result..m12)
View(hscrp)

# Joining the datasets on 'Subject ID'
df_c <- left_join(baseline, hscrp, by = 'Subject.ID')
View(df_c)
glimpse(df_c)


#######################################################################################

##Preprocessing (Data Cleaning)

#renaming columns in combined dataframe
colnames(df_c) <- tolower(colnames(df_c))
df_c <- df_c %>% rename(
  subject_id = `subject.id`,
  rand_group = `randomised.group`,
  acs_type = `type.of.acs...stratifying.variables`,
  cd8_prop = `proportion.of.cd8..temra.at.baseline...stratifying.variables`,
  gender = `gender...stratifying.variables`
)

#convert subject_id to factors in combined df
df_c$subject_id <- as.factor(df_c$subject_id)
df_c$rand_group <- as.factor(df_c$rand_group)
df_c$acs_type <- as.factor(df_c$acs_type)
df_c$cd8_prop <- as.factor(df_c$cd8_prop)
df_c$gender <- as.factor(df_c$gender)

#dropping columns in combined dataframe
df_c <- df_c %>% select(-aged.65.or.over.with.an.index.presentation.of.an.acute.coronary.syndrome..within.the.previous.6.months...inclusion.criteria,
                      -written.informed.consent...inclusion.criteria
)
View(df_c)
summary(df_c)

#Convert to long format for longitudinal analysis
df <- df_c %>%
  pivot_longer(
    cols = starts_with("hscrp.result"), 
    names_to = "time", 
    values_to = "hscrp",
    names_prefix = "hscrp result -m"
  ) %>%
  mutate(
    time = case_when(
      time == "hscrp.result..m0"  ~ 0,
      time == "hscrp.result..m6"  ~ 6,
      time == "hscrp.result..m12" ~ 12,
      TRUE ~ NA_real_ # Handles any unexpected values
    )
  )


#convert time to factor and insert as a new column long df
df <- df %>%
  mutate(timef = as.factor(df$time))

# Renaming '0' with 'baseline'
df$timef <- as.character(df$timef)
df$timef[df$timef == '0'] <- 'baseline'
df$timef[df$timef == '6'] <- '6 months'
df$timef[df$timef == '12'] <- '12 months'
df$timef <- as.factor(df$timef)

#ANALYSING OUTLIERS
ggplot(df_c, aes(x = hscrp.result..m0)) +
  geom_histogram(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP levels at baseline",
    x = "hsCRP",
    y = "Count"
  )

ggplot(df_c, aes(x = hscrp.result..m6)) +
  geom_histogram(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP levels at 6 months",
    x = "hsCRP",
    y = "Count"
  )
ggplot(df_c, aes(x = hscrp.result..m12)) +
  geom_histogram(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP levels at 12 months",
    x = "hsCRP 12",
    y = "Count"
  )


ggplot(df_c, aes(x = hscrp.result..m0)) +
  geom_boxplot(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP levels at baseline",
    x = "hsCRP  baseline",
    y = "Count"
  )

ggplot(df_c, aes(x = hscrp.result..m6)) +
  geom_boxplot(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP levels at 6 months",
    x = "hsCRP 6 Months",
    y = "Count"
  )
ggplot(df_c, aes(x = hscrp.result..m12)) +
  geom_boxplot(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP levels at 12 months",
    x = "hsCRP 12 months",
    y = "Count"
  )

#before winsorsing
ggplot(df, aes(x = hscrp)) +
  geom_boxplot(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP Levels before winsorisation",
    x = "hsCRP",
    y = "Count"
  )

winsorise <- function(x, lower_quantile = 0.05, upper_quantile = 0.95) {
  lower_bound <- quantile(x,lower_quantile,na.rm = TRUE)
  upper_bound <- quantile(x,upper_quantile,na.rm = TRUE)
  x[x < lower_bound] <- lower_bound
  x[x > upper_bound] <- upper_bound
  return (x)
}

df$hscrp_w <- winsorise(df$hscrp, lower_quantile = 0.05, upper_quantile = 0.95)
dfg <- df
View(dfg)
#after winsorising
ggplot(df, aes(x = hscrp_w)) +
  geom_boxplot(fill = 'blue') +
  labs(
    title = "Distribution of hsCRP Levels after winsorisation",
    x = "hsCRP",
    y = "Count"
  )
View(df)
##Taking log of hscrp levels
df <- df %>%
  mutate(hscrp_logw = log(hscrp_w))
dfg <- df
df <- df %>%
  mutate(hscrp_log = log(hscrp))
#after logging hsCRP
ggplot(df, aes(x = hscrp_logw)) +
  geom_boxplot(fill = 'blue') +
  labs(
    title = "Distribution of logged hsCRP levels",
    x = "hsCRP",
    y = "Count"
  )


#testing time as a factor and continuous variable
lmm_tf <- lmer(hscrp_log ~ rand_group + timef + ( 1 | subject_id), data = df)
lmm_tc <-lmer(hscrp_log ~ rand_group + time + (1 | subject_id), data = df)
anova(lmm_tf,lmm_tc)
View(df)

##testing impact of winsorising
ww <- lmer(hscrp_logw ~ rand_group + timef + (1 | subject_id), data = df)
summary(ww)

wtw <- lmer(hscrp_log ~ rand_group + timef + (1 | subject_id), data = df)
summary(wtw)

AIC(ww, wtw)
BIC(ww,wtw)


View(df)

######################################################################################################

##########MISSING DATA ANALYSIS

#vISUALISE MISSING DATA
gg_miss_upset(df_c)

colSums(is.na(df_c))


#gender missing data
missing_data <- df_c %>%
  group_by(gender) %>%
  summarise_all(~ sum(is.na(.))) %>%
  gather(key = "variable", value = "missing_count", -gender)

ggplot(missing_data, aes(x = gender, y = missing_count)) +
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Missing Values by Gender",
    x = "Gender",
    y = "Number of Missing Values"
  )
#cd8_prop missing data
missing_data_cd8 <- df_c %>%
  group_by(cd8_prop) %>%
  summarise_all(~ sum(is.na(.))) %>%
  gather(key = "variable", value = "missing_count", -cd8_prop)

ggplot(missing_data_cd8, aes(x = cd8_prop, y = missing_count)) +
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Missing Values by cd8_prop",
    x = "cd8_prop",
    y = "Number of Missing Values"
  )

#rand_group missing data
missing_data_rand_group <- df_c %>%
  group_by(rand_group) %>%
  summarise_all(~ sum(is.na(.))) %>%
  gather(key = "variable", value = "missing_count", -rand_group)

ggplot(missing_data_rand_group, aes(x = rand_group, y = missing_count)) +
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Missing Values by rand_group",
    x = "rand_group",
    y = "Number of Missing Values"
  )

#acs missing data
missing_data_acs <- df%>%
  group_by(acs_type) %>%
  summarise_all(~ sum(is.na(.))) %>%
  gather(key = "variable", value = "missing_count", -acs_type)

View(missing_data_acs)
ggplot(missing_data_acs, aes(x = acs_type, y = missing_count)) +
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Missing Values by ACS Type",
    x = "ACS Type",
    y = "Number of Missing Values"
  )

colSums(is.na(df_c[df_c$gender == 'Male',]))
colSums(is.na(df_c[df_c$gender == 'Female',]))

colSums(is.na(df_c[df_c$cd8_prop == 'High >45%',]))
colSums(is.na(df_c[df_c$cd8_prop == 'Low ≤45%',]))

colSums(is.na(df_c[df_c$rand_group == 'Placebo',]))
colSums(is.na(df_c[df_c$rand_group == 'TA65',]))

colSums(is.na(df_c[df_c$acs_type == 'STEMI',]))
colSums(is.na(df_c[df_c$acs_type == 'NSTEMI',]))

#####
colSums(is.na(df[df$gender == 'Male',]))
colSums(is.na(df[df$gender == 'Female',]))

colSums(is.na(df[df$cd8_prop == 'High >45%',]))
colSums(is.na(df[df$cd8_prop == 'Low ≤45%',]))

colSums(is.na(df[df$rand_group == 'Placebo',]))
colSums(is.na(df[df$rand_group == 'TA65',]))

colSums(is.na(df[df$acs_type == 'STEMI',]))
colSums(is.na(df[df$acs_type == 'NSTEMI',]))

colSums(is.na(df_c))

#######HANDLING MISSING VALUES

###MULTIPLE IMPUTATION
results <- list()

m_values <- c(5, 10, 20, 30, 40)

df$timef <- factor(df$timef, levels = c("baseline", "6 months", "12 months"))
for (m_value in m_values) {
  imp <- mice(df, method = 'pmm', m = m_value, seed = 500)
  fit <- with(imp, lmer(hscrp_logw ~ rand_group + timef + (1 | subject_id)))  # Adjust the formula as per your model
  pooled_results <- pool(fit)
  results[[as.character(m_value)]] <- summary(pooled_results)
}

print(pooled_results)

estimates <- data.frame(
  m = rep(m_values, each = 4),  # 'm' values
  term = rep(c('(Intercept)', 'rand_groupTA65', 'timef6','timef12'), times = length(m_values)),
  estimate = unlist(lapply(results, function(x) x$estimate)),
  std.error = unlist(lapply(results, function(x) x$std.error))
)

ggplot(estimates, aes(x = m, y = estimate, color = term, group = term)) +
  geom_line() +
  geom_point() +
  labs(title = "Stability of Estimates Across Different m Values for LMER",
       x = "Number of Imputations (m)",
       y = "Estimate") +
  theme_minimal()

df_imp <- mice(df, m = 20, method = 'pmm', seed = 500)
df_i <- complete(df_imp,1)
View(df)
#LINEAR MIXED MODELS TO IMPUTE MISSING VALUES
df_t <- df
df_t$timef <- factor(df$timef, levels = c("baseline", "6 months", "12 months"))
lmm <- lmer(hscrp_logw ~ rand_group + timef + (1 | subject_id), data = df_t, na.action = na.omit)
summary(lmm)
df_t$hscrp_logw[is.na(df_t$hscrp_logw)] <- predict(lmm, newdata=df_t[is.na(df_t$hscrp_logw), ])
df_t$hscrp_w <- exp(df_t$hscrp_logw)
View(df_t)

#testing lmm imputation vs multiple imputation
lmm <- lmer(hscrp_logw ~ rand_group + timef  + (1 | subject_id), data = df_t)
summary(lmm)

mi <- lmer(hscrp_logw ~ rand_group + timef + (1 | subject_id), data = df_i)
summary(mi)

# AIC comparison
AIC(lmm,mi)

# BIC comparison
BIC(lmm, mi)

summary(df$hscrp_w)  # Original data with missing values
summary(df_i$hscrp_w)  # After multiple imputation
summary(df_t$hscrp_w) # Mixed model imputation

boxplot(df$hscrp_w, df_i$hscrp_w, df_t$hscrp_w, 
        names=c("Original", "Multiple Imputation", "Mixed Model Imputation"), 
        main="Boxplot Comparison of Imputation Methods", col=c("blue", "red", "green"))


#######################################################################################

##Exploratory Data Analysis (EDA)
summary(df_t)
View(df_t)

#Log distribution of hsCRP
ggplot(df_t[df_t$timef == 'baseline',], aes(x = hscrp_logw)) +
  geom_histogram(fill = 'blue') +
  labs(
    title = "Distribution of log hsCRP levels at baseline",
    x = "hsCRP",
    y = "Count"
  )

ggplot(df_t[df_t$timef == "6 months",], aes(x = hscrp_logw)) +
  geom_histogram(fill = 'blue') +
  labs(
    title = "Distribution of log hsCRP levels at baseline",
    x = "hsCRP",
    y = "Count"
  )

ggplot(df_t[df_t$timef == "12 months",], aes(x = hscrp_logw)) +
  geom_histogram(fill = 'blue') +
  labs(
    title = "Distribution of log hsCRP levels at baseline",
    x = "hsCRP",
    y = "Count"
  )

View(df)
######## hscrp values for TA65 and placebo

df_t$timef <- factor(df_t$timef, levels = c("baseline", "6 months", "12 months"))
ggplot(df_t, aes(x = timef, y = hscrp_w)) +
  geom_point( color = 'blue') +
  labs(title = "hsCRP trend over time" , x = "time") +
  theme_minimal()

df_ts <- df_t[df_t$rand_group == 'Placebo',]
df_ts <- df_ts[37:48,]
View(df_ts)
df_ts$timef <- factor(df_ts$timef, levels = c("baseline", "6 months", "12 months"))
ggplot(df_ts, aes(x = timef, y = hscrp_w, color = rand_group, group = subject_id)) +
  geom_line(color = 'red') +  # Line for each subject
  geom_point( color = 'red') +  # Points for each observation
  labs(title = "hsCRP Trends Over Time(Placebo)",
       x = "Time (Months)", 
       y = "hsCRP Levels (mg/L)") +
  theme_minimal() +
  theme(legend.title = element_blank())

df_ts <- df_t[df_t$rand_group == 'TA65',]
df_ts <- df_ts[37:48,]
df_ts$timef <- factor(df_ts$timef, levels = c("baseline", "6 months", "12 months"))
ggplot(df_ts, aes(x = timef, y = hscrp_w, color = rand_group, group = subject_id)) +
  geom_line( color = 'blue') +  # Line for each subject
  geom_point( color = 'blue') +  # Points for each observation
  labs(title = "hsCRP Trends Over Time(TA65)",
       x = "Time (Months)", 
       y = "hsCRP Levels (mg/L)") +
  theme_minimal() +
  theme(legend.title = element_blank())


## hsCRP by Randomised group
df$timef <- factor(df$timef, levels = c("baseline", "6 months", "12 months"))
ggplot(df_t, aes(x =timef, y = hscrp, fill = rand_group)) +
  geom_boxplot() +
  labs(title = "Boxplot of for each Treatment Group", x = "hsCRP Levels", y = "Value") +
  scale_y_log10() +
  scale_fill_manual(values = c("Placebo" = "red", "TA65" = "blue")) +
  theme_minimal()

View(df_t)
##DISTRIBUTION OF HSCRP ACROSS TIME FOR DIFFERENT GROUPS
#ACS_TYPE
mean_hscrp_at <- df_t[df_t$acs_type == 'NSTEMI',] %>%
  group_by(acs_type, timef, rand_group) %>%  # Include rand_group in grouping
  summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)

ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) +  # Line connecting the means
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels Over Time for NSTEMI",
       x = "Time (Months)",
       y = "Mean hsCRP Levels",
       color = "ACS Type") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)

mean_hscrp_at <- df_t[df_t$acs_type == 'STEMI',] %>%
  group_by(acs_type, timef, rand_group) %>%  # Include rand_group in grouping
  summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)

ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) +  # Line connecting the means
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels Over Time for STEMI",
       x = "Time (Months)",
       y = "Mean hsCRP Level",
       color = "ACS Type") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)

#GENDER

mean_hscrp_at <- df_t[df_t$gender == 'Male',] %>%
  group_by(gender, timef, rand_group) %>% 
  summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)

ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels for Male",
       x = "Time (Months)",
       y = "Mean hsCRP Level",
       color = "gender") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)

mean_hscrp_at <- df_t[df_t$gender == 'Female',] %>%
  group_by(gender, timef, rand_group) %>%  
  summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)

ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels for Female",
       x = "Time (Months)",
       y = "Mean hsCRP Level",
       color = "gender") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)

#cd8_prop
mean_hscrp_at <- df_t[df_t$cd8_prop == 'High >45%',] %>%
  group_by(cd8_prop, timef, rand_group) %>%  # Include rand_group in grouping
  summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)

ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels for High Cd8 Proportion",
       x = "Time (Months)",
       y = "Mean hsCRP Level",
       color = "gender") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)

mean_hscrp_at <- df_t[df_t$cd8_prop == 'Low ≤45%',] %>%
  group_by(cd8_prop, timef, rand_group) %>%  # Include rand_group in grouping
  summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)

ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels for Low Cd8 Proportion",
       x = "Time (Months)",
       y = "Mean hsCRP Level",
       color = "gender") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)


##MEAN HSCRP LEVELS rand_group
mean_hscrp_at <- df_t %>%
  group_by(rand_group, timef) %>%
  dplyr::summarize(mean_hscrp = mean(hscrp_w, na.rm = TRUE))

print(mean_hscrp_at)


ggplot(mean_hscrp_at, aes(x = timef, y = mean_hscrp, color = rand_group, group = rand_group)) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  labs(title = "Mean hsCRP Levels Randomised Groups",
       x = "Time (Months)",
       y = "Mean hsCRP Level",
       color = "Randomised Group") +
  scale_color_manual(values = c("Placebo" = "blue", "TA65" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin = mean_hscrp - 0.1 * mean_hscrp, ymax = mean_hscrp + 0.1 * mean_hscrp), width = 0.2)



#######################################################################################



##LINEAR MIXED EFFECT MODELS WITH MIXED MODEL IMPUTED DF


#BASELINE MODEL
base <- lmer(hscrp_logw ~ rand_group + timef + (1 | subject_id), data = df_t)
summary(base)

base_wt <- lmer(hscrp ~ rand_group + as.factor(time) + (1 | subject_id), data = dfb)
summary(base_wt)

AIC(base,base_wt)

# Calculate percent change
percent_change <- (exp(-0.03705) - 1) * 100
percent_change

plot(resid(base), col = 'blue')
qqnorm(resid(base), col = 'blue')
qqline(resid(base), col = 'red')

tidy(baseline)
ranef(baseline)
fixef(baseline)
confint(baseline)

coef_estimates <- 
  tidy(baseline,conf.int = TRUE) %>%
  filter(effect == "fixed")

print(coef_estimates)

ggplot(coef_estimates, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline( yintercept = 0, color = 'red') +
  geom_linerange() + geom_point() + coord_flip() + theme_minimal()


AIC(base,lmm1,lmm2,lmm3)

##LMM 1
lmm1 <-lmer(hscrp_logw ~ rand_group + timef + gender + (1 | subject_id), data = df_t)
summary(lmm1)

plot(resid(lmm1))
qqnorm(resid(lmm1))
qqline(resid(lmm1), col = 'red')


tidy(lmm1)
ranef(lmm1)
fixef(lmm1)
confint(lmm1)

coef_estimates <- 
  tidy(lmm1,conf.int = TRUE) %>%
  filter(effect == "fixed")

print(coef_estimates)

ggplot(coef_estimates, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline( yintercept = 0, color = 'red') +
  geom_linerange() + geom_point() + coord_flip() + theme_minimal()

anova(base, lmm1)


##LMM 2
lmm2 <-lmer(hscrp_logw ~ rand_group + timef + acs_type + (1 | subject_id), data = df_t)
summary(lmm2)

anova(base, lmm2)


##LMM 3
lmm3 <-lmer(hscrp_logw ~ rand_group + timef +  cd8_prop +(1 | subject_id), data = df_t)
summary(lmm3)
anova(base, lmm3)

##LMM 4
lmm4 <-lmer(hscrp_logw ~ rand_group + timef + cd8_prop + acs_type + (1 | subject_id), data = df_t)
summary(lmm4)
anova(base,lmm4)


#####
#lmm5
lmm5 <-lmer(hscrp_logw ~ rand_group + timef + gender + acs_type + (1 | subject_id), data = df_t)
summary(lmm5)
anova(base,lmm5)

#lmm6
lmm6 <-lmer(hscrp_logw ~ rand_group + timef + gender + acs_type + cd8_prop + (1 | subject_id), data = df_t)
summary(lmm6)

anova(base,lmm6)



AIC(base,lmm1,lmm2,lmm3,lmm4,lmm5,lmm6)

BIC(base,lmm1,lmm2,lmm3,lmm4,lmm5,lmm6)

percent_change <- (exp(-0.04641) - 1) * 100
percent_change




##EXPERIMENTING WITH INTERACTIONS
#LMM
lmm7 <- lmer(hscrp_logw ~ rand_group * timef + ( 1 | subject_id), data = df_t)
summary(lmm7)

anova(base,lmm5)

#LMM6 
lmm8 <- lmer(hscrp_logw ~ rand_group * timef + gender + ( 1 | subject_id), data = df_t)
summary(lmm8)

anova(base,lmm6)


#LMM9
lmm9 <- lmer(hscrp_logw ~ rand_group * timef + acs_type + ( 1 | subject_id), data = df_t)
summary(lmm9)

anova(base,lmm7)

#LMM10
lmm10 <- lmer(hscrp_logw ~ rand_group * timef + cd8_prop + ( 1 | subject_id), data = df_t)
summary(lmm10)

anova(base,lmm8)

#LMM11 
lmm11 <- lmer(hscrp_logw ~ rand_group * timef + cd8_prop + acs_type + ( 1 | subject_id), data = df_t)
summary(lmm11)

anova(base,lmm9)


#LMM12 
lmm12 <- lmer(hscrp_logw ~ rand_group * timef + gender + acs_type + ( 1 | subject_id), data = df_t)
summary(lmm12)

anova(base,lmm9)

#LMM13 
lmm13 <- lmer(hscrp_logw ~ rand_group * timef + gender + acs_type + cd8_prop +( 1 | subject_id), data = df_t)
summary(lmm13)

AIC(lmm7,lmm8,lmm9,lmm10,lmm11,lmm12,lmm13)
BIC(lmm7,lmm8,lmm9,lmm10,lmm11,lmm12,lmm13)


###PLOTTING LMM
plot(resid(base),col = 'blue', main = "Residuals Plot")
qqnorm(resid(base), col = 'blue')
qqline(resid(base), col = 'red')


tidy(base)
ranef(base)
fixef(base)
confint(base)

coef_estimates <- 
  tidy(lmm1,conf.int = TRUE) %>%
  filter(effect == "fixed")

print(coef_estimates)

ggplot(coef_estimates, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline( yintercept = 0, color = 'red') +
  geom_linerange() + geom_point() + coord_flip() + theme_minimal()



#########################################################



#### PLOTTING MODELS
tidy(lmm2)
ranef(lmm2)
fixef(lmm2)
confint(lmm2)

coef_estimates <- 
  tidy(lmm2,conf.int = TRUE) %>%
  filter(effect == "fixed")

print(coef_estimates)

ggplot(coef_estimates, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline( yintercept = 0, color = 'red') +
  geom_linerange() + geom_point() + coord_flip() + theme_minimal()

##EXPERIMENTING WITH GENERALISED LINEAR MIXED EFFECTS MODEL

# Using GLMM with Gamma distribution and log link
#GLMM1
glmm1 <- glmer(hscrp_w ~ rand_group + timef + (1 | subject_id),family = Gamma(link = "log"), data = df_t)
summary(glmm1)

View(df_t)

#GLMM2
glmm2 <- glmer(hscrp_w ~ rand_group * timef + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm2)

anova(glmm1,glmm2)


#GLMM3
glmm3 <- glmer(hscrp_w ~ rand_group + timef + gender + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm3)

#GLMM4
glmm4 <- glmer(hscrp_w ~ rand_group + timef + acs_type + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm4)

#GLMM5
glmm5 <- glmer(hscrp_w ~ rand_group + timef + cd8_prop + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm5)

#GLMM6
glmm6 <- glmer(hscrp ~ rand_group + timef + gender + acs_type + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm6)

#GLMM7
glmm7 <- glmer(hscrp ~ rand_group + timef + gender + acs_type + cd8_prop + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm7)

#GLMM8
glmm8 <- glmer(hscrp ~ rand_group + timef + acs_type + cd8_prop + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm8)

View(df_t)

##interactions
#GLMM9
glmm9 <- glmer(hscrp_w ~ rand_group * timef +(1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm9)

#GLMM10
glmm10 <- glmer(hscrp_w ~ rand_group * timef + gender + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm10)

#GLMM11
glmm11 <- glmer(hscrp_w ~ rand_group * timef + cd8_prop + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm11)

#GLMM12
glmm12 <- glmer(hscrp_w ~ rand_group * timef + gender + acs_type + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm12)

#GLMM13
glmm13 <- glmer(hscrp_w ~ rand_group * timef + acs_type + cd8_prop + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm13)

#GLMM14
glmm14 <- glmer(hscrp_w ~ rand_group * timef + gender + acs_type + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm14)

#GLMM15
glmm15 <- glmer(hscrp_w ~ rand_group * timef + gender + acs_type + cd8_prop + (1 | subject_id), family = Gamma(link = "log"), data = df_t)
summary(glmm15)

anova(glmm6,glmm6)

##Plotting Models

plot(resid(glmm1))
qqnorm(resid(glmm1))
qqline(resid(glmm1), col = 'red')


tidy(glmm1)
ranef(glmm1)
fixef(glmm1)
confint(glmm1)

coef_estimates <- 
  tidy(glmm1,conf.int = TRUE) %>%
  filter(effect == "fixed")

print(coef_estimates)

ggplot(coef_estimates, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline( yintercept = 0, color = 'red') +
  geom_linerange() + geom_point() + coord_flip() + theme_minimal()

anova(base, lmm1)










