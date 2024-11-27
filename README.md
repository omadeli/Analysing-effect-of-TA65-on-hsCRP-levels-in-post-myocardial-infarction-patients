# Analysing-effect-of-TA65-on-hsCRP-levels-in-post-myocardial-infarction-patients

# PROBLEM STATEMENT
The dataset used is from a clinical trial at the Newcastle Clinical Trials where a double-blinded randomised controlled pilot trial evaluated the use of TA-65 to reduce cell aging in patients following Myocardial Infarction. There were 90 patients (45 Male, 45 Female) all aged over 65. In the trial, 45 random patients were given the placebo while the others were given TA-65. Data collection was performed over a 12-month period, with regular follow-ups to monitor changes in immune cell counts and hsCRP levels. The outcome of the trial High-sensitivity C-reactive protein (hsCRP) levels were recorded at baseline, 6 months, and 12 months. The task at hand is to model this data using a Linear Mixed Effect Model in other study the effect at 6 months and 12 months


# KEY INSIGHTS TO INFORM DECISION MAKING

# DATA CLEANING AND PREPARATION
To facilitate data manipulation and analysis, several R packages were installed and loaded, including dplyr and tidyr for data wrangling, lme4 and lmerTest for fitting linear mixed-effects models, VIM and mice for missing data visualisation and imputation, ggplot2 for data visualisation, broom and broom.mixed for tidy output of model results.
After baseline and hsCRP datasets were loaded, any non-numeric values in hsCRP columns were identified and cleaned using regular expressions to remove unwanted characters and convert the values to numeric types. The baseline dataset was merged with the hsCRP dataset based on the common column ‘Subject ID’ to ensure all baseline information was retained, while hsCRP measurements from multiple time points where added. For longitudinal analysis the new combined dataset was reshaped from wide to long format using the pivot_longer() function in R in other to allow the hsCRP measurements at different time points to be analysed as repeated measures. The time variable was transformed to represent the actual points ( 0, 6 and 12 months). To enhance consistent referencing throughout the analysis column names were standardised to lowercase. Data types where checked and converted to the correct format, subject_id , rand_group, gender, acs_type, and cd8_prop were converted to factors.


# MODELLING DATA
The lmer package from R. Linear mixed models were fitted using the lmer() function.The models were incrementally built, starting with a simple random intercept model and progressing to more complex models. The models were employed to analyse the clinical trial specifically focusing on hsCRP levels measured at baseline, 6 months, and 12 months.
In our case, the response variable is the log-transformed hsCRP levels (hscrp_log), the fixed effects include treatment group (rand_group), time (timef), and other covariates such as gender, ACS type (acs_type), and CD8 proportion (cd8_prop), while the random effects account for the variability among individual subjects (subject_id).

The base model was tweaked to include interactions for the randomised group and time. revealing significant reductions in hsCRP levels at both 6 months (Estimate: -1.20233, p < 0.001) and 12 months (Estimate: -1.29079, p < 0.001) as expected. The interaction terms for the TA65 treatment, however, did not reach significance (6 months: Estimate: 0.05380, p =0.862; 12 months: Estimate: -0.08322, p = 0.788). The TA65 treatment effect’s estimate of -0.06766 in suggests a 6.6% decrease in hsCRP levels compared to the placebo group, though this was not statistically significant (p = 0.758).
# RESULT AND EVALUATION
