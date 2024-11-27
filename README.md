# Analysing-effect-of-TA65-on-hsCRP-levels-in-post-myocardial-infarction-patients

# PROBLEM STATEMENT
The dataset used is from a clinical trial at the Newcastle Clinical Trials where a double-blinded randomised controlled pilot trial evaluated the use of TA-65 to reduce cell aging in patients following Myocardial Infarction. There were 90 patients (45 Male, 45 Female) all aged over 65. In the trial, 45 random patients were given the placebo while the others were given TA-65. Data collection was performed over a 12-month period, with regular follow-ups to monitor changes in immune cell counts and hsCRP levels. The outcome of the trial High-sensitivity C-reactive protein (hsCRP) levels were recorded at baseline, 6 months, and 12 months. The task at hand is to model this data using a Linear Mixed Effect Model in other study the effect at 6 months and 12 months


# KEY INSIGHTS TO INFORM DECISION MAKING
In this study, we investigated the effect of TA65 on hsCRP levels, a marker of inflammation, over 6 and 12 months.Both the TA65 and placebo groups experienced significant reductions in hsCRP levels over time, with decreases of approximately 70% at 6 months and 72.5% at 12 months compared to baseline. At 6 months, the TA65 group had a 5.5% more improvement in hsCRP levels compared to the placebo group. At 12 months, the TA65 group had an 8% less improvement compared to the placebo group. However, these differences were not significant, meaning there is no clear evidence that TA65 worked better or worse than the placebo in reducing inflammation.

# DATA CLEANING AND PREPARATION
To facilitate data manipulation and analysis, several R packages were installed and loaded, including dplyr and tidyr for data wrangling, lme4 and lmerTest for fitting linear mixed-effects models, VIM and mice for missing data visualisation and imputation, ggplot2 for data visualisation, broom and broom.mixed for tidy output of model results.
After baseline and hsCRP datasets were loaded, any non-numeric values in hsCRP columns were identified and cleaned using regular expressions to remove unwanted characters and convert the values to numeric types. The baseline dataset was merged with the hsCRP dataset based on the common column ‘Subject ID’ to ensure all baseline information was retained, while hsCRP measurements from multiple time points where added. For longitudinal analysis the new combined dataset was reshaped from wide to long format using the pivot_longer() function in R in other to allow the hsCRP measurements at different time points to be analysed as repeated measures. The time variable was transformed to represent the actual points ( 0, 6 and 12 months). To enhance consistent referencing throughout the analysis column names were standardised to lowercase. Data types where checked and converted to the correct format, subject_id , rand_group, gender, acs_type, and cd8_prop were converted to factors.**Data in long format can be seen below**
![image](https://github.com/user-attachments/assets/afc9e7a6-74f6-40f0-83a7-469ab6e3d792)


## Summary Statistics
![image](https://github.com/user-attachments/assets/e5c890aa-b1e5-4e74-9ed1-81304dd705bf)
Analysing the summary statistics above, we can see that the placebo group and TA65 contains 45 observations which indicates that the data is nicely balanced based on treatment groups. On the other hand there are predominantly more males than females in the dataset indicating an imbalance in gender, NSTEMI and STEMI groups quite balanced with 47 and 43 observations each. Post-myocardial patients with a high cd8_prop level( > 45%) are more frequent than patients with low cd8_prop(<= 45%), indicating an imbalance in cd8_prop.Observing the mean and median of hsCRP levels at baseline, 6 months and 12 months, it can be seen that they differ significantly indicating the presence of outliers which can also be observed in the max values of each timepoint.


## Handling Outliers
Analysing the outliers through a boxplot, we can observe that there a hsCRP levels as high as 153 mg/L which could indicate patients with a seriously high level of inflammation due to serious infection or underlying disease. In line with the methodology, winsorsing was carried out to mitigate the impact of the outliers on the mixed model. The winsorised hsCRP levels were also log transformed to reduce the impact of the remaining outliers. Below is a boxplot showing distribution of hsCRP level before and after winsorising.
![image](https://github.com/user-attachments/assets/078c390b-78c6-4e07-8c70-46bd950485cb) ![image](https://github.com/user-attachments/assets/8c050875-5e76-4693-9501-94d8591dd0af) ![image](https://github.com/user-attachments/assets/e80166a2-200c-4d39-b313-375e44b00846)

## Handling Missing Values
The dataset contains 52 missing values which is 20% of the dataset. In the figure below the missing values pattern is observed. One interesting observation is that there were 18 patients who’s hsCRP levels where not recorded for in month 6 and 12 while 1 patient did not have a record for baseline and month 12.
![image](https://github.com/user-attachments/assets/9394c854-04c1-4775-b98c-b0b722f0c0a4) .

Observing the missing data we can see that patients with a high cd8 proportion and the male gender tend to have a higher probability of having a missing value, this suggests the data could be missing at random. According to the methodology multiple imputation was carried out as well as using a linear mixed model to impute missing values. The linear mixed model had randomised group and time specified as fixed effects with subject as random effects, the outputted variance and residual variance were 0.01066 (std.dev 0.1033) and 1.06490 (std.dev 1.0319) respectively. Multiple imputation was implemented using the mice package in R with 10, 20, 30 and 40 imputations and estimates were pooled. The estimates were observed to be showing signs of stability at 20 imputations, thus 20 imputations was used going forth in the experiment.

![image](https://github.com/user-attachments/assets/f975ab7f-4adc-4e3a-b53f-f5c5d27089a1)

# EXPLORATORY DATA ANALYSIS
## Exploring the distribution of hsCRP
Exploring the distribution of the hsCRP levels at baseline, 6 months and 12 months, it is observed that the distribution is heavily tailed (right skewed), this could be as a result of situations where hsCRP levels were high at baseline. Historgram of the distribution of hsCRP across baseline, 6months and 12 months are shown below.
![image](https://github.com/user-attachments/assets/773bbd34-85b8-4aab-baf9-677785f13e8c)
![image](https://github.com/user-attachments/assets/d036c915-77ef-4054-a117-aff862b58237)
![image](https://github.com/user-attachments/assets/eaf662a2-e8c2-4189-9f6c-d53bcdb161f6)

Stratifying by randomised group we can clearly visualise the changes in hsCRP levels over time in TA65 group and Placebo group. In Fig we can observe that hsCRP levels at baseline, both treatment groups show higher hsCRP levels, with the Placebo group exhibiting a wider interquartile range, suggesting greater variability in the response among these participants. Over time, a notable decline in hsCRP levels is observed for both groups; however, the TA65 group consistently shows lower hsCRP levels compared to the Placebo group at each time point. This trend is particularly pronounced at the 6-month and 12-month marks, where the median hsCRP levels for TA65 approach lower values, indicating a potential treatment effect. The presence of outliers in both groups suggests individual variability in responses, but the overall reduction in hsCRP levels may indicate the effectiveness of the TA65 treatment over time.

![image](https://github.com/user-attachments/assets/c78607ce-8a90-4181-942b-317bb77717e4)

## Random Sampling
Furthermore, taking a random sample from the placebo and TA65 group to have a more intimate visualisation of the trend, we can observe in fig… and fig . that hsCRP generally decrease with time for both groups.

![image](https://github.com/user-attachments/assets/9bf391d9-cabe-4b1f-8967-8452482854b6)

![image](https://github.com/user-attachments/assets/6694c94d-e57d-4348-9cc0-63fe79dc248d)

## Exploring the mean of each group
The mean of hsCRP for placebo group and TA65 group were explored for each fixed effect.
![image](https://github.com/user-attachments/assets/8f7f9df2-4ccc-4a58-b311-06f021832ddf)


# MODELLING DATA
The lmer package from R. Linear mixed models were fitted using the lmer() function.The models were incrementally built, starting with a simple random intercept model and progressing to interaction models. The models were employed to analyse the clinical trial specifically focusing on hsCRP levels measured at baseline, 6 months, and 12 months.
In our case, the response variable is the log-transformed hsCRP levels (hscrp_log), the fixed effects include treatment group (rand_group), time (timef), and other covariates such as gender, ACS type (acs_type), and CD8 proportion (cd8_prop), while the random effects account for the variability among individual subjects (subject_id).
ALinear Mixed Effect Model was fitted with the interaction terms specified for the randomised group and time.
# RESULT AND EVALUATION
The results revealed a significant reductions in hsCRP levels at both 6 months (Estimate: -1.20233, p < 0.001) and 12 months (Estimate: -1.29079, p < 0.001) as expected. The interaction terms for the TA65 treatment, however, did not reach significance (6 months: Estimate: 0.05380, p =0.862; 12 months: Estimate: -0.08322, p = 0.788). The TA65 treatment effect’s estimate of -0.06766 in suggests a 6.6% decrease in hsCRP levels compared to the placebo group, though this was not statistically significant (p = 0.758).
