# Analysing effect of TA65 on High-Sensitivity C-reative protien levels in post myocardial infarction patients

# PROBLEM STATEMENT
The dataset used is from a clinical trial at the Newcastle Clinical Trials where a double-blinded randomised controlled pilot trial evaluated the use of TA-65 to reduce cell aging in patients following Myocardial Infarction. There were 90 patients (45 Male, 45 Female) all aged over 65. In the trial, 45 random patients were given the placebo while the others were given TA-65(an oral telomerase activator). Data collection was performed over a 12-month period, with regular follow-ups to monitor changes in immune cell counts and hsCRP levels. The outcome of the trial, High-sensitivity C-reactive protein (hsCRP) levels were recorded at baseline, 6 months, and 12 months. The task at hand is to model this data using a Linear Mixed Effect Model in other study the effect of TA65 at 6 months and 12 months


# KEY INSIGHTS TO INFORM DECISION MAKING
* Both the TA65 and placebo groups experienced significant reductions in hsCRP levels over time, with decreases of approximately 70% at 6 months and 73.57% at 12 months compared to baseline(when patient first arrived)
* At 6 months, people in the TA65 group had results that were about 5.5% higher compared to the placebo group, while at 12 months, their results were about 8% lower compared to the placebo group. However, these differences aren’t meaningful because the data suggests they could simply be due to chance. This shows how important it is to consider both the size of the effect and whether it’s statistically reliable when evaluating a treatment’s impact.
* The statistical model showed a statistically insignificant effect of -6.54% for the TA65 group, meaning that patients who took TA65 had 6.54% lower hsCRP levels than patients who did not take TA65 but this could be due to luck as the pvalue was 0.54.

# DATA CLEANING AND PREPARATION
To facilitate data manipulation and analysis, several R packages were installed and loaded, including dplyr and tidyr for data wrangling, lme4 and lmerTest for fitting linear mixed-effects models, VIM and mice for missing data visualisation and imputation, ggplot2 for data visualisation, broom and broom.mixed for tidy output of model results.
After baseline and hsCRP datasets were loaded, any non-numeric values in hsCRP columns were identified and cleaned using regular expressions to remove unwanted characters and convert the values to numeric types. The baseline dataset was merged with the hsCRP dataset based on the common column ‘Subject ID’ to ensure all baseline information was retained, while hsCRP measurements from multiple time points where added. For longitudinal analysis the new combined dataset was reshaped from wide to long format using the pivot_longer() function in R in other to allow the hsCRP measurements at different time points to be analysed as repeated measures. The time variable was transformed to represent the actual points ( 0, 6 and 12 months). To enhance consistent referencing throughout the analysis column names were standardised to lowercase. Data types where checked and converted to the correct format, subject_id , rand_group, gender, acs_type, and cd8_prop were converted to factors.

**Data in long format can be seen below**
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
The plots consistently demonstrate a decrease in hsCRP levels over time for both treatment groups. The TA65 treatment group exhibits lower hsCRP levels in comparison to the Placebo group particulary in the NSTEMI group and male gender. This observation suggests that TA65 may have a beneficial impact on reducing inflammatory markers over time, potentially influencing clinical outcomes with male patients who are also in the NSTEMI group. Due to the variations in intercepts, it is important to interpret the statistical significance of these findings in conjunction with other analyses to reach definitive conclusions about the efficacy of TA65.

# MODELLING DATA
The lmer package from R. Linear mixed models were fitted using the lmer() function.The models were incrementally built, starting with a simple random intercept model and progressing to interaction models. The models were employed to analyse the clinical trial specifically focusing on hsCRP levels measured at baseline, 6 months, and 12 months.
In our case, the response variable is the log-transformed hsCRP levels (hscrp_log), the fixed effects include treatment group (rand_group), time (timef), and other covariates such as gender, ACS type (acs_type), and CD8 proportion (cd8_prop), while the random effects account for the variability among individual subjects (subject_id).
A Linear Mixed Effect Model was fitted with interaction terms specified for the randomised group and time.
# RESULT AND EVALUATION
This study examined the effect of the TA65 treatment on hsCRP levels, a biomarker of inflammation, over time compared to a placebo group. A linear mixed-effects model was employed to account for repeated measurements and individual variability. The dependent variable, hsCRP, was log-transformed to stabilize variance, allowing the effects to be interpreted as relative changes.

The analysis revealed that hsCRP levels significantly decreased over time in both the placebo and TA65 groups. At six months, hsCRP levels dropped by approximately 70% compared to baseline (
𝑝
<
0.001
p<0.001). By 12 months, the reduction reached approximately 73% (
𝑝
<
0.001
p<0.001). These findings indicate a substantial improvement in inflammatory profiles over time, regardless of treatment group.

At baseline, hsCRP levels in the TA65 group were about 6.5% lower than those in the placebo group. However, this difference was not statistically significant (
𝑝
=
0.758
p=0.758). Over time, the interaction between treatment group and time showed no meaningful differences in the pattern of hsCRP reduction. At six months, the TA65 group’s hsCRP levels were approximately 5.5% higher than those in the placebo group, but this difference was not statistically significant (
𝑝
=
0.862
p=0.862). Similarly, at 12 months, the TA65 group’s hsCRP levels were about 8% higher than the placebo group, with no statistical significance (
𝑝
=
0.788
p=0.788). These results suggest that the TA65 treatment did not significantly alter the natural reduction in hsCRP levels observed in both groups over time.

The model also revealed slight variability in baseline hsCRP levels across individuals, as indicated by the random effects. However, the majority of variability in the data was attributed to residual factors rather than differences between individuals.

These findings carry important implications. While hsCRP levels decreased substantially over time, this reduction appears to reflect general lifestyle or study conditions rather than any specific effect of TA65. The lack of a statistically significant effect of TA65 on hsCRP levels suggests that the treatment may not directly influence inflammation as measured by this biomarker. Any differences observed between the TA65 and placebo groups were small and likely due to random variation.

The analysis has several strengths, including the use of a linear mixed-effects model, which accounts for repeated measurements and individual differences, providing a nuanced understanding of changes over time. However, the study is not without limitations. It focused solely on hsCRP as a marker of inflammation, and additional biomarkers or clinical outcomes might offer a more comprehensive assessment of TA65’s effects. Additionally, the sample size may have limited the ability to detect smaller but potentially meaningful differences.

In conclusion, while hsCRP levels significantly decreased over time, the TA65 treatment did not demonstrate any additional benefits compared to the placebo. Future research should explore larger samples, longer follow-up periods, and alternative markers of inflammation to fully assess the potential effects of TA65 on inflammation and other clinical outcomes.
