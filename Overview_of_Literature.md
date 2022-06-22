# Missing Race and Ethnicity Data in Criminology Research

#criminology

Starting in Spring 2022, Ofer and I are working with Clare Strange (post-doc in criminology at Penn State - Pennsylvania Commision on Sentencing in the Criminal Justice Research Center) on developing an expository paper covering the lack of missing data analysis in criminology research with a specific focus on missing race and ethnicity data.

> In a nutshell, sentencing researchers tend to explore the “race effect” (i.e., the fact that non-White defendants tend to receive significantly harsher punishments than their White counterparts, even when accounting for relevant legal factors) using pooled data and a two-stage process: 1) modeling a binary incarceration outcome (e.g., receiving a jail or prison sentence vs. not) using logistic regression, and then 2) modeling a continuous sentence length outcome using OLS regression (with logged sentence length), or a negative binomial or Poisson model (the list of approaches goes on). Some researchers have implemented more sophisticated methods, but I think we are keeping things simple for the sake of illustration.
> 
> The reason this feels very relevant now is because there is a lot of “talk” about whether the “race effect” still exists (and what its true magnitude is). While most studies do show race effects, some don’t, and that has also fueled this debate. There are hundreds of articles on this topic—the ones that are attached are relatively recent, by high-profile researchers, published in high-quality criminological and sociological journals, and generally use the methods I described above.
\- Clare Strange (email 5/12/2022)

## Questions for Clare Strange

- When describing the data sets, some of the papers describe only retaining "the most serious offense per judicial proceeding". Why retain only these cases? 

## Questions of Interest to Criminology

- Is there a disparity gap between incarceration and arrests? What drives that gap? 
- How does race influence incarceration?
- How is sentencing impacted by intersectionality of race, gender, and age?
- What is the cumulative effect of race and ethnicity throughout the interaction with the criminal justice system from arrest to incarceration decision?

## Types of Data

- Arrest Data
	- Records of crime type and race/ethnicity
	- Uniform Crime Reports
	- Reported by FBI
- Sentencing Data
	- Records of demographics for incarcerated people
	- National Corrections Reporting Program, Florida Sentencing Guidelines Database (Florida Department of Corrections)
	- Binary: Was the offender sentenced?
	- Count/Numeric: How long is the sentence?

## Types of Analyses

- Logistic Regression
	- Univariate Response
	- Multivariate
		- Multiple decision points
		- Multiple categories
- Precision Matching
	- Related to propensity score matching, [matching](https://en.wikipedia.org/wiki/Matching_(statistics), variable by variable matching, etc
	> It is important to note that precision matching and propensity score methods (PSM) are different ways of dealing with the selection bias issue. What precision matching offers is just that, precise matching on key variables. [...] PSM seeks to derive equivalent control and experimental groups by matching cases based on their probability of being placed in the experimental group through the creation of propensity scores [...] based on the covariates available. [...] The important feature of this approach is that cases are matched on the propensity score rather than explicitly on the covariates included in its estimation. [...] In contrast, precision matching results in cases in the control and experimental groups that are exactly the same on all of the covariates included in the matching process. 
	\- Bales and Piquero 2012
-  Tobit Regression
- Two stage analysis
	- First logistic regression on incarceration decision followed by linear regression on the sentence length for only those who are to be incarcerated