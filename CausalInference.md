---
name: CausalInference
topic: Causal Inference
maintainer: Pan Zhao, Achim Ahrens, Philipp Bach, Julie Josse, Nick Huntington-Klein, Lucy D'Agostino McGowan
email: pz298@cam.ac.uk
version: 2026-03-04
source: https://github.com/cran-task-views/CausalInference/
---

### Overview

Causal inference can be seen as a subfield of statistical analysis. It
is used in various fields such as econometrics, epidemiology,
educational sciences, etc. With causal inference one addresses questions
about effects of a treatment, intervention, or policy on some target over
a given sample or population. Under certain identifiability and model assumptions, causal
inference can be carried out by fitting simple regression models or
combining several regression models in a specific way as will be
sketched out later. For observational data, additional untestable
assumptions have to be made to (non-parametrically) identify causal
effects.

There are no basic R functions that are direct implementations of standard causal
inference designs, but many methods - more or less complex - are implemented in
different packages on CRAN, which we structure into main topics:

- [Methods for experimental data](#rct)
- [Average treatment effect estimation and other univariate treatment effect estimates](#ate)
- [Heterogeneous treatment effect estimation](#hte)
- [Policy learning and dynamic treatment regimes](#policy)
- [Structural equation models (SEM), do-calculus causal discovery](#dag)
- [Specific types of data](#data)
- [Specific application fields](#applications)

Certain causal inference methods originated in specific fields such as
econometrics or clinical trials and remain most popular therein. In certain
cases, we therefore refer to other task views covering these methods in more
depth. More generally, in this task view we focus on causal analyses with
observational data.

If you think that we missed some important packages in this list, please
contact the maintainers.

### [Methods for randomized controlled trial (RCT) and other experimental data]{#rct}

-   *Construction of experimental designs* is implemented in
    `r pkg("blocksdesign")` (blocks for general factorial
    treatment designs), `r pkg("BCHM")` (Bayesian cluster
    hierarchical model design for multiple subgroup basket trials),
    `r pkg("Boptbd")` (Bayesian optimal block designs under
    linear mixed effects model), `r pkg("seqDesign")`
    (sequential design of randomized two-stage treatment efficacy trials
    with time-to-event endpoints). Many other tools and packages exist
    for designing experiments and clinical trials, we refer to the
    `r view("ExperimentalDesign")` and
    `r view("ClinicalTrials")` CRAN Task Views.
-   *Tests based on pairwise comparisons* are provided in
    `r pkg("BuyseTest")`.
-   *Regression models* where the causal estimand is defined as a
    coefficient of a regression model are implemented in the packages
    `r pkg("allestimates")`.
-   *Analysis methods for RCTs* are provided in `r pkg("experiment")`
    (various statistical methods) and `r pkg("eefAnalytics")`
    (Frequentist and Bayesian multilevel models).
-   *Posterior analysis tools* are implemented in
    `r pkg("cjoint")` (conjoint analysis for survey
    experiments).
-   Design and analysis of *two-stage preference trials* is implemented
    in `r pkg("preference")`.
-   In case of *non-compliance*, `r pkg("rpsftm")` uses
    g-estimation to estimate the causal effect of a treatment in a
    two-armed randomised control trial where non-compliance exists and
    is measured, under an assumption of an accelerated failure time
    model and no unmeasured confounders.
-   A time series causal inference model for RCT under *spillover
    effect* is implemented in `r pkg("SPORTSCausal")`.
-   `r pkg("cosa")` offers a means of performing power analysis for multi-level randomized trials that have continuous outcomes.

### [Average treatment effect estimation and other univariate treatment effect estimates]{#ate}

-   *Regression models* where the causal estimand is a regression
    parameter are implemented in `lm()` and `glm()` from stats, as well as in a number of more specialized packages such as `r pkg("fixest", priority = "core")`,
    `r pkg("estimatr")`, `r pkg("CausalGAM")` (using generalized additive
    models), `r pkg("sampleSelection")` (two-step and maximum
    likelihood estimation of Heckman-type sample selection models),
    `r pkg("BCEE")` (Bayesian causal effect estimation for
    binary or continuous treatment and outcomes),
    `r pkg("borrowr")` (Bayesian PATE estimation for multiple
    exchangeable data sources), `r pkg("causaldrf")` (average
    causal dose response functions), `r pkg("hdm")` (efficient
    estimators with uniformly valid confidence intervals, it assumes
    approximately sparse models for high-dimensional settings). Estimation in
    *fixed effects designs* is possible through `r pkg("fixest")` (linear and
    generalized linear fixed effects models and combined with instrumental
    variables), `r pkg("plm")` (for panel data), and `r pkg("alpaca")` (for
    high-dimensional k-way fixed effects).
-   *G-computation* and other *conditional outcome regression* based methods
    are supported in the packages `r pkg("gfoRmula")` (also
    for time-varying treatment and confounding),
    `r pkg("EffectLiteR")` (based on structural equation
    modeling), `r pkg("switchSelection")` (maximum likelihood or two-step
    estimation of endogenous switching regression models), and
    `r pkg("riskRegression", priority = "core")` (for survival
    outcomes with or without competing risks). For parametric models, g-computation is the same as estimating average marginal effects, which can be achieved using `r pkg("margins")`, `r pkg("marginaleffects")`, `r pkg("modelbased")`, and `r pkg("stdReg")`.
-   *Matching* methods are implemented in `r pkg("MatchIt", priority = "core")`, which provides wrappers for a number of popular methods including propensity score matching and subclassification, (coarsened) exact matching, full matching, and cardinality matching; more specialized matching methods are implemented in some of the packages below, some of which MatchIt depends on. `r pkg("MatchThem")` provides a wrapper for MatchIt with multiply-imputed data. `r pkg("Matching", priority = "core")` performs nearest neighbor and genetic matching and implements Abadie and Imbens-style matching imputation estimators. `r pkg("optmatch")` performs optimal matching using network flows; several other packages rely on the same infrastructure, including `r pkg("DiPs")` (near-fine matching with directional penalties), `r pkg("matchMulti")` (optimal matching for clustered data), `r pkg("rcbalance")` and `r pkg("rcbsubset")` (optimal matching for refined balance), and `r pkg("approxmatch")` (near-optimal matching for multi-category treatments). Other packages include `r pkg("cem")` (coarsened exact matching), `r pkg("stratamatch")` (matching and stratification in large datasets), `r pkg("FLAME")` (almost-matching-exactly via learned weighted Hamming distance), `r pkg("PanelMatch")` (matching with time-series cross-sectional data), `r pkg("kbal")` (kernel balancing), `r pkg("vecmatch")` for vector matching with multiple treatments, and `r pkg("CausalGPS")` (generalized propensity score matching for continuous treatments).
-   *Inverse propensity weighting* (IPW, also known as inverse probability of treatment weighting, IPTW) methods are implemented in `r pkg("WeightIt", priority = "core")`, which provides implementations and wrappers for several popular weighting methods for binary, multi-category, continuous, and longitudinal treatments. `r pkg("MatchThem")` provides a wrapper for WeightIt with multiply-imputed data. `r pkg("PSweight", priority = "core")` offers propensity score weighting and uncertainty estimation using M-estimation. `r pkg("inferference")` offers weighting methods in the context of interference. Several packages offer specialized methods of estimating balancing weights for various treatment types, which may or may not involve a propensity score: `r pkg("CBPS")` (generalized method of moments-based propensity score estimation for binary, multi-category, continuous, and longitudinal treatments), `r pkg("twang")` and `r pkg("twangContinuous")` (propensity score weighting using gradient boosting machines for binary, multi-category, continuous, and longitudinal treatments), `r pkg("sbw")` and `r pkg("optweight")` (optimization-based weights using quadratic programming), and `r pkg("ebal")` (entropy balancing). `r pkg("mvGPS")` estimates weights for multivariate treatments using WeightIt's infrastructure. *Matching-adjusted indirect comparison*, a relative of propensity score weighting when unit-level data is only available for some groups, is available in `r pkg("maicChecks")` and `r pkg("optweight")` (using the `optweight.svy()` function).
-   *Doubly robust methods* involve both a treatment and outcome model. Augmented IPW (AIPW) is implemented in `r pkg("AIPW")`, `r pkg("PSweight")`, `r pkg("DoubleML")`, `r pkg("grf")` (functions `causal_forest` followed by `average_causal_effect`), and `r pkg("causalweight")`. Targeted maximum likelihood estimation (TMLE, also known as targeted minimum loss-based estimation) is available in `r pkg("drtmle")`, `r pkg("tmle", priority = "core")`, `r pkg("ctmle")` (for TMLE with variable selection), `r pkg("ltmle")` (for longitudinal data), and `r pkg("AIPW")`.
-   *Difference in differences* methods are implemented in
    `r pkg("DRDID")` (doubly robust estimators with two
    choices for nuisance function estimation),
    `r pkg("bacondecomp")` (using the Goodman-Bacon
    decomposition to allow for variation in treatment timing),
    `r pkg("did")` (for cases with more than two periods and
    with variation in treatment timing), `r pkg("fixest")` (Sun & Abraham
    estimator), and in `r pkg("qte")`.
-   *Quantile treatment effects* can be estimated using the
    `r pkg("qte")`, `r pkg("Counterfactual")` and
    `r pkg("grf")` packages.
-   *Odds ratio* estimation and power calculation for the *Trend in
    Trend* model is implemented in `r pkg("TrendInTrend")`.
-   *Synthetic control* methods are implemented in
    `r pkg("Synth")` (using a group method for comparative
    case studies),`r pkg("microsynth")` (for micro- and
    meso-level data), and `r pkg("gsynth")` (extension to
    multiple treated units and variable treatment periods).
    `r pkg("tidysynth")` offers an easy-to-use syntax for using synthetic
    control methods. `r pkg("scpi")` implements multiple synthetic control estimators using lasso, ridge, simplex, and linear constraints, and constructs prediction intervals.
-   *Instrumental variable* methods are implemented in
    `r pkg("ivreg")`, `r pkg("ivmodel")`,
    `r pkg("bpbounds")` (nonparametric bounds on ATE),
    `r pkg("grf")`, `r pkg("fixest")`, `r pkg("estimatr")`, and
    `r pkg("DoubleML")` (function `DoubleMLIIVM`).
    `r pkg("ivmte")` provides a choice-theoretic
    interpretation to IV models using *Marginal Treatment Effects* to
    extrapolate from the compliers to estimate treatment effects for
    other subpopulations. `r pkg("LARF")` uses
    Local Average Response Functions for IV estimation of treatment
    effects with binary endogenous treatment and instrument.
    `r pkg("ivdesc")` gives descriptive statistics for the 
    complier, never-taker and always-taker subpopulations. 
    More details and a longer list of packages for
    IV methods can be found in `r view("Econometrics", "Instrumental variables")`
    in the `r view("Econometrics")` task view.
-   *Mediation analysis* can be performed with `r pkg("cit")`
    (likelihood-based tests), `r pkg("MultisiteMediation")`
    (multisite trials), `r pkg("DirectEffects")` (controlled
    direct effect when fixing a potential mediator to a specific value),
    `r pkg("medflex")` (natural effect models). `r pkg("causalweight")`, `r pkg("rmpw")`, and `r pkg("twangMediation")`
    implement weighted estimators for mediation.
    `r pkg("mediation", priority = "core")` and `r pkg("cfdecomp")`
    implement identification, inference and both
    `r pkg("mediation")` nad `r pkg("mediationsens")` provide
    sensitivity analysis for causal mediation effects.
    `r pkg("bama")` performs Bayesian mediation analysis with high-dimensional mediators. 
    `r pkg("paths")` uses an imputation approach to estimate
    path-specific causal effects along with a set of bias formulas for
    conducting sensitivity analysis. `r pkg("regmedint")`
    implements regression-based analysis with a treatment-mediator
    interaction term. `r pkg("bmem")` provides several different
    methods for mediation analysis in the case of missing data (listwise/pairwise 
    deletion, multiple imputation, two stage maximum likelihood) and power analysis 
    for mediation analysis. `r pkg("bmemLavaan")` mediation analysis with missing data and non-normal data.
-   Under *interference,* causal effect estimation can be achieved using
    `r pkg("inferference")` by inverse-probability weighted
    (IPW) estimators.
-   Diagnostics and visualization for *Multiplicative Interaction
    Models* are implemented in `r pkg("interflex")`.
-   `r pkg("InvariantCausalPrediction")` provides confidence
    intervals for causal effects, using data collected in different
    experimental or environmental conditions (with hidden variables),
    extensions to nonlinear models are implemented in
    `r pkg("nonlinearICP")`.
-   *Regression discontinuity design* (RDD) methods are implemented in
    `r pkg("rdrobust")` (offering robust confidence interval construction and
    bandwidth selection). A more detailed curated list of packages for
    RDD methods can be found in `r view("Econometrics", "Regression discontinuity design")`
    in the `r view("Econometrics")` task view.

In addition, `r pkg("causalsens")`, `r pkg("OVtool")`,
`r pkg("dstat")`, `r pkg("sensemakr")` and `r pkg("EValue")` provide functions for *sensitivity analyses*
(for unmeasured confounding, selection bias, measurement error),
and `r pkg("ui")` implements functions to derive uncertainty
intervals and conduct sensitivity analysis for missing data and
unobserved confounding. `r pkg("cobalt", priority = "core")` and `r pkg("tableone")`
generate balance tables and plots before and after covariate balancing, while `r pkg("BalanceCheck")` offers tests for balance between groups. `r pkg("Counternull")` provides methods for effect testing using randomization inference. `r pkg("WhatIf")` offers methods to assess overlap and extrapolation.

### [Heterogeneous treatment effect estimation]{#hte}

Some of the above mentioned packages can also be used for heterogeneous
treatment effect (HTE) estimation.

-   *Bayesian approaches* for heterogeneous causal effect estimation are
    available in `r pkg("bartCause")` (based on Bayesian Additive
    Regression Trees) and `r pkg("stochtree")` (which implements the Bayesian
    Causal Forest method, based on Bayesian Additive Regression Trees).
-   *Fisherian approaches* for an omnibus test of heterogeneity and decomposition of overall treatment effect heterogeneity into a systematic component explained by covariates and an idiosyncratic component is implemented in `r pkg("hettx")`.
-   *Efficacious treatment or population subset selection* exploiting
    treatment effect heterogeneity is implemented in
    `r pkg("FindIt")` and `r pkg("grf")`. This
    latter package supports missing covariate values using the Missing
    Incorporated in Attributes approach. Additionally the package
    `r pkg("subdetect")` provides a test for the existence of
    a subgroup with enhanced treatment effect.
-   Other approaches for personalized causal predictions are provided by
    `r pkg("EffectTreat")` (exploiting correlation-based
    expressions), and for randomized data by
    `r pkg("evalITR")` (it additionally allows for defining
    budget constraints) by `r pkg("SortedEffects")`
    (estimation and inference methods for sorted causal effects and
    classification analysis), and by `r pkg("CRE")` (exploration of heterogeneity patterns using an ensemble-of-trees approach).
-   `r pkg("stepp")` provides diagnostic plots to explore
    *treatment-covariate interactions* for survival or generalized
    linear models, applicable for continuous, binomial and count data
    arising from two or more treatment arms of a clinical trial.
-   There are several approaches to HTE estimation that come out of machine learning. Causal forests, as covered in `r pkg("grf")` are one. Another, using an ensemble-of-trees approach to discover heterogeneity, is in `r pkg("CRE")`.

### [Policy learning and dynamic treatment regimes]{#policy}

- *Direct methods* perform a direct classification task: outcome-weighted learning `r pkg("DTRlearn2")`; efficient augmentation/relaxation learning (EARL), residual learning, weighted learning, value-search methods based on Augmented Inverse Probability Weighted Estimators and Inverse Probability Weighted Estimators `r pkg("DynTxRegime")`.
- *Indirect methods*, which proceed through nuisance estimation to perform Q-learning `r pkg("DTRlearn2")`, `r pkg("DynTxRegime")`, `r pkg("DTRreg")`, Interactive Q-Learning `r pkg("DynTxRegime")` or double robust Q-Learning.
r pkg("polle") provides a unified framework for learning and evaluating finite stage policies based on observational data with methods such as doubly robust restricted Q-learning, policy tree learning, and outcome weighted learning. Flexible machine learning methods can be used to estimate the nuisance components and valid inference for the policy value is ensured via cross-fitting. The package wraps and extends some functionalities from other packages  `r pkg("DynTxRegime") `,  `r pkg("policytree") `,  `r pkg("grf") `,  `r pkg("DTRlearn2")`.
- *For sequential, multiple assignment, randomized trials (SMART)*
Set of tools for determining the necessary sample size in order to identify the optimal DTR `r pkg("smartsizer")`. Estimators for general K-stage DTRs from SMARTs `r pkg("DTRlearn2")`.
- *With variable selection*
Outcome weighted learning with variable selection via penalization `r pkg("DTRlearn2")`. Estimation of individualized treatment rules from observational and randomized data with options for variable-selection and gradient boosting based estimation, and for outcome model augmentation for continuous, binary, count, and time-to-event outcomes `r pkg("personalized")`. Penalized regression that can simultaneously estimate the optimal treatment strategy and identify important variables for either censored or uncensored continuous response `r pkg("OTRselect")`.
- *Other methods using quantiles*:
Marginal quantile, marginal mean, and mean absolute difference and doubly-robust methods for quantile-optimal treatment regime `r pkg("quantoptr")`, estimation of mean- and quantile-optimal treatment regimes from censored data `r pkg("QTOCen")`.
- *Other approaches*
Learn optimal policies via doubly robust empirical welfare maximization over trees `r pkg("policytree")`. Doubly-robust causal effect estimated for modified treatment policies, dynamic treatment regimes and static interventions `r pkg("lmtp")`. G-estimation, dynamic weighted OLS, several variance estimation approaches with possibilities to handle survival outcomes and continuous treatment variables `r pkg("DTRreg")`. Methods using kernel smoothing to examine optimal linear regimes `r pkg("DTRKernSmooth")`. Single-Index Models with Multiple-Links for, respectively, experimental and observational data `r pkg("simml")`, `r pkg("simsl")`


### [Structural equation models (SEM), do-calculus causal discovery]{#dag}

-   *Identifiability* is addressed by `r pkg("causaleffect")`
    and `r pkg("dosearch")` providing algorithms to decide
    whether a causal effect is identifiable (non-parametric
    identifiability) and by `r pkg("CausalQueries")` that
    calculates arbitrary estimands for a given causal model. 
    `r pkg("SEMID")` implements SEM-based routines to check identifiability or non-identifiability of linear SEMs.
    `r pkg("causaloptim")` provides tight bounds for a user
    defined DAG, query and constraints using a symbolic linear
    optimizer.
-   *Causal structure learning* is possible with functions from
    `r pkg("pcalg", priority = "core")`: PC, for observational
    data without hidden variables, FCI and RFCI, for observational data
    with hidden variables, and GIES, for a mix of observational and
    interventional data without hidden variables;
    `r pkg("pcalg")` also allows to do *causal
    inference using graphical models* (the IDA algorithm, the Generalized
    Backdoor Criterion - GBC, the Generalized Adjustment Criterion - GAC).
    Incorporating background knowledge is also possible. Many algorithms and
    methods for general and specific graphical models exist, we refer to the
    `r view("GraphicalModels")` and `r view("Psychometrics")` CRAN Task Views
    for a comprehensive overview.
-   *Estimation of causal effects* is possible in
    `r pkg("CIEE")` using estimating equations derived from a
    DAG and in `r pkg("InvariantCausalPrediction")` using
    adjustment sets derived from conditional independence tests that
    leverage causal invariances across environments.
-   *Causal networks estimation* is implemented in
    `r pkg("CompareCausalNetworks")`.
-   `r pkg("generalCorr")` computes generalized correlations,
    partial correlations and *plausible causal paths*.

In addition, `r pkg("dagitty", priority = "core")` provides
methods to define different types of graphical models (cpdags, pdag,
ect.) and to identify adjustment sets (a web-based graphical
environment is also available: [DAGitty](http://dagitty.net)). The package `r pkg("ggdag")` produces plots of these causal diagrams from within R. The package package `r pkg("causalHyperGraph")` instead provides plots of causal hypergraphs for configurational comparative methods, the .


### [Specific types of data]{#data}

-   *Longitudinal data / time series and censored data*: Causal effect
    estimation for time series is implemented in
    `r pkg("CausalImpact")` (using a Bayesian approach) and
    `r pkg("CausalMBSTS")` (for multivariate responses).
-   *GWAS and SNPs*: `r pkg("CKAT")` implements kernel based
    methods to jointly test genetic main effect and gene-treatment
    interaction effects for a set of SNPs.
-   *Example data sets* to run frequent example problems from causal
    inference textbooks are accessible through the
    `r pkg("causaldata")` package.
-   Weighted, two-mode, and longitudinal networks analysis is
    implemented in `r pkg("tnet")`
-   Latent treatment effect estimation in text corpora is implemented in `r pkg("texteffect")`.  

### [Specific application fields]{#applications}

-   Behavior change sciences use specialized analyses and visualization
    tools implemented in `r pkg("behaviorchange")`.
-   Evaluation of biomarkers and estimation of treatment-biomarker
    effects can be done using tools from `r pkg("bhm")` (for
    biomarker-treatment effects).
-   Qualitative Comparative Analysis type methods are implemented in
    `r pkg("cna")`.
-   *Mendelian randomization methods* used to examine causal effects
    related to certain genes are implemented in
    `r pkg("MendelianRandomization")` and
    `r pkg("MRPC")` (PC algorithm with the principle of
    Mendelian Randomization).
-   Causal inference approaches in genetic systems exploit quantitative
    trait loci (QTL) genotypes to infer causal relationships among
    phenotypes: functions to simultaneously infer causal graphs and
    genetic architecture (acyclic and cyclic) are implemented in
    `r pkg("qtlnet")`. `r pkg("cophescan")` takes a Bayesian approach to causal inference in genetic systems, addressing gene/trait confounding that arises due to linkage disequilibrium.
-   `r pkg("tools4uplift")` uplift modeling aims at predicting
    the causal effect of an action such as a marketing campaign on a
    particular individual.
-   *Coincidence analysis* through configurational comparative methods
    is provided by `r pkg("cna")`.



### Links
-   [DAGitty: draw and analyze causal diagrams](http://dagitty.net)
