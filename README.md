# Hypothesis Testing in Gaussian Graphical Models

This repository contains source code for conducting hypothesis testing in Gaussian graphical models, including a Monte Carlo goodness-of-fit (MC-GoF) test and graphical conditional randomization test (G-CRT) introduced by the following papers:
Xiaotong Lin, Weihao Li, Fangqiao Tian and Dongming Huang. (2025). Goodness-of-Fit Tests via Co-sufficient Sampling. Journal of the Royal Statistical Society Series B, qkaf032.
[paper on [arXiv](https://arxiv.org/abs/2312.01815)].

## Goodness-of-Fit (GoF) Tests

GoF testing is an indispensable part of statistical inference in GGMs. Suppose the rows of the observed data **X** are independent and identically distributed (i.i.d.) samples from some population $P$. Given a graph $G$, the GoF testing problem aims to test the null hypothesis that $$H_0: P \in \mathcal{M}_G$$

We sample exchangeable copies to construct goodness-of-fit tests for Gaussian graphical models and then propose several test statistics. In simulations investigating power properties, our proposed methods excel over other benchmarks especially when the signals are weak and diffuse.

## Conditional Randomization Tests (CRTs)

Conditional independence test is to test the null hypothesis that the response $Y$ is conditionally independent of the covariate variables $X_\mathcal{T}$ given the other covariate variables $X_\mathcal{S}$, where $X_\mathcal{T}$ and $X_\mathcal{S}$ are the sub-vectors of $X$ corresponding to $\mathcal{T}$ and $\mathcal{S}$, respectively. The null hypothesis can be expressed mathematically as:

$$Y \bot X_\mathcal{T} | X_\mathcal{S}$$

Our methodology is based on a combination of the conditional randomization test (CRT) and our algorithm for generating exchangeable copies. Notably, our method imposes less stringent conditions on the distribution of X among existing CRT methods. We discuss different strategies for constructing test statistic functions and demonstrate in simulations that our approach could be much more powerful than alternatives when the sub-vector of X contains many elements.

## Files:

-   Instructions and notations are presented in the R files.

    -   `Function.R`：Functions for all the algorithms and simulations in the paper.

    -   `Example.Rmd`：Examples on implementation of our methods and the benchmark methods for GoF test and CRT test.

    -   `GoF_simulations.R`：Contains the setups and instructions for generating data and implementation of GoF test in our simulation studies. The setups inlcude: 1. Type I error control. 2. Power comparison (band graphs, hub graphs and Erdös-Rényi random graphs). 3. Additional simulations on GoF tests for GGMs (consider an experimental setup in Verzelen and Villers [2009] for a balanced comparison with existing methods). 4. Additional simulations on GoF tests with prior information (using statistic functions PRC-w and ERC-w).

    -   `CRT_simulations.R`：Contains the setups and instructions for generating data and implementation of CRT test in our simulation studies. The setups inlcude: 1. Gaussian linear regression. 2. Logistic regression. 3. Nonlinear regression. 4. Nonlinear binary regression.

    -   `real_data_analysis.R`：Contains the following sections:

        -   Average Daily Precipitation in the United States: 1. Create graph. 2. Simulation Type-I Check. 3. Subsampling.
        -   Dependence of Fund Return: 1. Define X and Y. 2. Estimation by Glasso. 3. GoF Test. 4. Type I error check. 5. CIT Test.
        -   Breast Cancer Relapse: 1. Type I Error. 2. CIT Test. 3. Noise Experiments.

-   Folder `GGM-Test-Version0`: R code for two tests (MC-GoF and G-CRT)

    -   `CSS-GGM.R`：Sampling exchangeable copies
    -   `GGMTesting.R`：Computing p-value
    -   `Gof`: Monte Carlo goodness-of-fit (MC-GoF) test and simulation
    -   `CRT`: Graphical conditional randomization test (G-CRT) and simulation
    -   `Application`: Real-world examples
        -   `US`: Average daily precipitation in the United States
        -   `Stock`: Dependence of fund return
        -   `Breast`: Breast cancer relapse
        
-   Folder `GoF-Additional-Test`: R code for additional simulations of GoF

    -   `Selection_of_L`：Examining the effect of the parameter $L$ in the algorithm *Sampling Exchangeable Copies for a GGM*.
    -   `More_Complicated_Graph_Structures`: Exploring more complicated graph structures.
    
-   Folder `CRT-Additional-Test`: R code for additional simulations of CRT
 
    -   `CRT-Oversized`: Studying the power when a moderately oversized supergraph is used in the G-CRT
    -   `CRT-Univariate`: Examining the case where $X_\mathcal{T}$ is univariate.
    -   `CRT-UniformT`: Studying the performance under violations of the Gaussian assumption.
    
Rendered tutorials demonstrating the usage of the code are available at: [MC-GoF_test_for_GGMs](https://tfq-acd.github.io/MC-GoF_test_for_GGMs/), [GoF_simulation](https://tfq-acd.github.io/GoFsimulation/), [Graphical_conditional_randomization_test](https://tfq-acd.github.io/CRT/), [G-CRT_simulation](https://tfq-acd.github.io/CRTsimulation/).

-   Folder `Tutorials`: The corresponding R code used in the webpage tutorials.
    -   [`MC-GoF_test_for_GGMs.Rmd`](https://tfq-acd.github.io/MC-GoF_test_for_GGMs/): Key steps involved in MC-GoF.
    -   [`GoF_simulation.Rmd`](https://tfq-acd.github.io/GoFsimulation/): MC-GoF simulation example: take band graph for example to evaluate the statistical power and demonstrate the theoretically valid Type I error control of our MC-GoF test with numerical comparisons to established methods.
    -   [`Graphical_conditional_randomization_test.Rmd`](https://tfq-acd.github.io/CRT/): Key steps involved in G-CRT.\
    -   [`G-CRT_simulation.Rmd`](https://tfq-acd.github.io/CRTsimulation/): G-CRT simulation example: take the linear regression model with a low-dimension setting for example to evaluate the statistical power of our proposed G-CRT.\
    -   `CIT_function_for_webpage.R`: Utility functions for MC-GoF tutorials
    -   `GOF_function_for_webpage.R`: Utility functions for G-CRT tutorials

## Reference Paper:

-   Xiaotong Lin, Weihao Li, Fangqiao Tian and Dongming Huang. *Goodness-of-fit tests for high-dimensional Gaussian graphical models via exchangeable sampling.* Journal of the Royal Statistical Society Series B: Statistical Methodology, qkaf032. URL: <https://doi.org/10.1093/jrsssb/qkaf032>.
-   Xiaotong Lin, Jie Xie, Fangqiao Tian and Dongming Huang. *Testing Multivariate Conditional Independence Using Exchangeable Sampling and Sufficient Statistics* Journal of the Royal Statistical Society Series B, under review. URL: <https://arxiv.org/abs/2504.06685>.
