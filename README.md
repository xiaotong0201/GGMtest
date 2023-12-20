# Hypothesis Testing in Gaussian Graphical Models
This repository contains source code for conducting hypothesis testing in Gaussian graphical models, including a Monte Carlo goodness-of-fit (MC-GoF) test and graphical conditional randomization test (G-CRT)  [paper on [arXiv](https://arxiv.org/abs/2312.01815)].

- Folder `GGM-Test-Version0`: R code for two tests (MC-GoF and G-CRT)
  - `CSS-GGM.R`：Sampling exchangeable copies
  - `GGMTesting.R`：Computing p-value
  - `Gof`: Monte Carlo goodness-of-fit (MC-GoF) test and simulation
  - `CRT`: Graphical conditional randomization test (G-CRT) and simulation
  - `Application`: Real-world examples
    - `US`: Average daily precipitation in the United States
    - `Stock`: Dependence of fund return
    - `Breast`: Breast cancer relapse

Rendered tutorials demonstrating the usage of the code are available at the [page]. 

- Folder `Tutorials`: The corresponding R code used in the webpage tutorials.
  - `MC-GoF_test_for_GGMs.Rmd`: Key steps involved in MC-GoF[]
  - `GoF_simulation.Rmd`: MC-GoF simulation example: take band graph for example to evaluate the statistical power and demonstrate the theoretically valid Type I error control of our MC-GoF test with numerical comparisons to established methods[]
  - `Graphical_conditional_randomization_test.Rmd`: Key steps involved in G-CRT[]
  - `G-CRT_simulation.Rmd`: G-CRT simulation example: take the linear regression model with a low-dimension setting for example to evaluate the statistical power of our proposed G-CRT[]
  - `CIT_function_for_webpage.R`: Utility functions for MC-GoF tutorials
  - `GOF_function_for_webpage.R`: Utility functions for G-CRT tutorials

## Goodness-of-Fit Tests 

We sample exchangeable copies to construct goodness-of-fit tests for Gaussian graphical models and then propose several test statistics. In simulations investigating power properties, our proposed methods excel over other benchmarks when the signals are weak and diffuse.

## Conditional Randomization Tests

Our methodology is based on a combination of the conditional randomization test [CRT](https://academic.oup.com/jrsssb/article/80/3/551/7048447) and our algorithm for generating exchangeable copies in Section 2. Notably, our method imposes less stringent conditions on the distribution of X among existing CRT methods. We discuss different strategies for constructing test statistic functions and demonstrate in simulations that our approach could be much more powerful than alternatives when the sub-vector of X contains many elements.

## Reference Paper:   
Xiaotong Lin, Fangqiao Tian and Dongming Huang. *Hypothesis
Testing in Gaussian Graphical Models: Novel Goodness-of-Fit Tests
and Conditional Randomization Tests.* 4 december 2023. DOI:
10.48550/arXiv.2312.01815. arXiv: 2312.01815[stat]. URL:
http://arxiv.org/abs/2312.01815.