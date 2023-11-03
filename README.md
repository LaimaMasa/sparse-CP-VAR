# Sparse change point vector autoregression, or sparse CP-VAR(p)

## Abstract

When modeling long economic and financial time series, several problems of parameter estimation arise, often due to the ignorance of structural breaks. To overcome these problems, change point methods can be incorporated into the models.

In this thesis, we implement the sparse change point vector autoregressive model (CP-VAR(p)) using the Just Another Gibbs Sampler (JAGS) program. This model is known for its rich parameterization, which can lead to a large number of parameters when dealing with multiple regimes. To mitigate this problem, [Dufays et al., 2021] proposed shrinkage priors to shrink the parameters of irrelevant regimes to zero. The aim of this study is to use the sparse CP-VAR(p) model to detect the number of change points in mean parameters in time series.

To evaluate the robustness of the model in detecting change points, we first perform a simulation study on a bivariate CP-VAR(1) data generation process using JAGS. We can identify two out of three change points in a simulated time series of T = 500. However, our model had slightly poorer change point detection performance in the shorter time series of T =120. Also, our model can identify the absence of the change points in most cases. Unfortunately, the computational power of our machine limits the scope of a more extensive simulation study.

We then apply the model to two financial and economic datasets that include monthly observations between 2013 and 2022. The first dataset includes log returns of stylized stock indices (value and growth) and changes in the Consumer Price Index (CPI), while the second dataset includes variables from energy markets and two macro variables. In both experiments, we set the lag to p=1 and find that no change point is detected with the methodology used. However, we identify some interesting and significant variables.

In summary, our study demonstrates the effectiveness of our sparse CP-VAR(p) model in detecting change points on simulated data. We also show the capabilities of our model in identifying the change points and significant parameter values in time series of financial and economic data.
