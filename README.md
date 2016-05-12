# SpatioTemporalCovariance
Reproducible research code for covariance estimation methods regularized using Kronecker product structure. Author: Kristjan Greenewald.


Robust Kronecker Product PCA for Spatio-Temporal Covariance Estimation
===

Authors: Kristjan Greenewald, Alfred O. Hero

Abstract
---

Kronecker PCA involves the use of a space vs. time Kronecker product decomposition to estimate spatio-temporal covariances. In this work the addition of a sparse correction factor is considered, which corresponds to a model of the covariance as a sum of Kronecker products of low (separation) rank and a sparse matrix. This sparse correction extends the diagonally corrected Kronecker PCA of [Greenewald et al 2013, 2014] to allow for sparse unstructured “outliers” anywhere in the covariance matrix, e.g. arising from variables or correlations that do not fit the Kronecker model well, or from sources such as sensor noise or sensor failure. We introduce a robust PCA-based algorithm to estimate the covariance under this model, extending the rearranged nuclear norm penalized LS Kronecker PCA approaches of [Greenewald et al 2014, Tsiligkaridis et al 2013]. An extension to Toeplitz temporal factors is also provided, producing a parameter reduction for temporally stationary measurement modeling. High dimensional MSE performance bounds are given for these extensions. Finally, the proposed extension of KronPCA is evaluated on both simulated and real data coming from yeast cell cycle experiments. This establishes the practical utility of robust Kronecker PCA in biological and other applications.


Kronecker Sum Decompostions of Space-Time Data
===

Authors: Kristjan Greenewald, Theodoros Tsiligkaridis, Alfred O. Hero III

Abstract
---

In this paper we consider the use of the space vs. time Kronecker product decomposition in the estimation of covariance matrices for spatio-temporal data. This decomposition imposes lower dimensional structure on the estimated covariance matrix, thus reducing the number of samples required for estimation. To allow a smooth tradeoff between the reduction in the number of parameters (to reduce estimation variance) and the accuracy of the covariance approximation (affecting estimation bias), we introduce a diagonally loaded modification of the sum-of-kronecker products representation in [Tsiligkaridis et al. 2013]. We derive an asymptotic Cramer-Rao bound (CRB) on the minimum attainable mean squared predictor coefficient estimation error for unbiased estimators of Kronecker structured covariance matrices. We illustrate the accuracy of the diagonally loaded Kronecker sum decomposition by applying it to the prediction of human activity video. 


