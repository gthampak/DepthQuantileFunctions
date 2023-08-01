# DepthQuantileFunctions
Mathematics Thesis (2022-23) Advised by Gabe Chandler (Pomona College)

Many machine learning and statistical methods rely on distance as a metric for similarity and dissimilarity. However, due to the curse of dimensionality and data sparsity in high dimensional ambient space, distance metrics are not always the best distinguishing measure. Depth Quantile Functions(DQF) map data points to real-valued single variable functions on the [0, 1] interval using ideas of statistical depth. DQFs have been shown to be effective and interpretable as part of algorithms for multiscale geometric feature extraction for high-dimensional and non-Euclidean data with applications in both supervised (classification) and unsupervised (anomaly detection) learning settings. The novelty of this thesis extending the uses of Depth Quantile Functions to clustering. The algorithm outlined and implemented is an user interactive clustering technique.

Depth Quantile Functions (DQFs) are innovated by Gabe Chandler (Pomona College Mathematics and Statistics) and Wolfgang Polonik (UC Davis Statistics). It was been shown to work well in [surpervised learning](https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-2/Multiscale-geometric-feature-extraction-for-high-dimensional-and-non-Euclidean/10.1214/20-AOS1988.full) and [anomaly detection](https://arxiv.org/abs/2201.06682). My work (advised by Gabe Chandler) extends the use of DQFs to unsupervised settings.

---

### Repository Organization

Work done as part of MATH190 and MATH191 (Senior Theses Class) can be found in the folder titled [ThesisWork](https://github.com/gthampak/DepthQuantileFunctions/tree/main/ThesisWork).

Repository has been restructured for organization and readability post thesis work (summer 2023).

Each layer of numbers and characters are used to create a tree-like hierachical structure with topics and subtopics.

- 00 series - Functions to generate simulated datasets
- 01 series - Functions from the dqfAnomaly R package (Chandler). The code can be found [here](https://github.com/GabeChandler/AnomalyDetection).
- 02 series - Functions used to create `dqf.subset` objects, which is used to calculate average depth quantile functions of subsets of data quickly.
- 03 series - Functions that calculate average depth quantile functions from `dqf.subset` object.
- 04 series - The `dqf.clustering` function and its helper functions.
- 10 series - Transformation functions for average depth quantile functions, used to highlight patterns and differences between DQFs.
- 20 series - Functions used to evaluate patterns and differences between DQFs.
- 40 series - Functions that compile comprehensive sets of DQF transformations and metrics for the user.
- 50 series - Evaluating effectiveness of using DQFs for clustering.

---
### Using DQFClustering

Instructions for how to download and use the R package for `dqf.clustering` can be found [here](https://github.com/gthampak/DQFClustering).

---
### TODO:
- More testing with Real Data.
- Make it easier to identify indices of data point of interest from dqf and dqf diagnostic plots.
- Implement adaptivity.
- Explore more ways that dqf.clustering beats single linkage hierarchical clustering.
