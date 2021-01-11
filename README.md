# Data Mining project 2020/2021 [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dmeoli/DataMiningUniPi/master)

This project was developed during the Data Mining course @ 
[Department of Computer Science](https://www.di.unipi.it/en/) @ [University of Pisa](https://www.unipi.it/index.php/english) 
under the supervision of prof. 
[Anna Monreale](http://pages.di.unipi.it/amonreale/).

A project consists in data analysis based on the use of data mining tools.
The project has to be performed by a team of 2/3 students. It has to be performed by
using Python. The guidelines require to address specific tasks and results must be
reported in a unique paper. The total length of this paper must be max 20 pages of text
including figures. The students must deliver both: paper and well commented Python
notebooks.

## Task 1: Data Understanding and Preparation

### Task 1.1: Data Understanding
Explore the dataset with the analytical tools
studied and write a concise “data understanding” report describing data
semantics, assessing data quality, the distribution of the variables and the
pairwise correlations.

Subtasks of DU:
- Data semantics
- Distribution of the variables and statistics
- Assessing data quality (missing values, outliers)
- Variables transformations & generation
- Pairwise correlations and eventual elimination of redundant variables

### Task 1.2: Data Preparation 
Improve the quality of your data and prepare it by extracting new features interesting for describing the customer profile and his
purchasing behavior. These indicators have to be extracted for each customer.
Indicators to be computed are:
- I: the total number of items purchased by a customer during the period of
observation
- Iu: the number of distinct items bought by a customer in the period of
observation
- Imax: the maximum number of items purchased by a customer during a
shopping session
- E: the Shannon entropy on the purchasing behavior of the customer

It is mandatory that each team defines additional indicators leading to the
construction of a customer profile that can lead to an interesting analysis of
customer segmentation.

Once, the set of indicators will be computed the team has to explore the new
features for a statistical analysis (distributions, outliers, visualizations,
correlations).

## Task 2: Clustering Analysis
Based on the customer’s profile explore the dataset using various clustering techniques.
Carefully describe your decisions for each algorithm and which are the advantages
provided by the different approaches.

Subtasks:
- Clustering analysis by K-means:
    - Identification of the best value of k
    - Characterization of the obtained clusters by using both analysis of
the k centroids and comparison of the distribution of variables within
the clusters and that in the whole dataset
    - Evaluation of the clustering results
- Analysis by density-based clustering:
    - Study of the clustering parameters
    - Characterization and interpretation of the obtained clusters
- Analysis by hierarchical clustering:
    - Compare different clustering results got by using different version of
the algorithm
    - Show and discuss different dendrograms using different algorithms
- Final evaluation of the best clustering approach and comparison of the clustering
obtained

Explore the opportunity to use alternative clustering
techniques in the library: https://github.com/annoviko/pyclustering/

## Task 3: Predictive Analysis
Consider the problem of predicting for each customer a label that defines if (s)he is a
*high-spending* customer, *medium-spending* customer or *low-spending* customer.
The students need to:
- Define a customer profile that enables the above customer classification. Please,
reason on the suitability of the customer profile, defined for the clustering
analysis. In case this profile is not suitable for the above prediction problem you
can also change the indicators.
- Compute the label for any customer. Note that, the class to be predicted must be
nominal.
- Perform the predictive analysis comparing the performance of different models
discussing the results and discussing the possible preprocessing that they
applied to the data for managing possible problems identified that can make the
prediction hard. Note that the evaluation should be performed on both training
and test set.
  
## Task 4: Sequential Pattern Mining
Consider the problem of mining frequent sequential patterns. To address the task:
- Model the customer as a sequence of baskets
- Apply the sequential pattern mining algorithm
- Discuss the resulting patterns

Note that you can decide to filter out some customers.

The algorithm presented during the Python lecture on sequential
pattern mining does not consider the opportunity to define temporal constraints.
As an optional point, the students can extend the algorithm and analysis considering one
or more constraints.