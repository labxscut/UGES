---
title: A unified genetic and epigenetic model to predict breast cancer
  intrinsic subtypes using large DNA-level multi-omics data and
  hierarchical learning
---

Xintong Chang, Jiemin Xie, Hongyu Duan, Keyi Li, Xuemei Liu, Yunhui
Xiong, Xiangqi Bai, Kaida Ning, Li C. Xia*

*Correspondence author: Li C. Xia (<lcxia@scut.edu.cn>)

***Abstract*--- Breast cancer subtyping presents a significant clinical
and scientific challenge. The prevalent expression-based Prediction
Analysis of the Microarray of 50 genes (PAM50) system and its
Immunohistochemistry (IHC) surrogate tests showed substantial
inconsistencies and did not apply to the rapidly progressing circulating
tumor DNA screenings. We developed Unified Genetic and Epigenetic
Subtyping (UGES), a new intrinsic subtype classifier, by integrating
large-scale DNA-level omics data with a hierarchy learning algorithm.
Our benchmarks showed that both multi-step hierarchical learning and
using all DNA-level alteration data are crucial, improving the overall
AUC score by over 8.3% compared to the one-step multi-classification
method. Based on these insights, we developed UGES, a three-step
classifier based on 50831 DNA features of 2065 samples, including
mutations, copy number** **aberrations, and methylations. UGES achieved
an overall AUC score of 0.963 and greatly improved the clinical
stratification of real-world patients, as each subtype strata's survival
difference became statistically more significant, P=9.7e-55 (UGES) vs.
2.2e-47 (PAM50). Finally, UGES identified 52 subtype-specific DNA
biomarkers that can be targeted in early screening technology to expand
the time window for precision care. The UGES code is freely available at
<https://github.com/labxscut/UGES>.**

***Index Terms*--- Breast cancer subtype classifier; Hierarchical
learning; Subtype-specific biomarkers; Multi-omics big-data; Cancer
early screening.**

# I. INTRODUCTION

BREAST cancer is a malignant, complex, and
highly heterogeneous tumor that overwhelmingly impacts women. Precise
identification of the cancer's subtypes is crucial to optimize treatment
and improve outcome for patients \[1\]. However, breast cancer subtyping
remains a difficult scientific and clinical challenge \[2-4\]. The most
prevalent subtyping system, termed PAM50, can discriminate among four
main intrinsic subtypes: Basal-like (**Basal**), Her2-enriched
(**Her2**), Luminal A (**LumA**), and Luminal B (**LumB**). This gene
expression-based system was first established in the Prediction Analysis
of Microarray of 189 patients based on the expression of 50 signature
genes (i.e., **PAM50** subtyping) \[5\]. Clinically, breast cancer
intrinsic subtyping was surrogated to Immunohistochemistry tests
(**IHC** subtyping) specifically targeting estrogen receptor (ER),
progesterone receptor (PR), human epidermal growth factor 2 (Her2), and
Ki67 proteins \[6, 7\].

However, our analysis, along with abundant external evidence, revealed
pervasive and significant inconsistencies between PAM50 and IHC subtypes
\[8-10\]. The discrepancies were observed in both the TCGA (U.S.) and
METABRIC (U.K.) populations. In the METABRIC dataset, 599 out of 1,086
samples (55.16%) showed discordant PAM50 and IHC subtypes (see
**Supplementary Table S1**). Independent studies support these findings.
For instance, a study conducted at the Tri-Service General Hospital in
Taiwan analyzed 372 breast cancer patients. It found that the IHC
subtyping misclassified up to 61.11% of PAM50 Her2-enriched patients as
Luminal B \[9\]. Another study from the Samsung Medical Center (n=607)
reported classification discrepancies in 38.4% of patients. The
misclassification between Luminal-A and Luminal-B types was associated
with significantly poorer overall survival (p\<0.001) \[10\]. These
findings highlighted the dire clinical consequence of such
misclassifications.

The inconsistencies between PAM50 and IHC subtyping arise from several
factors. First, expression-based subtyping is inherently variable due to
the high variability of mRNA expression, which is often influenced by
cellular, temporal, and spatial noises \[11-13\]. As a result, PAM50's
classification based on the expression of only 50 genes frequently
deviates from the actual intrinsic subtypes, as proven by numerous
studies \[14, 15\]. Second, the IHC surrogate for PAM50 is
systematically biased. While using only four biomarkers is convenient in
clinical settings, it overlooks rich genome-wide patterns essential for
accurate subtyping \[16-18\]. In fact, intrinsic subtypes are inherently
determined by a broader range of genetic and epigenetic factors, *i.e.*,
mutations, copy number aberrations, and methylation alterations
\[19-21\]. Last but not least, both PAM50 and IHC subtyping
underestimated inter-tumoral heterogeneity. Developed in the early 2000s
with only small sample sizes, the methods missed important features that
only become identifiable with large-scale omics data \[5, 22\].

To overcome these challenges and improve subtyping accuracy and
consistency, we proposed to identify novel genetic and epigenetic
determinants of breast cancer intrinsic subtypes using large-scale data.
Based on the identified features, we developed a DNA-level classifier
(**Fig. 1a**). Unlike mRNA and protein changes, DNA alterations are more
fundamental, discrete, robust, and easier to ascertain \[23\]. According
to the central dogma, genetic information flows from DNA to RNA and then
to protein, dictating that DNA-level intrinsic patterns would determine
expression and proteomic phenotypes. Additionally, DNA molecules are
structurally stable and less prone to transcriptional and translational
noises, making DNA alterations more reliable for subtyping \[24\].
Notably, early breast cancer detection using the genomic and epigenomic
sequencing of circulating tumor DNA, either targeted or genome-wide, is
picking up the pace \[25\]. By integrating these sequencing advances,
our DNA-level classifier enables breast cancer subtyping at an early
stage, even upon diagnosis, thus significantly extending the time window
for precision care.

The integration of multi-omics profiles for cancer subtyping has been
widely explored in studies of various cancers \[26-28\]. DNA
alterations, such as mutations \[29\], copy number aberrations
(**CNAs**) \[30, 31\], and methylation changes \[32, 33\], have been
used individually or in combination with transcriptomic features for
subtyping. For example, Curtis et al. proposed a system of
transcriptional processes influenced by methylation and CNAs,
demonstrating that CNAs are key drivers of inter-tumor heterogeneity in
breast cancer \[34\]. Islam et al. developed a deep neural network
classifier based on expression and CNA data \[35\], while Lin et al.
created a classifier incorporating expression, CNAs, and methylation
data \[36\]. Choi et al. proposed moBRCA-net, which integrates
expression, methylation, and miRNA data using attention mechanisms for
breast cancer subtyping \[37\]. Most recently, Rajpal et al. introduced
XAI-MethylMarke, a method combining deep learning and explainable AI to
identify a set of 52 DNA methylation biomarkers for breast cancer
subtype classification \[38\].

Improving upon these studies, we developed a DNA-level classifier by
leveraging the two largest multi-omics datasets -- **The Cancer Genome
Atlas (TCGA)** and **METABRIC** cohorts (n=2,065). This classifier,
termed **UGES** (Unified Genetic and Epigenetic Subtyping for Breast
Cancer), is shown in **Fig. 1c**. Methodologically, instead of using the
conventional one-step multi-class learning approach, we introduced a
novel hierarchical learning algorithm designed to reflect the underlying
tumor evolution process. Our stepwise hierarchical classification method
identifies intrinsic subtypes incrementally, effectively capturing
nuanced similarities and differences between subtypes. Compared to the
one-step multi-class strategy, our hierarchical method better addresses
class confusion by optimizing classification tasks at all sublevels.
This significantly reduces overall misclassification rates.

This approach proved crucial, improving overall AUC by at least 8.3%
compared to the state-of-the-art one-step methods. Additionally, when
evaluated against two other competitive hierarchical approaches, the
**UGES** algorithm demonstrated superior performance, achieving the best
precision-recall performance (AUC=0.963) and the highest prognostic
prediction power. As an application, we used the UGES-classified
subtypes from the **TCGA** and **METABRIC** datasets to conduct a
survival analysis. The results demonstrated its clinical relevance,
showing significantly improved prediction for patients compared to
PAM50. In addition, we implemented an effective Lasso-Logistic model for
supervised learning, which soft-selected significant DNA alteration.
The follow-up differential analysis identified 52 key DNA alterations
that potentially define breast cancer intrinsic subtypes.

# II. Materials and Methods

## A. Data collection and pre-processing 

The Cancer Genome Atlas (**TCGA**) \[39\] and Molecular Taxonomy of
Breast Cancer International Consortium (**METABRIC**) \[34\] datasets,
including mutation, CNA, methylation, transcription, and clinical data,
were downloaded from the cBioPortal Cancer Genomics Platform
(<https://www.cbioportal.org>) \[40\] and the TCGA data portal
(<https://portal.gdc.cancer.gov/>).

We summarized all raw DNA alterations to the gene level. For mutation,
we mapped the downloaded mutation profiles to the gene level, resulting
in a binary matrix where the $i,j$-th cell is an indicator of the $i$-th
gene that was mutated at least once in the $j$-th patient. For CNA, we
used the downloaded absolute copy number values, as estimated by the
TCGA consortium and mapped to the gene level by cBioPortal. For
methylation, we mapped the raw methylation array profiles
(HumanMethylation450 BeadChip, Illumina) to the gene level using the
*FDb.InfiniumMethylation.hg19* package in *R*. In the case of
many-to-one mapping between the CpG sites and a gene, we used the
average methylation level for the gene.

We intersected DNA features from both datasets, retaining only shared
methylation features. DNA alteration features with over 50% missing
values were excluded, and missing values in the remaining features were
imputed using the 10-nearest neighbor method with the *impute.knn()*
function from the *impute* R package \[41, 42\]. This resulted in the
inclusion of 50,831 gene-level DNA features, including 16,770 mutations,
25,594 CNAs, and 8,467 methylation features.

## B. Lasso-Logistic regression

We applied the Lasso-Logistic regression model using the *glmnet*
package in R \[43\] with 10-fold cross-validation for classification.
The Least Absolute Shrinkage and Selection Operator (**Lasso**) handles
high-dimensional, multi-modal data while avoiding overfitting by
regularizing non-informative features to zero \[44\]. The model's single
hyperparameter λ controls the regularization level. To determine the
optimal λ value, we used the *cv.glmnet()* function and selected the
model with the highest accuracy.

To describe this multinomial model, we used the following annotations.
Suppose the response variable has $K$ levels
$G = \{ 1,\ 2,\ \ldots,\ K\}$ and each sample $x_{i}$ in our study has
m = 50,831 features, i.e.,
$x_{i} = \{ x_{i,1},\ x_{i,2},\ldots,x_{i,m}\}$, then we model:

$$\Pr\left( G = k \middle| X = x \right) = \frac{e^{\beta_{0k} + \beta_{k}^{T}x}}{\sum_{l = 1}^{K}e^{\beta_{0l} + \beta_{l}^{T}x}}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (1)$$

Let $Y$ be the $N \times K$ indicator response matrix, with elements
$y_{il} = I(g_{i} = l)$. Then the Lasso elastic net penalized negative
log-likelihood function becomes:

$$
l(\{\beta_{0k}, \beta_k\}_1^K) =-[\frac{1}{N}\sum_{i = 1}^{N}(\sum_{i = 1}^{K}y_{il}(\beta_{0k}+x_i^T\beta_k)-log(\sum_{i = 1}^{K}e^\beta_{0l}+x_i^T\beta_l))]+\lambda\sum_{j = 1}^{p}|\beta_j| \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (2)$$


where $\beta$ is a $p \times K$ matrix of coefficients, $\beta_{k}$ is
the $k$-th column (for outcome category $k$), $\beta_{j}$ the $j$-th row
(vector of $K$ coefficients for variable $j$), $\lambda$ is the penalty
parameter, and *N* is the total number of samples.

## C. Hierarchical learning algorithm

We extended the hierarchical learning algorithm proposed by Bengio et
al. \[45\] to construct the optimal multi-step classification training
framework (see **Algorithm 1**).

-------------------------------------------------------------------------
 Algorithm 1: the UGES training framework                              
-------------------------------------------------------------------------
**Input:**

$$
\mathbf{X_{\text{train}}} = (X_1, X_2, \dots, X_n),\quad
\mathbf{Y_{\text{train}}} = (Y_1, Y_2, \dots, Y_n),\quad
k = 4
$$

The subtype hierarchy: T=\{\}.

**Repeat until k=1:**

 1\. Train by cross-validation **k** One-vs-Rest classifiers $${\bar{f}_1}^{(k)}, \ldots, {\bar{f}_k}^{(k)}$$ using Lasso-Logistic regression.                                      

2\. Compute the pairwise classification error, the confusion matrix $$\bar{C}_{ij}$$ and the symmetrized confusion matrix **A**:

i. ${\bar{C}_{ij} = |{(x,y_i) \in X_{test}:argmax_r {\bar{f}_r}^{(k)} (x) =j}|$

ii. $A = \frac{1}{2} ( \bar{C} + \bar{C}^T )$

3\. Merge the two most similar subtypes to an internal subtype **c** =**a $$\cup$$ b** s.t.**(a,b)** = **$$argmax_{(i,j)}A_{i,j}$$** and add **c** to **T**.

4\. k=k-1; update all sample labels of **a** and **b** to **c** in **Y**. 
**Output:** 
The learned subtype hierarchy **T** and its associated multi-step hierarchical classifiers $${\bar{f}\.}^{\.}$$ s.

-----------------------------------------------------------------------

The training framework is designed to build hierarchical classifiers by
iteratively merging subtypes based on their classification similarity.
The process begins with the initialization of an empty subtype
hierarchy, denoted as *T*. In each iteration, a set of *k* one-vs-rest
classifiers is trained using Lasso-Logistic regression, followed by
cross-validated to assess classifier performance. The pairwise
classification error between subtypes is calculated to generate a
confusion matrix that captures the frequency of misclassifications
between subtype pairs. The confusion matrix is symmetrized to account
for reciprocal errors, providing a standardized metric for measuring
subtype similarity.

After computing the confusion matrix, the framework identifies the two
most similar subtypes based on maximal similarity and merges them into a
new internal major subtype. This newly formed major subtype is added to
the learned hierarchy, and the sample labels corresponding to the merged
subtypes are updated accordingly. This process continues iteratively,
merging the most similar subtypes at each step, until only one subtype
remains. The final output consists of the learned hierarchical structure
of subtypes and the associated multi-step hierarchical classifiers,
forming a set of classifiers that operate all sub-levels. This approach
also enables learning of the most probably subtype relationships at the
same time.

Thus, the framework UGES consists of both the subtype hierarchy ***T***
and the set of multi-step hierarchical classifiers ***f***'s. In
application, UGES guides the decision-making process by navigating
through a sequence of hierarchical steps, progressively narrowing down
the possible subtypes. Starting from the root of the hierarchy, a sample
is passed through a series of classifiers in ***f***'s, each of which is
trained to distinguish between fewer but more similar subtypes at lower
levels of the hierarchy. This multi-step process not only improves
classification accuracy by focusing on finer distinctions between
closely related subtypes, but also reduces the complexity of the
classification task, making it more interpretable and computationally
efficient.

## D. Relative importance 
Genetic and epigenetic determinants of intrinsic subtypes were
identified by the Lasso-Logistic regression, enabling the quantification
of effects for all three types of DNA alterations. Lasso-Logistic
regression allows for high-dimensional variable compression and
selection. Consequently, the corresponding regression coefficients of
the selected DNA alterations by Lasso were summarized to assess the
importance. In total, 236, 367, and 530 DNA alterations were selected
from the B-vs-(H,LA,LB), H-vs-(LA,LB), and LA-vs-LB sub-classifiers of
UGES, respectively (**Supplementary Table. S2**).

We applied relative importance to measure the total effect of each DNA
alteration class \[46-48\]. A straightforward strategy for estimating
the relative importance of a DNA alteration class is to calculate the
sum of regression coefficients within that class and express the sum as
a percentage. Suppose $f = \{ V,C,M\}$ is the union set of all mutation,
CNA, and methylation alterations. Then, the relative importance
${RI}_{i}$ of alteration $i$ can be calculated by:

$$\ \ {RI}_{i} = \frac{S_{i}}{S},\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)$$

where $S_{i} = \sum_{i\  \in \ f}^{}\left| c_{i} \right|$ is the sum of
all absolute regression coefficients for a certain DNA alteration class
and $S = \ \sum_{f}^{}S_{f}$.

## E. Survival analysis

Univariate and multivariate survival analyses were performed in R using
the *survival* and *survminer* packages, with Overall Survival (OS) as
the primary endpoint. For the multivariate analysis, Cox regression was
used to assess the effects of UGES subtypes and age with the *coxph()*
function. A forest plot was generated to show survival outcomes
according to UGES subtypes and age, using a resampled balanced dataset
(304 samples for each subtype) as the reference. P-values were adjusted
using the Benjamini-Hochberg (BH) method to control for False Discovery
Rate, with significance determined based on a BH-corrected threshold of
0.05.

## F. Differential analyses
Differential analyses were conducted in R on a set of 112 DNA alteration
features, including the top 10% of alterations selected by
sub-classifier. For mutation data, the Chi-Square Test was applied using
the *chisq.test()* function. For copy number aberration data, the
Kruskal-Wallis Test was performed using the *kruskal.test()* function.
For methylation data, the Welch ANOVA Test was performed using the
*oneway.test()* function, with parameter *var.equal = FALSE*. Features
were considered significantly different if the BH-corrected test p-value
was \< 0.01. All the statistical tests were two-sided.

# III. Results

## A. Cohort characteristics and data separation

We downloaded and analyzed DNA-level and clinical data for a total of
2065 samples, including 931 TCGA samples and 1134 METABRIC samples
(**Fig. 1b**). We assigned these samples to three breast cancer cohorts:
TCGA-only cohort, METABRIC-only cohort, and combined TCGA+METABRIC
cohort. The descriptive statistics of demographics (age) and relevant
clinical variables (ER status, PR status, Her2 status) were collected
for all cohorts (see **Table. 1**).

Univariate survival analyses of clinical variables also yielded
confirmatory statistics suggesting that all the cohorts represented the
general breast cancer population and were suitable for integrative
analysis. Elderly patients had a statistically significant worse
survival rate than expected \[49, 50\], with a relative increase of
death risk of approximately 4% (HR = 1.039, 95% CI = \[1.029, 1.049\],
P \<2e-16), demonstrating that the clinical variables were generally
consistent across the cohorts. Patients with positive ER status or
negative Her2 status also had a statistically significantly improved
survival rate than expected \[51, 52\], with a relative reduction of
death risk by \>25% or \>64% (positive ER status: HR = 0.749, 95%
CI = \[0.563, 0.997\], P=0.0474; positive Her2 status: HR = 1.643, 95%
CI = \[1.244, 2.169\], P=0.0005).

We divided the TCGA+METABRIC cohorts into training and test sets (**Fig.
1c**). To ensure balanced representation of the four intrinsic subtypes
in the training set, we randomly selected 200 samples from each subtype,
resulting in a total of 800 samples for training, with the remaining
1,265 samples used as the test set. We trained and cross-validated
DNA-level classifiers based on all strategies with the training set and
then applied the models to the test set. We used Receiver Operating
Characteristic (**ROC**) curves and Area Under Curve (**AUC**) scores as
performance evaluation metrics. For technical details of data
processing, see **Materials and Methods**.

## B. UGES \-- a sufficient and accurate DNA-level classifier

We first developed a sufficient yet trimmed DNA-level classifier for
distinguishing intrinsic subtypes. One drawback of multi-omics subtyping
is its mixed and undifferentiated molecular features, including genomic,
epigenomic, and transcriptomic data. In contrast, a restrained yet
effective DNA-level classifier can improve model interpretability,
applicability, and robustness. To achieve this, it is essential to
determine the optimal classification strategy and the specific DNA-level
features to include. Both the TCGA and METABRIC datasets provide
comprehensive DNA-level alterations, including mutations, CNAs, and
methylations. Before using the complete set of DNA-level alterations to
evaluate candidate hierarchical classification, we assessed the
informativeness of individual DNA-level alteration categories in
combination with four potential strategies.

The first candidate strategy was a three-step hierarchical
classification approach derived from the results of the proposed
hierarchical learning algorithm (see **Algorithm 1** in **Materials and
Methods**). The approach was encoded as **(B,(H,(LA,LB)))** by applying
the standard *Newick* tree annotation (**Fig. 2a**) \[53\], where **B,
H, LA,** and **LB** denote the **B**asal, **H**er2, **L**um-**A,** and
**-B** subtypes, respectively. In this strategy, the Basal subtype was
separated from all other subtypes at the first step. Then, Her2 was
separated at the second step. Finally, Luminal-A and Luminal-B were
distinguished in the last step.

We also applied three other classification strategies (**Fig. 2b**). The
first is the most common simple one-step multi-class classification
\[54\], denoted as **(B,H,LA,LB)**, which classifies four intrinsic
subtypes all at once. The second strategy, **(B,H,(LA,LB))**, is a
two-step classification inspired by the clinical decision process. In
this approach, any given patient was first tested for the Basal and Her2
biomarkers, followed by LumA and LumB sub-markers for further subtyping
\[54, 55\]. This strategy is reasonable, as the gene expression profile
of Luminal-like subtypes are considered more like each other compared to
other intrinsic subtypes \[56, 57\]. Finally, to address the
underrepresentation of Basal and Her2 subtypes in the cohorts, we
introduced a new strategy **((B,H),(LA,LB))**. In this approach, Basal
and Her2 subtypes are grouped at the first classification step to boost
their sample size, thereby balancing the training set. This adjustment
improves learning efficiency and model generalizability \[58\].

Overall, we found that the self-learned (B,(H,(LA,LB))) model using our
hierarchical learning algorithm, outperformed the other strategies.
Indeed, all multi-step strategies had higher performance compared to the
one-step approach, highlighting the importance of incorporating
hierarchical classification into the task. As shown in **Fig. 2c**, the
training and testing overall AUC scores (averaged over all subtypes) for
(B,(H,(LA,LB))) were 0.984 and 0.963, respectively. In comparison, the
scores for ((B,H),(LA,LB)) were 0.991 and 0.952, for (B,H,(LA,LB)) were
0.984 and 0.956, and for (B,H,LA,LB), they were only 0.965 and 0.889.
The fact that all AUC scores for these strategies were greater than 0.9
supports our conclusion that DNA-level information is sufficient and
robust for identifying intrinsic subtypes accurately.

Notably, all multi-step strategies grouped Luminal-like subtypes at the
initial step, resulting in a substantial improvement in the overall AUC
by at least 6.3% compared to the one-step strategy model, which was
underlying PAM50 and many other multi-omics subtyping methods. This
result gave strong empirical evidence of the common genetic and
epigenetic basis of Luminal-like tumors as they diverged further from
Basal and Her2 subtypes. Such commonality was completely not exploited
in one-step multi-class subtyping. We also compared the three strategies
using the ROC curves of the combined set by the one-vs-one method
(**Supplementary Fig. S1**). It was obvious that (B,H,LA,LB) was the
worst at distinguishing LumA and LumB (AUC=0.828) compared to the
multi-step classifiers AUC=0.948 in (B,(H,(LA,LB))), AUC=0.941 in
((B,H),(LA,LB)), AUC=0.951 in (B,H,(LA,LB)). The latter improves the
separation accuracy by 11.3% in distinguishing Luminal A and B subtypes,
which is the most difficult case.

We found that the (B,(H,(LA,LB))) strategy is the best for identifying
Her2 and Luminal-like subtypes. In terms of the test set, its AUC scores
are 0.987, 0.974, 0.96, and 0.933 for Basal, Her2, LumA, and LumB
subtypes, respectively, while those for ((B,H),(LA,LB)) are 0.991,
0.956, 0.953 and 0.909; for (B,H,(LA,LB)) are 0.992, 0.94, 0.96 and
0.933; and for (B,H,LA,LB) are 0.988, 0.935, 0.852 and 0.782,
respectively (**Fig. 2c**). In fact, (B,(H,(LA,LB))) gained performance
considerably (AUC +0.018) in identifying the Her2 subtype, and performed
best in identifying LumA and LumB subtypes, with only a slight loss in
identifying the Basal cases (AUC -0.005). Considering the high
performance of (B,(H,(LA,LB))), and it was self-learned by the
hierarchical learning algorithm, we selected and implemented it into the
final UGES classifier.

We also found that multi-omics classifiers using all DNA-level
alterations significantly outperformed those relying on a single type of
alteration. For example, in the test set, the multi-omics classifier
(B,(H,(LA,LB))) achieved an overall AUC of 0.963. However, the AUCs for
mutation-only, CNA-only, and methylation-only classifiers were only 0.7,
0.877, and 0.935, respectively, all of which were notably lower.
Subtype-wise, the multi-omics classifier (B,(H,(LA,LB))) also
consistently outperformed the corresponding single-feature classifiers
(**Fig. 2d**). The same pattern was observed for classifiers based on
other strategies (**Supplementary Fig. S2**), except for the
methylation-only (B,H,LA,LB) classifier, which performed better
(AUC=0.943) than the multi-omics classifier (AUC= 0.935) in identifying
Her2 samples. This exception may be attributed to model overfitting
specific to the methylation training data, as discussed later (see
**Discussion**). Based on these findings, we decided to include all
mutations, CNAs, and methylation features in the final classifier. We
finally built UGES as a three-step all-feature (B,(H,(LA,LB)))
classifier, which achieves the highest overall AUC (0.963) in our
benchmark analysis.

Additionally, we trained UGES on TCGA and tested it on METABRIC and vice
versa. When UGES model was trained solely on METABRIC and tested on
TCGA, the overall AUC on the test set was 0.735 (**Supplementary Fig.
S3**). In the reverse scenario, when UGES was trained on TCGA and tested
on METABRIC, the overall AUC improved to 0.846 (**Supplementary Fig.
S3**). These results are acceptable, as training a model on a single
cohort typically results in lower performance when tested on an
independent cohort. This is primarily due to the batch effect, including
differences in patient demographics, sample processing, and sequencing
protocols between cohorts. Such differences can introduce domain shifts
that limit a model's ability to generalize across datasets while trained
solely on a single source. The observed performance decline in both
scenarios reflects the inherent limitations of cross-cohort validation,
emphasizing the importance of integrative multi-source training.
Exposing the model to a broader range of variability and clinical
characteristics allows it to address better the challenges posed by
dataset heterogeneity.

We further compared UGES with leading methylation and/or transcription
profile based methods, such as XAI-MethylMarker \[38\] and moBRCA-net
\[37\] (**Supplementary Fig. S4**). With the entire test set, which
included both TCGA and METABRIC samples, UGES model outperformed
XAI-MethylMarker across all metrics, including AUC, accuracy, F1-score,
precision, and recall (**Supplementary Fig. S4a**). Specifically, UGES
achieved an AUC of 0.963, compared to XAI-MethylMarker's 0.898. Further,
to make a fair comparison with moBRCA-net \[37\], we used only the TCGA
part of the test set. UGES outperformed moBRCA-net in all three metrics:
accuracy (0.9154 vs. 0.891), weighted F1-score (0.921 vs. 0.887), and
Matthews Correlation Coefficient (**MCC**) (0.8487 vs. 0.831)
(**Supplementary Fig. S4b**). These results underscore the effectiveness
of the UGES model's DNA-based design, supporting its potential as a
robust and effective solution for clinical applications such as liquid
biopsy.

## C. methylation has the main effect in defining subtypes 

We studied the relative importance of each DNA feature based on
UGES-defined intrinsic subtypes. Biologically, genetic and epigenetic
alterations can have different impacts on intrinsic subtypes via varied
biological mechanisms \[59, 60\], including gene functional disruption,
gene dosage changes, and genome accessibility changes. The availability
of full DNA-level data on a large scale allowed us to quantify these
effects based on the UGES subtyping and to explore subtype-wise early
detection and precise treatment options. Relative importance is a metric
that calculates the relative contribution of predictive variables (see
**Materials and Methods**). The relative importance of DNA features in
the UGES sub-classifiers *B-vs-(H,LA,LB)*, *H-vs-(LA,LB)*, and
*LA-vs-LB* are shown in **Table. 2**.

We found that DNA methylation changes were consistently the most
important alteration in defining all sub-classifiers. The relative
importance of methylation features was 53.97%, 55.08% and 51.62% for
*B-vs-(H,LA,LB)*, *H-vs-(LA,LB)*, and *LA-vs-LB* sub-classifiers,
respectively, all over 50%. In general, CNA is a remote second in
importance, accounting for 33.35%, 25.15%, and 19.76% of the total
importance in the respective models, and mutation is in third place,
accounting for 12.68%, 19.78%, and 28.63% of the total in the respective
models. Consistent with our quantitative assessment here, the
methylation-feature-only classifiers had the highest AUC scores for all
subtypes (see **Fig. 2d)**, suggesting that methylation change as a
whole is probably the most prominent feature in defining the intrinsic
status of cancerous cells in breast cancer, which were also reported in
previous studies alphabet in smaller cohort sizes \[61, 62\].

Either CNA or mutation features could be the second most important,
depending on the sub-classifier. CNA features had the second largest
relative importance in the *B-vs-(H,LA,LB)* and *H-vs-(LA,LB)*
sub-classifiers, 33.35% and 25.15%, respectively, but the smallest
relative importance in the *LA-vs-LB* sub-classifier, at 19.76%.
Mutation features had the second largest relative importance in the
*LA-vs-LB* sub-classifier at 28.63% but the smallest in the
*B-vs-(H,LA,LB)* and *H-vs-(LA,LB)* sub-classifiers, at 12.68% and
19.78%, respectively. Notably, the relative importance of mutation
features was the second highest (28.63%) in the *LA-vs-LB*
sub-classifier but the lowest in the other sub-classifiers, suggesting
that mutations were crucial in distinguishing LumA and LumB subtypes, as
reported by previous research \[63\].

## D. UGES-defined subtypes show improved clinical relevance

We found that UGES subtypes were substantially different from PAM50
subtypes. The alluvial plot and the confusion matrix showed that subtype
reassignment occurred for 11.57% (239 in 2065 samples) (**Figs. 3ab**).
The reassignment better resolves the classification of LumA and LumB
subtypes. Among Basal, Her2, LumA, and LumB PAM50 subtypes, 5.79%,
7.63%, 15.64%, and 9.87% changed subtypes by UGES, with the most
frequent change being 11.25% LumA (PAM50) to LumB (UGES). This is
plausible because Luminal-like types are genetically and epigenetically
more alike than others and are the most challenging to resolve.

UGES increased the model sensitivity to detect the most difficult Her2
intrinsic subtype (+28.92% after reassignment). 94.44% of the Basal
(PAM50) subtype changed to Her2 (UGES) subtype (17 out of 18 switchings)
and 67.31% of LumB (PAM50) subtype changed to Her2 (UGES) subtype (35
out of 52 switchings), while the frequency of change to other subtypes
was much lower, with 0% and 5.56% of Basal (PAM50) subtype changed to
LumA and LumB (UGES) subtypes, and 7.69% and 25% of LumB (PAM50) subtype
changed to Basal and LumA (UGES) subtypes, respectively.

We then demonstrated that UGES subtypes were more predictive of
patients' overall survival as evidence of their validity and clinical
relevance. We found that UGES subtypes showed improved prognostic value
with better overall survival (OS) stratification compared to PAM50
subtypes (**Fig. 3cd**). The results were based on both univariate and
multivariate survival analyses (see **Materials and Methods**). Based on
univariate survival analysis, in five out of six pairwise comparisons
(each comparing two subtypes), we found that the survival of UGES
subtypes was more distinguishable from each other than PAM50 subtypes
(**Fig. 3c**). For example, the statistical significance of survival
differences among these UGES subtypes increased (smaller p-values),
except between Her2 and LumB subtypes, proving that the UGES subtypes
are more clinically relevant. Moreover, under UGES, the survival
differences between the LumA and Basal subtypes became statistically
significant (p-value decreased from 0.07 to 0.04). In general, UGES
subtyping resulted in sharper differences between clinical intrinsic
subtype groups.

In multivariate survival analysis, where age was added as a risk factor,
UGES subtyping also showed improved clinical relevance. Using a
resampled balanced sample set as reference, UGES subtypes were also
found to better stratify the intrinsic subtypes for their survival
difference (P = 9.7e-55), as compared to PAM50 subtypes (P = 2.2e-47)
(**Fig. 3d**). For example, the significance of the differences between
Her2 and the other subtypes had increased, changing from statistically
insignificant to significant (p-value decreased from 0.0557 to 0.0048).
Moreover, the LumA subtype had a higher overall survival rate (Hazard
Ratio \< 1) compared to the other subtypes, verifying the better
prognosis of LumA patients, as reported by others \[64\]. These findings
proved that UGES subtypes yielded more clinically distinguishable
patient subgroups, influential for more precise clinical management.

## E. UGES identified signature subtype-specific alterations

Unraveling genetic and epigenetic heterogeneity and identifying
subtype-specific markers of breast cancer are key to developing targeted
therapies that will lead to improved patient outcomes \[65\]. Here, we
conducted differential analyses using 112 genetic and epigenetic
features selected as the top 10% most important features in each UGES
sub-classifier (see **Materials and Methods**). As a result, we
identified 52 subtype-delineating signature alterations. These signature
DNA features are presented in **Table. 3** grouped by their
corresponding differential subtypes. Heatmap plots of signature DNA
alteration profiles in **Fig. 4** also showed patterns consistent with
the results given in **Table. 3**.

The heat patterns are, in general, consistent with the relative
importance of DNA alteration classes previously discussed (**Fig. 4**).
For example, the B and (H,LA,LB) subtypes differed little in mutation
profiles but significantly in CNA and methylation profiles. According to
the relative importance analysis, CNA and methylation alterations
totally accounted for 87.02% of the feature importance for the
*B-vs-(H,LA,LB)* sub-classifier.

We observed four distinct signature heat patterns of respective DNA
alteration classes, which are explained as follows: First, gene mutation
patterns showed substantial differences between the subtype groups, with
a total of nine highly mutated signature genes, including TP53, PIK3CA,
GPR124, HRAS, SETD2, ERBB2, PNN, ATXN1, and CP (**Fig. 4a**). A few of
these alterations were well-known cancer-driving events. For example,
the TP53 gene, a tumor suppressor gene enabling tumors to be more
aggressive, is typically mutated at the early stage of cancer to confer
genetic instability to cancer cells and is associated with poor
prognosis \[66\]. We observed that the TP53 gene was highly mutated in
Basal and Her2 subtypes (**pattern ① in Fig. 4a**) and less mutated in
Luminal-like subtypes. This suggested that patients with Basal and Her2
subtypes may have a worse prognosis, as proven by many studies \[1, 6,
7\]. PIK3CA mutation was identified as (H,LA,LB) marker, which is most
common in Luminal-like breast cancer, followed by Her2 breast cancer and
triple-negative breast cancer (**pattern ② in Fig. 4a**), consistent
with previous studies \[67, 68\]. We also found HRAS, SETD2, and ERBB2
mutations to be significantly associated with Her2 breast cancer,
consistent with previous reports \[69-71\].

UGES also identified highly mutated signature genes as novel biomarkers
of intrinsic subtypes, such as GPR124, PNN, ATXN1, and CP. Notably,
mutations in the ATXN1 gene were newly identified by UGES to be
associated with the LumB subtype. The major phenotypical difference
between LumA and LumB subtypes is that LumB tumors express a high level
of Her2 and the proliferation marker Ki-67. The Notch signaling pathway
regulates Ki-67 \[72, 73\]. Studies showed that ATXN1 is a negative
regulator of the Notch pathway, and ATXN1 knockdown enhanced tumor
invasion by activating the Notch pathway. Therefore, a mutated ATXN1
gene may lose its inhibitor role of the Notch pathway, and the abnormal
expression of the Notch pathway leads to the overexpression of Ki-67 in
Her2 cancer cells \[74, 75\]. GPR124 and CP gene mutations were also
associated with tumor progression and invasion in breast cancer
previously \[76,77\]. However, this is different with specific intrinsic
subtypes. While our results suggest that these gene mutations are
potential biomarkers of their specific intrinsic subtypes, more evidence
is required to prove their functional roles.

Second, CNA alterations were found effective in identifying Basal and
Her2 subtypes. As shown in **Fig. 4b** and **Table. 3**, Basal, Her2,
LumA, and LumB subtypes have distinct CNA patterns. For example, PGAP3
and MIR4728 gene copy numbers were highly amplified in the Her2 subtype
(**pattern ③ in Fig. 4b**), as compared to the Basal subtype. This is
consistent with previous studies. For example, PGAP3 might impact on the
Her2 subtype by being co-amplified with Her2 \[78\]; MIR4728 gene is
located within an intron of the ERBB2 gene, a well-known oncogene
strongly associated with the Her2 subtype \[79\]. Notably, UGES newly
identified MIR5684 copy number variations deletion as a novel biomarker
of the Her2 subtype, while more experimental evidence is required to
confirm the relation.

Finally, in addition to genetic differences, epigenetic alterations,
such as hyper- or hypo-methylation changes, also exhibited
subtype-specific patterns (**Fig. 4c**). For example, as shown in
**pattern ④ of Fig. 4c**, nine genes (ADCY4, CRYAB, DNM3, HOXA11, INS,
KCNH8, LAMB3, PRTN3, UCN) showed hypomethylation in the Basal subtype.
Hypomethylation generally elevates a gene's expression, and many of
those genes were shown to be highly expressed in triple-negative breast
cancer cells (TNBC), accounting for \>80% of the Basal subtype. For
example, CRYAB was known to enhance cell migration in TNBC cells,
specifically \[80\]; and LAMB3 was known to mediate the apoptotic,
proliferative, invasive, and metastatic abilities of Basal breast cancer
cells \[81\]. Notably, UGES newly found AIRE and PRSS41 are
hypermethylated in the Basal subtype, which are potential biomarkers
requiring further experimental evidence.

# IV. Discussion

In this study, we constructed a DNA-level classifier UGES using only DNA
alterations to identify breast cancer intrinsic subtypes. Through the
development, we have gained several insights. First, a multi-step
strategy with the grouping of LumA and LumB subtypes is key to the
better performance of all classifiers. As the test set result showed
(**Fig. 2c**), the best AUC scores for single and overall subtypes were
all obtained from such multi-step strategies. (B,(H,(LA,LB))) had the
best overall AUC scores, while (B,H,(LA,LB)) obtained the highest AUC
score, 0.992, for identifying the Basal case. The AUC scores of the
one-step (B,H,LA,LB) strategy were, however, consistently the lowest, in
particular having trouble resolving the LumA and LumB subtypes, with
AUC=0.852 and 0.782, respectively. One potential explanation is that the
one-step classifier may need more accuracy to ignore the natural
hierarchy of intrinsic subtypes. In fact, in the past, breast cancers
have been broadly classified based on their gene expression profiles
into Luminal- and Basal-type tumors. It is only recently that the
Luminal type was further divided into two sub-groups \[57\]. This
demonstrates the fact that Luminal-like subtypes differ more from Basal
and Her2 subtypes than each other. Thus, we conclude that the
hierarchical stepwise classification method is a more realistic approach
to identifying intrinsic subtypes. The inaccuracy of the PAM50
classifier is perhaps also attributable to the fact that it took a
four-classification, one-step approach.

Second, all types of DNA alterations are informative for defining
subtypes. As shown in **Fig. 2d** and **Supplementary Fig. S2**, the
overall AUC scores of single-feature classifiers were all over 0.6 in
testing for the four strategies. Specifically, the average AUC scores of
0.698, 0.854, and 0.919 for mutation, CNA, and methylation features,
respectively. It suggested that the features are all informative.
Notably, the AUC scores of methylation-only classifiers were comparable
to the all-feature classifiers, with the lowest score being 0.879 in the
one-step multi-classification strategy. In fact, methylation changes are
becoming widely accepted as biomarkers for early breast cancer detection
\[82\]. Biologically, promoter hypermethylation is a more frequent event
than mutations in the process of carcinogenesis \[83\], with estimates
varying from 600 to 1000 aberrantly methylated genes per tumor \[84\].
This significant contribution of variability may explain the strong
effect of methylation in defining intrinsic subtypes and its high
relative importance, as observed.

Even though methylation has a dominant influence on breast cancer
subtyping, mutations and copy number alterations (CNAs) still contribute
significantly to certain subtypes. For instance, somatic mutations in
genes such as **TP53** and **PIK3CA** are crucial in the Basal-like and
Her2-positive subtypes, serving as key driver alterations. Similarly,
CNAs play a pivotal role in these subtypes, with **PGAP3** copy number
gain showing high relevance in Her2-positive cases, indicating that
specific chromosomal changes are important in defining these subtypes.
While methylation, particularly in promoter regions, is a key marker for
early detection and subtype classification, especially in Luminal A and
Luminal B subtypes, more than relying on methylation data may be
required \[85-87\]. Our findings suggest that integrating mutations,
CNAs, and methylation provides a more comprehensive approach, improving
the model\'s accuracy and generalization ability while reducing the risk
of overfitting that can occur when using a single feature type.

For example, the methylation-only (B,H,LA,LB) classifier performed
better in identifying the Her2 subtype (AUC = 0.943) compared to the
all-feature classifier (AUC = 0.935). However, if we looked into their
specific training and testing performance, the overall AUC scores for
the training and test sets of the methylation-only (B,H,LA,LB)
classifier were 1 and 0.879, while these for the all-feature (B,H,LA,LB)
classifier, were 0.965 and 0.889, as shown in **Supplementary Fig.
S2c**. The larger difference in training and test set performance of the
methylation-only classifier demonstrated that it may be overfitted.
Thus, only using methylation features to analyze breast cancer intrinsic
subtypes would result in the loss of model generalizability.

Third, our findings provide a foundation for the development of more
specific and mutually exclusive DNA biomarkers for distinguishing breast
cancer subtypes. For instance, in our dataset, among the 716 samples
with **TP53** mutations, only 199 also had **PIK3CA** mutations;
similarly, of the 824 samples with **PIK3CA** mutations, only 199 had
**TP53** mutations. This reveals a significant mutual exclusivity
between **PIK3CA** and **TP53** mutations in breast cancer, a pattern
consistent with previous research \[88\]. Such mutually exclusive
mutation patterns highlight the potential to use these alterations as
specific biomarkers to refine subtype classification and enhance
targeted therapeutic strategies.

Fourth, the 52 DNA biomarkers we identified hold significant clinical
potential. Notably, the ATXN1 mutation shows potential specificity for
the LumB subtype. Since LumA and LumB subtypes differ in treatment
strategies and prognosis, ATXN1 could serve as a valuable marker in
early screening. Additionally, the association between the MIR5684 copy
number loss and the Her2 subtype suggests that MIR5684 could enhance
early screening for Her2-positive patients. Furthermore,
hypermethylation of the AIRE and PRSS41 genes may serve as biomarkers
for the Basal subtype detection.

Beyond the signature alterations in well-known genes, many of our newly
identified biomarkers have the potential to become emerging targets in
breast cancer precision medicine. Some lesser-known but promising
targets include GPR124, which influences tumor growth and metastasis in
triple-negative breast cancer (TNBC) by regulating the Wnt signaling
pathway. Developing small molecule inhibitors or antibodies against
GPR124 could be a potential avenue for future research. Additionally,
since ATXN1 is a negative regulator of the Notch signaling pathway,
targeting ATXN1 mutations or inhibiting the Notch pathway could provide
new treatment strategies for LumB patients, particularly given their
poor response to chemotherapy. The discovery of these biomarkers lays a
good foundation for further research.

Another potential application for UGES is in facilitating the
development of early subtyping using ctDNA liquid biopsy. Although
promising, ctDNA-based breast cancer screening faces several technical
challenges. First, cancer heterogeneity presents a challenge because
molecular subtyping of breast cancer involves a complex interplay of
genes and mutations. Since ctDNA data can only capture a portion of this
information, it may not fully represent the heterogeneity of the tumor,
limiting the accuracy of models trained on such data. Second, low tumor
burden in the early stages of cancer results in low ctDNA concentration
and mutation frequencies, making it difficult to distinguish genuine
cancer-related mutations from background noise. Third, there needs to be
more sufficient ctDNA datasets, and the generalizability of trained
models could be better. High-quality labeled data, such as whole-genome
methylation, mutation, or copy number ctDNA data, is difficult to obtain
due to the scarcity of early-stage cancer samples. These limitations
affect the accuracy of detection.

The insights we gained from building the UGES model can help address
these challenges. First, the data used to train UGES comes from solid
tumor samples, which provide higher DNA yield and quality, allowing for
the more comprehensive capture of intra-tumor epi/genetic heterogeneity.
Multi-patient populations, such as those represented by TCGA and
METABRIC, enhance the capture of inter-tumor epi/genetic heterogeneity.
Second, the high tumor burden and quality of these tumor samples ensure
lower background noise, which minimizes false positives and false
negatives during classifier training. Third, we carefully integrated
data into a large omics dataset with 50,831 features, including
mutations, copy number aberrations, and methylation patterns, from 2,065
samples. This greatly enriches the training data and improves our
model\'s ability to manage tumor heterogeneity across subtypes. Lastly,
we employed a novel multi-step hierarchical learning algorithm, rather
than a traditional single-step multiclass classification approach. This
strategy enables the model to learn the subtle differences between
subtypes at each step, significantly enhancing its generalizability
while identifying subtype-specific biomarkers. These biomarkers could
potentially be used for early breast cancer screening.

Finally, our study does have several limitations worth noting. We used
PAM50 subtypes as the initial training labels to construct the UGES
classifier despite the fact that accurately determining the true
intrinsic subtypes of breast cancer experimentally is challenging \[4\].
We opted for the PAM50 label because it is one of the most
internationally recognized for classifying intrinsic subtypes. However,
unlike the typical approach that uses expression data, UGES was
developed using only DNA alterations. This method yielded highly
accurate results and produced clinically more distinguishable subtype
groups, demonstrating that DNA alterations can serve as sufficient and
reliable features. Nevertheless, UGES may not capture some emerging DNA
mutations or features. While our current sample size is adequate to
distinguish the major subtypes, smaller subtypes have been difficult to
identify due to insufficient data. As new subtype discoveries emerge
\[89-91\] and more data become available, future models should
incorporate these refined subtypes and additional feature information.
Also, we have not yet considered chromosomal-level structural changes
\[92, 93\] and mid-size insertion-deletions, which challenge detection
by next-generation sequencing-based technologies \[94\] but are
important to understanding breast cancers. Potential bias from
confounding clinical factors, such as tumor size and lymph node
metastasis, could influence survival outcomes and impact the
interpretation of subtypes and survival. Future analyses should consider
additional clinical variables, such as tumor size and treatment
modalities, to further adjust for these confounders.

# Data Availability Statement

The TCGA and METABRIC datasets are available in GDC
(<https://gdc.cancer.gov/about-data/publications/pancanatlas>) and
cBioPortal (<http://www.cbioportal.org/>). Analytical codes are provided
at <https://github.com/labxscut/UGES>.

\[1\] X. F. Dai, T. Li, Z. H. Bai, Y. K. Yang, X. X. Liu, J. L. Zhan,
and B. Z. Shi, "Breast cancer intrinsic subtype classification, clinical
use and future trends," *American Journal of Cancer Research,* vol. 5,
no. 10, pp. 2929-2943, 2015.

\[2\] C. Horr, and S. A. Buechler, "Breast Cancer Consensus Subtypes: A
system for subtyping breast cancer tumors based on gene expression,"
*Npj Breast Cancer,* vol. 7, no. 1, Oct 12, 2021.

\[3\] D. M. Wolf, C. Yau, J. Wulfkuhle, L. Brown-Swigart, R. I.
Gallagher, P. R. E. Lee, Z. Zhu, M. J. Magbanua, R. Sayaman, N.
O\'Grady, A. Basu, A. Delson, J. P. Coppe, R. X. Lu, J. Braun, S. M.
Asare, L. Sit, J. B. Matthews, J. Perlmutter, N. Hylton, M. C. Liu, P.
Pohlmann, W. F. Symmans, H. S. Rugo, C. Isaacs, A. M. DeMichele, D. Yee,
D. A. Berry, L. Pusztai, E. F. Petricoin, G. L. Hirst, L. J. Esserman,
L. J. V. Veer, and I.-S. Investigators, "Redefining breast cancer
subtypes to guide treatment prioritization and maximize response:
Predictive biomarkers across ten cancer therapies," *Cancer Cell,* vol.
40, no. 6, pp. 609-+, Jun 13, 2022.

\[4\] F. Schettini, F. Braso-Maristany, N. M. Kuderer, and A. Prat, "A
perspective on the development and lack of interchangeability of the
breast cancer intrinsic subtypes," *Npj Breast Cancer,* vol. 8, no. 1,
Jul 19, 2022.

\[5\] J. S. Parker, M. Mullins, M. C. Cheang, S. Leung, D. Voduc, T.
Vickery, S. Davies, C. Fauron, X. He, Z. Hu, J. F. Quackenbush, I. J.
Stijleman, J. Palazzo, J. S. Marron, A. B. Nobel, E. Mardis, T. O.
Nielsen, M. J. Ellis, C. M. Perou, and P. S. Bernard, "Supervised risk
predictor of breast cancer based on intrinsic subtypes," *Journal of
Clinical Oncology,* vol. 27, no. 8, pp. 1160-7, Mar 10, 2009.

\[6\] A. Goldhirsch, W. C. Wood, A. S. Coates, R. D. Gelber, B.
Thurlimann, H. J. Senn, and P. Members, "Strategies for subtypes-dealing
with the diversity of breast cancer: highlights of the St Gallen
International Expert Consensus on the Primary Therapy of Early Breast
Cancer 2011," *Annals of Oncology,* vol. 22, no. 8, pp. 1736-1747, Aug,
2011.

\[7\] A. Goldhirsch, E. P. Winer, A. S. Coates, R. D. Gelber, M.
Piccart-Gebhart, B. Thurlimann, H. J. Senn, and m. Panel, "Personalizing
the treatment of women with early breast cancer: highlights of the St
Gallen International Expert Consensus on the Primary Therapy of Early
Breast Cancer 2013," *Annals of Oncology,* vol. 24, no. 9, pp. 2206-23,
Sep, 2013.

\[8\] R. R. L. Bastien, A. Rodriguez-Lescure, M. T. W. Ebbert, A. Prat,
B. Munarriz, L. Rowe, P. Miller, M. Ruiz-Borrego, D. Anderson, B. Lyons,
I. Alvarez, T. Dowell, D. Wall, M. A. Segui, L. Barley, K. M. Boucher,
E. Alba, L. Pappas, C. A. Davis, I. Aranda, C. Fauron, I. J. Stijleman,
J. Palacios, A. Anton, E. Carrasco, R. Caballero, M. J. Ellis, T. O.
Nielsen, C. M. Perou, M. Astill, P. S. Bernard, and M. Martin, "PAM50
Breast Cancer Subtyping by RT-qPCR and Concordance with Standard
Clinical Molecular Markers," *BMC Medical Genomics,* vol. 5, Oct 4,
2012.

\[9\] Y. T. Chang, Z. J. Hong, J. C. Yu, W. Z. Lin, T. Y. Huang, H. H.
Tsai, A. C. Feng, K. F. Hsu, C. C. Huang, C. M. Chu, C. M. Liang, G. S.
Liao, \"Advancing breast cancer subtyping: optimizing
immunohistochemical staining classification with insights from
real-world Taiwanese data,\" *American Journal of Cancer Research*, vol.
13, no. 11, pp. 5719-5732, Nov. 2023. PMID: 38058819; PMCID:
PMC10695790.

\[10\] H. K. Kim, K. H. Park, Y. Kim, S. E. Park, H. S. Lee, S. W. Lim,
J. H. Cho, J. Y. Kim, J. E. Lee, J. S. Ahn, Y. H. Im, J. H. Yu, and Y.
H. Park, "Discordance of the PAM50 Intrinsic Subtypes Compared with
Immunohistochemistry-Based Surrogate in Breast Cancer Patients:
Potential Implication of Genomic Alterations of Discordance," *Cancer
Research and Treatment,* vol. 51, no. 2, pp. 737-747, Apr, 2019.

\[11\] J. M. Raser, and E. K. O\'Shea, "Noise in gene expression:
Origins, consequences, and control," *Science,* vol. 309, no. 5743, pp.
2010-2013, Sep 23, 2005.

\[12\] J. S. van Zon, M. J. Morelli, S. Tanase-Nicola, and P. R. ten
Wolde, "Diffusion of transcription factors can drastically enhance the
noise in gene expression," *Biophysical Journal,* vol. 91, no. 12, pp.
4350-4367, Dec, 2006.

\[13\] A. Sanchez, and I. Golding, "Genetic Determinants and Cellular
Constraints in Noisy Gene Expression," *Science,* vol. 342, no. 6163,
pp. 1188-1193, Dec 6, 2013.

\[14\] E. R. Paquet, and M. T. Hallett, "Absolute Assignment of Breast
Cancer Intrinsic Molecular Subtype," *JNCI-Journal of the National
Cancer Institute,* vol. 107, no. 1, Jan, 2015.

\[15\] M. K. Seo, S. Paik, and S. Kim, "An Improved, Assay Platform
Agnostic, Absolute Single Sample Breast Cancer Subtype Classifier,"
*Cancers,* vol. 12, no. 12, Dec, 2020.

\[16\] S. J. L. Payne, R. L. Bowen, J. L. Jones, and C. A. Wells,
"Predictive markers in breast cancer - the present," *Histopathology,*
vol. 52, no. 1, pp. 82-90, Jan, 2008.

\[17\] C. Chen, J. Peng, H. S. Xia, G. F. Yang, Q. S. Wu, L. D. Chen, L.
B. Zeng, Z. L. Zhang, D. W. Pang, and Y. Li, "Quantum dots-based
immunofluorescence technology for the quantitative determination of HER2
expression in breast cancer," *Biomaterials,* vol. 30, no. 15, pp.
2912-2918, May, 2009.

\[18\] A. R. Crowe, and W. Yue, "Semi-quantitative Determination of
Protein Expression using Immunohistochemistry Staining and Analysis: An
Integrated Protocol," *Bio Protoc,* vol. 9, no. 24, Dec 20, 2019.

\[19\] Z. Herceg, and P. Hainaut, "Genetic and epigenetic alterations as
biomarkers for cancer detection, diagnosis and prognosis," *Molecular
Oncology,* vol. 1, no. 1, pp. 26-41, Jun, 2007.

\[20\] H. Takeshima, and T. Ushijima, "Accumulation of genetic and
epigenetic alterations in normal cells and cancer risk," *Npj Precision
Oncology,* vol. 3, Mar 6, 2019.

\[21\] S. Byler, S. Goldgar, S. Heerboth, M. Leary, G. Housman, K.
Moulton, and S. Sarkar, "Genetic and Epigenetic Aspects of Breast Cancer
Progression and Therapy," *Anticancer Research,* vol. 34, no. 3, pp.
1071-1077, Mar, 2014.

\[22\] M. C. U. Cheang, S. K. Chia, D. Voduc, D. X. Gao, S. Leung, J.
Snider, M. Watson, S. Davies, P. S. Bernard, J. S. Parker, C. M. Perou,
M. J. Ellis, and T. O. Nielsen, "Ki67 Index, HER2 Status, and Prognosis
of Patients With Luminal B Breast Cancer," *Jnci-Journal of the National
Cancer Institute,* vol. 101, no. 10, pp. 736-750, May 20, 2009.

\[23\] F. Crick, "Central Dogma of Molecular Biology," *Nature,* vol.
227, no. 5258, pp. 561-563, 1970/08/01, 1970.

\[24\] J. D. Watson, and F. H. C. Crick, "Molecular Structure of Nucleic
Acids: A Structure for Deoxyribose Nucleic Acid," *Nature,* vol. 171,
no. 4356, pp. 737-738, 1953/04/01, 1953.

\[25\] I. Hoijer, A. Emmanouilidou, R. Ostlund, R. van Schendel, S.
Bozorgpana, M. Tijsterman, L. Feuk, U. Gyllensten, M. den Hoed, and A.
Ameur, "CRISPR-Cas9 induces large structural variants at on-target and
off-target sites in vivo that segregate across generations," *Nature
Communications,* vol. 13, no. 1, Feb 2, 2022.

\[26\] A. Paziewska, M. Dabrowska, K. Goryca, A. Antoniewicz, J.
Dobruch, M. Mikula, D. Jarosz, L. Zapala, A. Borowka, and J. Ostrowski,
"DNA methylation status is more reliable than gene expression at
detecting cancer in prostate biopsy," *British Journal of Cancer,* vol.
111, no. 4, pp. 781-789, Aug 12, 2014.

\[27\] X. Shao, N. Lv, J. Liao, J. B. Long, R. Xue, N. Ai, D. H. Xu, and
X. H. Fan, "Copy number variation is highly correlated with differential
gene expression: a pan-cancer study," *BMC Medical Genet,* vol. 20, no.
1, Nov 9, 2019.

\[28\] M. Rossello-Tortella, A. Bueno-Costa, L. Martinez-Verbo, L.
Villanueva, and M. Esteller, "DNA methylation-associated dysregulation
of transfer RNA expression in human cancer," *Molecular Cancer,* vol.
21, no. 1, Feb 12, 2022.

\[29\] J. Y. Kim, K. Park, H. H. Jung, E. Lee, E. Y. Cho, K. H. Lee, S.
Y. Bae, S. K. Lee, S. W. Kim, J. E. Lee, S. J. Nam, J. S. Ahn, Y. H. Im,
and Y. H. Park, "Association between Mutation and Expression of TP53 as
a Potential Prognostic Marker of Triple-Negative Breast Cancer," *Cancer
Res Treat,* vol. 48, no. 4, pp. 1338-1350, Oct, 2016.

\[30\] J. R. Pollack, T. Sorlie, C. M. Perou, C. A. Rees, S. S. Jeffrey,
P. E. Lonning, R. Tibshirani, D. Botstein, A. L. Borresen-Dale, and P.
O. Brown, "Microarray analysis reveals a major direct role of DNA copy
number alteration in the transcriptional program of human breast
tumors," *Proceedings of the National Academy of Sciences of the United
States of America,* vol. 99, no. 20, pp. 12963-12968, Oct 1, 2002.

\[31\] H. K. Solvang, O. C. Lingjaerde, A. Frigessi, A. L.
Borresen-Dale, and V. N. Kristensen, "Linear and non-linear dependencies
between copy number aberrations and mRNA expression reveal distinct
molecular pathways in breast cancer," *BMC Bioinform,* vol. 12, May 24,
2011.

\[32\] F. Jiao, S. Y. Bai, Y. Ma, Z. H. Yan, Z. Yue, Y. Yu, X. Wang, and
J. Wang, "DNA Methylation of Heparanase Promoter Influences Its
Expression and Associated with the Progression of Human Breast Cancer,"
*Plos One,* vol. 9, no. 3, Mar 14, 2014.

\[33\] B. Gyorffy, G. Bottai, T. Fleischer, G. Munkacsy, J. Budczies, L.
Paladini, A. L. Borresen-Dale, V. N. Kristensen, and L. Santarpia,
"Aberrant DNA methylation impacts gene expression and prognosis in
breast cancer subtypes," *International Journal of Cancer,* vol. 138,
no. 1, pp. 87-97, Jan 1, 2016.

\[34\] C. Curtis, S. P. Shah, S. F. Chin, G. Turashvili, O. M. Rueda, M.
J. Dunning, D. Speed, A. G. Lynch, S. Samarajiwa, Y. Y. Yuan, S. Graf,
G. Ha, G. Haffari, A. Bashashati, R. Russell, S. McKinney, A. Langerod,
A. Green, E. Provenzano, G. Wishart, S. Pinder, P. Watson, F. Markowetz,
L. Murphy, I. Ellis, A. Purushotham, A. L. Borresen-Dale, J. D. Brenton,
S. Tavare, C. Caldas, S. Aparicio, and M. Grp, "The genomic and
transcriptomic architecture of 2,000 breast tumours reveals novel
subgroups," *Nature,* vol. 486, no. 7403, pp. 346-352, Jun 21, 2012.

\[35\] M. Mohaiminul Islam, S. Huang, R. Ajwad, C. Chi, Y. Wang, and P.
Hu, "An integrative deep learning framework for classifying molecular
subtypes of breast cancer," *Comput Struct Biotechnol J,* vol. 18, pp.
2185-2199, 2020.

\[36\] Y. Q. Lin, W. Zhang, H. S. Cao, G. Y. Li, and W. Du, "Classifying
Breast Cancer Subtypes Using Deep Neural Networks Based on Multi-Omics
Data," *Genes,* vol. 11, no. 8, Aug, 2020.

\[37\] J. M. Choi and H. Chae, \"moBRCA-net: a breast cancer subtype
classification framework based on multi-omics attention neural
networks,\" BMC Bioinformatics, vol. 24, no. 1, p. 169, 2023, doi:
10.1186/s12859-023-05273-5.

\[38\] S. Rajpal, A. Rajpal, A. Saggar, A. K. Vaid, V. Kumar, M.
Agarwal, and N. Kumar, \"XAI-MethylMarker: Explainable AI approach for
biomarker discovery for breast cancer subtype classification using
methylation data,\" *Expert Syst. Appl.*, vol. 225, p. 120130, 2023.

\[39\] R. Bonneville, M. A. Krook, E. A. Kautto, J. Miya, M. R. Wing, H.
Z. Chen, J. W. Reeser, L. B. Yu, and S. Roychowdhury, "Landscape of
Microsatellite Instability Across 39 Cancer Types," *JCO Precis Oncol,*
vol. 1, 2017.

\[40\] E. Cerami, J. Gao, U. Dogrusoz, B. E. Gross, S. O. Sumer, B. A.
Aksoy, A. Jacobsen, C. J. Byrne, M. L. Heuer, E. Larsson, Y. Antipin, B.
Reva, A. P. Goldberg, C. Sander, and N. Schultz, "The cBio cancer
genomics portal: an open platform for exploring multidimensional cancer
genomics data," *Cancer Discov,* vol. 2, no. 5, pp. 401-4, May, 2012.

\[41\] T. Fleischer, A. Frigessi, K. C. Johnson, H. Edvardsen, N.
Touleimat, J. Klajic, M. L. H. Riis, V. D. Haakensen, F. Warnberg, B.
Naume, A. Helland, A. L. Borresen-Dale, J. Tost, B. C. Christensen, and
V. N. Kristensen, "Genome-wide DNA methylation profiles in progression
to in situ and invasive carcinoma of the breast with impact on gene
transcription and prognosis," *Genome Biology,* vol. 15, no. 8, 2014.

\[42\] X. F. Yang, L. Gao, and S. H. Zhang, "Comparative pan-cancer DNA
methylation analysis reveals cancer common and specific patterns,"
*Brief Bioinform,* vol. 18, no. 5, pp. 761-773, Sep, 2017.

\[43\] T. Hastie, and J. Qian, "Glmnet vignette," *Retrieved June,* vol.
9, no. 2016, pp. 1-30, 2014.

\[44\] R. Tibshirani, "Regression shrinkage and selection via the
lasso," *J R Stat Soc Series B Stat Methodol,* vol. 58, no. 1, pp.
267-288, 1996.

\[45\] S. Bengio, J. Weston, and D. Grangier, "Label embedding trees for
large multi-class tasks," *Advances in neural information processing
systems,* vol. 23, 2010.

\[46\] J. Johnson, and J. LeBreton, "History and Use of Relative
Importance Indices in Organizational Research," *Organizational Research
Methods,* vol. 7, pp. 238-257, 07/01, 2004.

\[47\] S. Tonidandel, and J. M. LeBreton, "Relative Importance Analysis:
A Useful Supplement to Regression Analysis," *Journal of Business and
Psychology,* vol. 26, no. 1, pp. 1-9, Mar, 2011.

\[48\] M. Brinkmann, L. M. Fricke, L. Diedrich, B. P. Robra, C. Krauth,
and M. Dreier, "Attributes in stated preference elicitation studies on
colorectal cancer screening and their relative importance for
decision-making among screenees: a systematic review," *Health Economics
Review,* vol. 12, no. 1, Sep 22, 2022.

\[49\] K. Meister, and J. Morgan, *Risk factors for breast cancer*: Am
Cncl on Science, Health, 2000.

\[50\] A. C. Society, "Breast cancer facts & figures 2019--2020," *Am
Cancer Soc*, pp. 1-44, 2019.

\[51\] E. E. Lower, E. L. Glass, D. A. Bradley, R. Blau, and S.
Heffelfinger, "Impact of metastatic estrogen receptor and progesterone
receptor status on survival," *Breast Cancer Research and Treatment,*
vol. 90, no. 1, pp. 65-70, Mar, 2005.

\[52\] V. R. Grann, A. B. Troxel, N. J. Zojwalla, J. S. Jacobson, D.
Hershman, and A. I. Neugut, "Hormone receptor status and survival in a
population-based cohort of patients with breast carcinoma," *Cancer,*
vol. 103, no. 11, pp. 2241-2251, Jun 1, 2005.

\[53\] J. Felsenstein, *Inferring phylogenies*: Sinauer associates
Sunderland, MA, 2004.

\[54\] A. Prat, P. Galván, B. Jimenez, W. Buckingham, H. A. Jeiranian,
C. Schaper, M. Vidal, M. Álvarez, S. Díaz, and C. Ellis, "Prediction of
Response to Neoadjuvant Chemotherapy Using Core Needle Biopsy Samples
with the Prosigna AssayProsigna ROR Score Predicts Chemosensitivity,"
*Clinical Cancer Research,* vol. 22, no. 3, pp. 560-566, 2016.

\[55\] A. Prat, P. Galván, W. Buckingham, M. Vidal, S. Díaz, P.
Nuciforo, S. Ferree, B. Adamo, S. Ramon y Cajal, and V. Peg, "Abstract
P6-01-06: Feasibility of the PROSIGNA® multigene test in core biopsies
and comparison to corresponding surgical breast cancer sections,"
*Cancer Research,* vol. 75, no. 9_Supplement, pp. P6-01-06-P6-01-06,
2015.

\[56\] C. Jones, A. Mackay, A. Grigoriadis, A. Cossu, J. S. Reis-Filho,
L. Fulford, T. Dexter, S. Davies, K. Bulmer, E. Ford, S. Parry, M.
Budroni, G. Palmieri, A. M. Neville, M. J. O\'Hare, and S. R. Lakhani,
"Expression profiling of purified normal human luminal and myoepithelial
breast cells: Identification of novel prognostic markers for breast
cancer," *Cancer Research,* vol. 64, no. 9, pp. 3037-3045, May 1, 2004.

\[57\] P. J. Keller, L. M. Arendt, A. Skibinski, T. Logvinenko, I.
Klebba, S. M. Dong, A. E. Smith, A. Prat, C. M. Perou, H. Gilmore, S.
Schnitt, S. P. Naber, J. A. Garlick, and C. Kuperwasser, "Defining the
cellular precursors to human breast cancer," *Proceedings of the
National Academy of Sciences of the United States of America,* vol. 109,
no. 8, pp. 2772-2777, Feb 21, 2012.

\[58\] L. L. Cheng, J. P. Zhang, J. Yang, and J. Ma, "An Improved
Hierarchical Multi-Class Support Vector Machine with Binary Tree
Architecture," *Icicse: 2008 International Conference on Internet
Computing in Science and Engineering, Proceedings*, pp. 106-109, 2008.

\[59\] D. H. Roukos, "Effect of genetic cancer risk assessment on
surgical decisions at breast cancer diagnosis - Invited critique,"
*Archives of Surgery,* vol. 138, no. 12, pp. 1329-1329, Dec, 2003.

\[60\] M. Kumaran, C. E. Cass, K. Graham, J. R. Mackey, R. Hubaux, W.
Lam, Y. Yasui, and S. Damaraju, "Germline copy number variations are
associated with breast cancer risk and prognosis," *Scientific Reports,*
vol. 7, Nov 7, 2017.

\[61\] Y. Gao, M. Widschwendter, and A. E. Teschendorff, "DNA
Methylation Patterns in Normal Tissue Correlate more Strongly with
Breast Cancer Status than Copy-Number Variants," *EBioMedicine,* vol.
31, pp. 243-252, 2018/05/01/, 2018.

\[62\] S. Yamashita, T. Kishino, T. Takahashi, T. Shimazu, H. Charvat,
Y. Kakugawa, T. Nakajima, Y.-C. Lee, N. Iida, and M. Maeda, "Genetic and
epigenetic alterations in normal tissues have differential impacts on
cancer risk among tissues," *Proceedings of the national academy of
sciences,* vol. 115, no. 6, pp. 1328-1333, 2018.

\[63\] X. Li, J. Zhou, M. Xiao, L. Zhao, Y. Zhao, S. Wang, S. Gao, Y.
Zhuang, Y. Niu, S. Li, X. Li, Y. Zhu, M. Zhang, and J. Tang, "Uncovering
the Subtype-Specific Molecular Characteristics of Breast Cancer by
Multiomics Analysis of Prognosis-Associated Genes, Driver Genes,
Signaling Pathways, and Immune Activity," *Front Cell Dev Biol,* vol. 9,
pp. 689028, 2021.

\[64\] O. Yersal, and S. Barutca, "Biological subtypes of breast cancer:
Prognostic and therapeutic implications," *World J Clin Oncol,* vol. 5,
no. 3, pp. 412-24, Aug 10, 2014.

\[65\] K. S. Johnson, E. F. Conant, and M. S. Soo, "Molecular Subtypes
of Breast Cancer: A Review for Breast Radiologists," *Journal of Breast
Imaging,* vol. 3, no. 1, pp. 12-24, 2020.

\[66\] L. A. Donehower, T. Soussi, A. Korkut, Y. X. Liu, A. Schultz, M.
Cardenas, X. B. Li, O. Babur, T. K. Hsu, O. Lichtarge, J. N. Weinstein,
R. Akbani, and D. A. Wheeler, "Integrated Analysis of TP53 Gene and
Pathway Alterations in The Cancer Genome Atlas (vol 28, pg 1370, 2019),"
*Cell Rep,* vol. 28, no. 11, pp. 3010-3010, Sep 10, 2019.

\[67\] R. Arsenic, A. Lehmann, J. Budczies, I. Koch, J. Prinzler, A.
Kleine-Tebbe, C. Schewe, S. Loibl, M. Dietel, and C. Denkert, "Analysis
of PIK3CA mutations in breast cancer subtypes," *Appl Immunohistochem
Mol Morphol,* vol. 22, no. 1, pp. 50-6, Jan, 2014.

\[68\] E. J. Anderson, L. E. Mollon, J. L. Dean, T. L. Warholak, A.
Aizer, E. A. Platt, D. H. Tang, and L. E. Davis, "A Systematic Review of
the Prevalence and Diagnostic Workup of PIK3CA Mutations in HR+/HER2-
Metastatic Breast Cancer," *Int J Breast Cancer,* vol. 2020, pp.
3759179, 2020.

\[69\] M. de Oliveira Taveira, S. Nabavi, Y. Wang, P. Tonellato, F. J.
Esteva, L. C. Cantley, and G. M. Wulf, "Genomic characteristics of
trastuzumab-resistant Her2-positive metastatic breast cancer," *Journal
of Cancer Research and Clinical Oncology,* vol. 143, no. 7, pp.
1255-1262, 2017.

\[70\] T. Yamaguchi, H. Mukai, S. Yamashita, S. Fujii, and T. Ushijima,
"Comprehensive DNA methylation and extensive mutation analyses of
HER2-positive breast cancer," *Oncology,* vol. 88, no. 6, pp. 377-384,
2015.

\[71\] A. Majumder, M. Sandhu, D. Banerji, V. Steri, A. Olshen, and M.
M. Moasser, "The role of HER2 and HER3 in HER2-amplified cancers beyond
breast cancers," *Scientific Reports,* vol. 11, no. 1, pp. 9091,
2021/04/27, 2021.

\[72\] K. Holli-Helenius, A. Salminen, I. Rinta-Kiikka, I. Koskivuo, N.
Brück, P. Boström, and R. Parkkola, "MRI texture analysis in
differentiating luminal A and luminal B breast cancer molecular
subtypes-a feasibility study," *BMC Medical Imaging,* vol. 17, no. 1,
pp. 1-9, 2017.

\[73\] M. Yanagawa, K. Ikemot, S. Kawauchi, T. Furuya, S. Yamamoto, M.
Oka, A. Oga, Y. Nagashima, and K. Sasaki, "Luminal A and luminal B (HER2
negative) subtypes of breast cancer consist of a mixture of tumors with
different genotype," *BMC research notes,* vol. 5, pp. 1-8, 2012.

\[74\] R. Yerushalmi, R. Woods, P. M. Ravdin, M. M. Hayes, and K. A.
Gelmon, "Ki67 in breast cancer: prognostic and predictive potential,"
*The Lancet Oncology,* vol. 11, no. 2, pp. 174-183, 2010.

\[75\] X. Tong, H. Gui, F. Jin, B. W. Heck, P. Lin, J. Ma, J. D.
Fondell, and C.-C. Tsai, "Ataxin-1 and Brother of ataxin-1 are
components of the Notch signaling pathway," *EMBO reports,* vol. 12, no.
5, pp. 428-435, 2011.

\[76\] Y. Wang, S.-G. Cho, X. Wu, S. Siwko, and M. Liu, "G-protein
coupled receptor 124 (GPR124) in endothelial cells regulates vascular
endothelial growth factor (VEGF)-induced tumor angiogenesis," *Current
molecular medicine,* vol. 14, no. 4, pp. 543-554, 2014.

\[77\] F. Chen, B. Han, Y. Meng, Y. Han, B. Liu, B. Zhang, Y. Chang, P.
Cao, Y. Fan, and K. Tan, "Ceruloplasmin correlates with immune
infiltration and serves as a prognostic biomarker in breast cancer,"
*Aging (Albany NY),* vol. 13, no. 16, pp. 20438, 2021.

\[78\] Z. Yu, Z. Wang, X. Yu, and Z. Zhang, "RNA-Seq-Based Breast Cancer
Subtypes Classification Using Machine Learning Approaches,"
*Computational Intelligence and Neuroscience,* vol. 2020, pp. 4737969,
2020/10/29, 2020.

\[79\] D. C. Schmitt, L. Madeira da Silva, W. Zhang, Z. Liu, R. Arora,
S. Lim, A. M. Schuler, S. McClellan, J. F. Andrews, A. G. Kahn, M. Zhou,
E. Y. Ahn, and M. Tan, "ErbB2-intronic MicroRNA-4728: a novel tumor
suppressor and antagonist of oncogenic MAPK signaling," *Cell Death &
Disease,* vol. 6, no. 5, pp. e1742-e1742, 2015/05/01, 2015.

\[80\] L. Yang, K. Higashisaka, Y. Haga, H. Tsujino, K. Nagano, and Y.
Tsutsumi, "Alpha-crystallin B chains enhance cell migration in
basal-like two triple-negative breast cancer cells," *Pharmazie,* vol.
77, no. 2, pp. 45-47, Feb 1, 2022.

\[81\] H. Zhang, Y.-z. Pan, M. Cheung, M. Cao, C. Yu, L. Chen, L. Zhan,
Z.-w. He, and C.-y. Sun, "LAMB3 mediates apoptotic, proliferative,
invasive, and metastatic behaviors in pancreatic cancer by regulating
the PI3K/Akt signaling pathway," *Cell Death & Disease,* vol. 10, no. 3,
pp. 230, 2019/03/08, 2019.

\[82\] K. Suijkerbuijk, P. Van Diest, and E. Van der Wall, "Improving
early breast cancer detection: focus on methylation," *Annals of
Oncology,* vol. 22, no. 1, pp. 24-29, 2011.

\[83\] K. E. Schuebel, W. Chen, L. Cope, S. C. Glöckner, H. Suzuki,
J.-M. Yi, T. A. Chan, L. V. Neste, W. V. Criekinge, S. v. d. Bosch, M.
van Engeland, A. H. Ting, K. Jair, W. Yu, M. Toyota, K. Imai, N. Ahuja,
J. G. Herman, and S. B. Baylin, "Comparing the DNA Hypermethylome with
Gene Mutations in Human Colorectal Cancer," *PLOS Genetics,* vol. 3, no.
9, pp. e157, 2007.

\[84\] T. Ushijima, and K. Asada, "Aberrant DNA methylation in contrast
with mutations," *Cancer Science,* vol. 101, no. 2, pp. 300-305, 2010.

\[85\] S. Kamalakaran, V. Varadan, H. E. G. Russnes, D. Levy, J.
Kendall, A. Janevski, M. Riggs, N. Banerjee, M. Synnestvedt, E.
Schlichting, R. Karesen, K. S. Prasada, H. Rotti, R. Rao, L. Rao, M. H.
E. Tang, K. Satyamoorthy, R. Lucito, M. Wigler, N. Dimitrova, B. Naume,
A. L. Borresen-Dale, and J. B. Hicks, "DNA methylation patterns in
luminal breast cancers differ from non-luminal subtypes and can identify
relapse risk independent of other clinical variables," *Molecular
Oncology,* vol. 5, no. 1, pp. 77-92, Feb, 2011.

\[86\] X. K. Ma, L. Yu, P. Z. Wang, and X. F. Yang, "Discovering DNA
methylation patterns for long non-coding RNAs associated with cancer
subtypes," *Computational Biology and Chemistry,* vol. 69, pp. 164-170,
Aug 2017.

\[87\] M. List, A.-C. Hauschild, Q. Tan, T. A. Kruse, J. Baumbach, and
R. Batra, "Classification of Breast Cancer Subtypes by combining Gene
Expression and DNA Methylation Data," *J Integr Bioinform,* vol. 11, no.
2, pp. 1-14, 2014.

\[88\] A. Thulin, C. Andersson, E. Werner Rönnerman, S. De Lara, C.
Chamalidou, A. Schoenfeld, A. Kovács, H. Fagman, F. Enlund, and B. K.
Linderholm, "Discordance of PIK3CA and TP53 mutations between breast
cancer brain metastases and matched primary tumors," *Scientific
Reports,* vol. 11, no. 1, pp. 23548, 2021/12/07, 2021.

\[89\] J. I. Herschkowitz, K. Simin, V. J. Weigman, I. Mikaelian, J.
Usary, Z. Hu, K. E. Rasmussen, L. P. Jones, S. Assefnia, S.
Chandrasekharan, M. G. Backlund, Y. Yin, A. I. Khramtsov, R. Bastein, J.
Quackenbush, R. I. Glazer, P. H. Brown, J. E. Green, L. Kopelovich, P.
A. Furth, J. P. Palazzo, O. I. Olopade, P. S. Bernard, G. A. Churchill,
T. Van Dyke, and C. M. Perou, "Identification of conserved gene
expression features between murine mammary carcinoma models and human
breast tumors," *Genome Biol,* vol. 8, no. 5, pp. R76, 2007.

\[90\] K. Dias, A. Dvorkin-Gheva, R. M. Hallett, Y. Wu, J. Hassell, G.
R. Pond, M. Levine, T. Whelan, and A. L. Bane, "Claudin-Low Breast
Cancer; Clinical & Pathological Characteristics," *PLoS One,* vol. 12,
no. 1, pp. e0168669, 2017.

\[91\] B. D. Lehmann, B. Jovanović, X. Chen, M. V. Estrada, K. N.
Johnson, Y. Shyr, H. L. Moses, M. E. Sanders, and J. A. Pietenpol,
"Refinement of Triple-Negative Breast Cancer Molecular Subtypes:
Implications for Neoadjuvant Chemotherapy Selection," *PLoS One,* vol.
11, no. 6, pp. e0157368, 2016.

\[92\] L. C. Xia, J. M. Bell, C. Wood-Bouwens, *et al.*,
\"Identification of large rearrangements in cancer genomes with barcode
linked reads,\" *Nucleic Acids Research*, vol. 46, no. 4, pp. e19-e19,
2018.

\[93\] J. M. Bell, B. T. Lau, S. U. Greer, *et al.*, \"Chromosome-scale
mega-haplotypes enable digital karyotyping of cancer aneuploidy,\"
*Nucleic Acids Research*, vol. 45, no. 19, pp. e162-e162, 2017.

\[94\] L. C. Xia, S. Sakshuwong, E. S. Hopmans, *et al.*, \"A
genome-wide approach for detecting novel insertion-deletion variants of
mid-range size,\" *Nucleic Acids Research*, vol. 44, no. 15, pp.
e126-e126, 2016.

**Xintong Chang** is currently an undergraduate senior majoring in
Mathematics and Applied Mathematics (Statistics) at South China
University of Technology.

**Jiemin Xie** is a graduate student at the School of Mathematics, South
China University of Technology.

**Hongyu Duan,** a Ph.D. candidate at South China University of
Technology, specializing in protein function prediction. Passionate
about AI applications in bioinformatics, he aims to leverage machine
learning to advance biological discovery and contribute to computational
biology innovation.

**Keyi Li** received his M.S. from Columbia University, NY 10027, USA
(February 2025). He specialized in statistics and bioinformatics.

**Xuemei Liu** received her PhD in Theoretical Physics from Sun Yat-sen
University, Guangzhou 510275, China (December 2008). She specialized in
computational modeling and statistical methods with applications in
bioinformatics and complex systems. She is an Associate Professor in the
School of Physics and Optoelectronics, South China University of
Technology, Guangzhou 510640, China.

**Yunhui Xiong**, professor, Ph.D. His main research interests include
image and video processing, geometric modeling and processing, 3D
reconstruction, 3D printing and biostatistics. He is an Associate
Professor in the School of Mathematics, South China University of
Technology, Guangzhou 510000, China.

**Xiangqi Bai** received a B.S. degree in Mathematics from Sichuan
University in 2015 and a Ph.D. in Systems Theory from the Academy of
Mathematics and Systems Science, Chinese Academy of Sciences, in 2021.
She is currently a postdoctoral researcher at the School of Medicine,
Stanford University. Her research focuses on developing computational
and statistical algorithms for analyzing complex biological data, with
applications in spatial data, cell-free DNA methylation, and single-cell
multi-omics datasets.

**Kaida Ning** received her PhD in Bioinformatics from the University of
Southern California, California 90007, USA (May 2020). She specialized
in statistical modeling for large-scale medical and biological data. She
is a research scientist at Peng Cheng Laboratory, Shenzhen 518055, China

**Li C. Xia** received his PhD in Bioinformatics from the University of
Southern California, California 90007, USA (May 2013). He specialized in
statistical modeling and algorithm development for biomedical data
sciences. He is a professor in the Department of Statistics and
Financial Mathematics, School of Mathematics, South China University of
Technology, Guangzhou 510000, China.

**Fig. 1.** **The overall** **study design, analysis workflow, and
dataset description:** (a) The principles of DNA-level subtype
classifier; (b) Distribution of the four intrinsic subtypes in the TCGA
(n = 931) and METABRIC (n = 1134) breast cancer datasets (c) The
conceptual design of UGES.
![Figure 1](https://github.com/labxscut/UGES/blob/main/Figures/Figure_1.png)

**Fig. 2. Evaluation of the classification strategies on DNA-level
features:** (a) The self-learned classification hierarchy
(B,(H,(LA,LB))); (b) Competitive classification hierarchies:
((B,H),(LA,LB)), (B,H,(LA,LB)), and (B,H,LA,LB); ROC and AUC of (c) all
strategies with full DNA features; (d) the (B,(H,(LA,LB))) strategy with
mutation, CNA, and methylation features.
![Figure 2](https://github.com/labxscut/UGES/blob/main/Figures/Figure_2.png)

**Fig. 3. UGES subtypes show improved clinical relevance:** (a) Alluvial
plot and (b) Confusion matrix between UGES and PAM50 subtyping for the
combined TCGA+METABRIC cohort; (c) Changes of p-values when switching
from PAM50 to UGES subtypes; (d) Forest plot of multivariate model
hazard ratios.
![Figure 3](https://github.com/labxscut/UGES/blob/main/Figures/Figure_3.png)

**Fig. 4. Heatmap of the 52 signature subtype-delineating alterations:**
from top to bottom are heatmaps of signature (a) Mutations, (b) Copy
Number Aberrations, and (c) Methylation features.
![Figure 4](https://github.com/labxscut/UGES/blob/main/Figures/Figure_4.png)

    This work was supported by the National Science Foundation of
    China (12571529 to L.C. Xia; 62472180 to Y. Xiong), the Guangdong
    Basic and Applied Basic Research Foundation (2022A1515-011426,
    2024A1515-010699 to L.C. Xia), and the Major Key Project of
    Pengcheng Laboratory (PCL2025AS212-3, PCL2024A02-2).

    Xintong Chang, Jiemin Xie, Hongyu Duan and Yunhui Xiong are at the
    School of Mathematics, South China University of Technology,
    Guangzhou 510000, Guangdong, China.

    Keyi Li is at the School of Professional Studies, Columbia
    University, NY 10027, USA

    Xuemei Liu is at the School of Physics and Optoelectronics, South
    China University of Technology, Guangzhou 510000, Guangdong, China.

    Xiangqi Bai is at the Division of Oncology, Department of Medicine,
    School of Medicine, Stanford University, CA 94305, USA.

    Kaida Ning is at the Peng Cheng Laboratory, Shenzhen, 518000,
    Guangdong, China.

    Li C. Xia is at the Department of Statistics and Financial
    Mathematics, School of Mathematics, South China University of
    Technology, Guangzhou 510000, China (e-mail: lcxia@scut.edu.cn).
