# effect-size-correlation
### Simulation study

#### Reminder for simulation v1

1. Issue about UKBB imputation:

   - check whether "raw genotypes" were imputed in "imputed genotypes"
     - see how imputation can be done via the same UKBB procedure
     - make sure all the wholes are filled in the "imputed genotypes"
   - for holes in "raw genotypes", can treat them as empty genotypes, impute and then retain the holes

2. true SNP effect size correlation—we coined the term "genetic effect correlation" to refer to it— is our target of estimation

3. Empirical correlation from GWAS effect size estimate should be a most direct estimate of "genetic correlation", but LD will hinge the estimation. One way is to disentangle **"GWAS effect size correlation estimate" from LD,** the other way is to use external biological knowledge or information to boost the power of estimation

#### Conception

1. How do we decompose the **"genetic effect correlation"**? Let's say two SNP has strong effect correlation, which means that they are functionally similar across all the traits: e.g. affecting the same gene, having similar function on cell type. Further understanding is that two SNPs not only show **co-causality** behavior, but also have similar pattern for directionality of effect——either beneficial or harmful simultaneously. It looks like a two step procedure.
2. **Genetic effect correlation** is SNP-wise, meaning that the effect size of two SNPs is assumed to be correlated for the same trait, instead of being treated IID. Its relationship to **genetic correlation** is that one is column wise correlation , the other is row-wise correlation of the effect size matrix $B_{q\times p}$ (q=trait number, p=SNP number)
3. There are two interpretations for **co-causality**: we need to construct a binary **true causality matrix $\pi_{q\times p}$** for q traits and p SNPs, and correspondingly a **co-causality matrix $C_{p\times p}$**. The first is that two columns, respectively referring to SNP 1 and SNP 2, has correlation, so that $C_{1,2}=\text{P(SNP 1 is causal, SNP B is causal)-P(SNP 1 is causal)P(SNP 2 is causal)}=cov(\pi_{.1},\pi_{.2})\neq 0$. The second is that two rows, respectively as trait A and trait B, has more overlapping level of causal SNPs than that is expected by random. Co-causality only means that two SNPs is functionally related, so that they are either causal simultaneously for this trait or non-causal for the other trait. The directionality of effect size is not yet under consideration. We can say that if we denote any non-zero effect size as one, then the "genetic effect correlation" is exactly the co-causality between two SNPs.

#### Tentative tasks

##### 1) Solving conditional expectation

Find a way to **estimate $E(\vec \beta_J^{(k)}|\vec{\hat{\beta}^{(k)}_J},\Sigma_J,LD_J,\Pi_J)$ for some trait k and all J SNPs**, which is our ultimate goal. Given genetic effect correlation matrix, LD matrix, GWAS effect size vector, conventional statistical method like **multivariate penalized regression** is an alternative, based on some nice assumption such as *joint normality*.

##### 2) Assuming posterior joint normality: estimating $\Sigma$ from GWAS summary statistics

Assuming joint normality of effect size, i.e. $\vec \beta_J |\pi_J\equiv1,\pi_J\in\Pi_J\sim MVN(\vec 0, \Sigma_J)$, given LD matrix, find the analytical solution for the conditional expectation. Also, try to solve the analytical relationship between $cor(\beta_1,\beta_2|\beta_1\neq0,\beta_2\neq0) =\Sigma_{1,2}$ and $cor(\hat\beta_1,\hat\beta_2)$, so as to estimate $\Sigma_J$ from GWAS effect size estimate across a large number of traits.

Substituting the MVN assumption to some multivariate point-normal assumption, by introducing a co-causality matrix  $C_{p\times p}$ . Find the analytical relationship between  $C_J$ and  $\Sigma_J$  between which a **nonzero/post-causality effect size correlation** may be needed, which is that $cor(\beta_1,\beta_2|\beta_1\neq0 \and \beta_2\neq0)=E(\beta_1\beta_2|\beta_1\neq0 \and \beta_2\neq0)-E(\beta_1|\beta_1\neq0 \and \beta_2\neq0)E(\beta_2|\beta_1\neq0 \and \beta_2\neq0)$

Then find the analytical solution for the conditional expectation. If unable, find a Bayesian solution and a MCMC implementation.

##### 3) Penalized regression to estimate B

Suppose we have used some penalized regression method like LASSO/Lassosum or regularized prior Bayesian method like FINEMAP to estimated the joint effect size. But two issues with these method are that: first, they are biased, second, they are underpowered, because their penalty/regularization/sparsity assumption tends to shrink small non-zero effects to zeros. **Therefore, we may believe that the non-zero effects given by the methods are highly credible, but zeros are not credible enough.** This scenario is similar to that "drop-out" problem in single-cell transcriptomics. 



##### 4) GCN to post-finemapping re-analysis

**Graph Convolution Network** is a learning algorithm designed for this kind of **semi-supervised** problem with some **"missing labels"**. GCN is a branch under Graph Neural Network. **GNN**, based on **graph embedding** that maps nodes into d-dimensional Euclidean Space $\R^n$ so that similar nodes have similar embeddings such as small cross product, can be applied to unsupervised, supervised, semi-supervised problem. It is popular in three fields of biology:

- single cell sequencing "imputation": 
  - features=binary single-cell transcription status matrix for all genes across all cell samples
  - known labels=a gene's 1s in the feature, unknown = a gene's 0s in the feature
  - Or, known labels=a cell's 1s in the feature, unknown = a cell's 0s in the feature
  - adjacent matrix: (from database) GTEx tissue specific co-expression, PPI
- gene expression (high or low) binary prediction: 
  - features=quantitative tissue expression level for all genes across all patients $N\times P$
  - known labels = a gene's binary expression level for this tissue (e.g. this gene is over/under-expressed in liver than other tissues) across part of the patients  $N_1\times1, N_1\leq N$; unknown labels =that are unmeasured for other patients
  - Or, known labels = a patient's binary expression level for this tissue across part of the genes  $P_1\times1, P_1<P$; unknown labels =that are unmeasured for other genes
  - adjacent matrix: (from database) GTEx tissue specific co-expression, PPI  $P \times P$
- gene priorization / disease association prediction:
  - features = quantitative expression level for all genes  $P\times N$
  - known labels = disease association profile for a subset of genes from database  $P_1\times q, P_1<P$
  - unknown labels = disease association profile for other genes compliment to the subset  $(P-P_1)\times q$
  - adjacent matrix: (from database) PPI, gene co-expression; (from data) correlation profile of the feature

**GCN may be a promising option in our problem, with $\hat B_{q\times p}$ from fine mapping methods' output as node feature for p SNPs, the columns with all 0s in $\hat B$ as nodes with missing labels, and a pre-specified adjacent matrix that is exactly or adjusted-form of our "genetic correlation matrix". **



##### 5) summary on ways to estimate $\Sigma$

**Functional similarity=Co-causality+Posterior effect size correlation**

Anyway, we will need to find a way to estimate $\Sigma_J$. 

- First, only by analytical solution using only $\hat B$ and $LD$ with some statistical assumption.
- Second, only by external information, such as **functional annotation** for all the SNPs, and then do correlation analysis, or logistic regression. The purpose is to estimating causality matrix $\Pi$.
- Third, some **weighted average** or regression integration of method of the first and the second approach's estimate.
- Fourth, as there are many gene-level interaction/functional similarity information available in database, find a way to estimate through the mediation of gene. Biologically, two SNPs located in the same gene must be co-causal (functional related), and if some genes are functionally related, we may say the SNPs within different genes are also like so. The problem will be on the SNPs that are non-genic, or can not be mapped to any gene.

Two more points are noteworthy. First, biologically, co-causality (functional relation) is easier to discover and estimate than functional similarity (genetic effect correlation). **Second, a fundamental assumption is that co-causality and genetic effect correlation is evenly distributed across all traits**; or rather, two SNPs shares a fixed co-causality/genetic correlation regardless of the trait category/domain we are investigating in. This is not a very strong (?) assumption as we are not saying two SNPs have the same magnitude of effect size across the same trait. Even with a fixed genetic effect correlation, two SNPs can still be simultaneously more related to one trait domain (larger $P(\beta_1\neq0 \and \beta_2\neq0)$ within this domain), and less related to another domain (smaller $P(\beta_1\neq0 \and \beta_2\neq0)$).

#### Lasso power test by simulation

##### A supervised scenario with perfect co-causality

###### Data generation scheme

sample size = **10000**

SNP number = **451** of chromosome one block one

using genotyped SNPs, some remaining holes imputed by median, then generate phenotype

causal SNP proportion = **20%**

causal SNP number = **93**

missing genotyped first imputed, then generate phenotype

**For one trait,**

$\vec y=X\vec \beta + \vec \epsilon$, 

$\vec \beta$ not i.i.d., but follows a joint point-normal distribution

with $\pi=0.01$, $1-\pi$ of them are zero, while the rest $\pi$ of them follows:

$\vec\beta \sim MVN(\vec 0,\Sigma) $

$\Sigma$ is the "genetic effect correlation" among the 1% causal SNPs.

$\epsilon \sim N(0, \frac{Var(X\beta)(1-h^2)}{h^2})$

$h^2=0.6$

$Var(X\beta)/Var(y)=1/(1+\frac{h^2}{(1-h^2)})=h^2$

**Replicate, and generate 100 traits, by setting the causal SNPs to perfectly (100%) overlap**

###### Lasso fine mapping

- Do lasso fine mapping, using CV+1lamda.se from *cv.glmnet*

- Obtain a $\hat B_{100\times 451}$

| N=N~SNP~*N~trait~=45100 | Null                            | Non-Null                               |       |                   |
| ----------------------- | ------------------------------- | -------------------------------------- | ----- | ----------------- |
| Exclude (Not Reject)    | TN=32783                        | FN=4836                                | 37619 |                   |
| Include (Reject H~0~)   | FP=3017                         | TP=4464                                | 7481  | **FD rate =0.40** |
|                         | 35800                           | 9300                                   |       |                   |
|                         | **FP rate (type1 error)=0.084** | FN rate (type 2) =0.52, **POWER=0.48** |       |                   |

**Conclusion: **

1. **the power is not bad in this scenario. FDR is not controlled.**
2. Knock-down method may be used to control FDR.

##### Using GCN to do post-fine mapping update

###### setting

- Training input: $\hat B_{73\times 450}$ as feature, SNP1's **non-zero** Lasso effect estimate' as label $\hat\beta_{73\times 1}$, a graph defined based on $\Sigma$

- Supervised learning, fitting model weights

- model:  the paper **Analysis of Gene Interaction Graphs as Prior Knowledge for Machine Learning Models** : [gene-graph-analysis](https://github.com/Bertinus/gene-graph-analysis)

  https://github.com/mila-iqia/gene-graph-conv

- Prediction input:  $\hat B_{27\times 450}$ 

- Prediction output: an updated label for SNP1's zero Lasso effect estimate $\hat\beta_{27\times 1}$

###### Compare the Power for SNP1 by GCN and simple Lasso

##### Univariate Lasso power spectrum

###### The same data generation scheme (independent trait, correlated SNP effect)

- alter $h^2$ and $\pi$ , therefore changing per-SNP heritability $h^2/(m\times \pi)$, m = 451 here.
- do univariate Lasso for each trait

###### Power against per SNP heritability spectrum



#### Analytical Formulation

Assume data

$\bold{Y_{n\times q}=X_{n\times p}B_{p\times q}+E_{n\times q}}$

$y_{nq\times 1}=vec(\bold{Y}) $

$\hat \beta^*(\bold{K})_{pq\times 1}$

\ is\ estimator\ for\ 

$\beta=vec(\bold{B})=(\beta_1^T,...,\beta^T_q)^T$

$\bold{K_{q\times q}}\succ0$



\ is\ the\ ridge\ matrix for $

##### 1) Basic Formulation under perfect co-causality and overlap

###### Basic assumptions in a  2 traits ✖ 2 SNPs setting by a linear mixed model

Here we actually assume that we only have perfectly-overlapping 2 causal SNPs for 2 traits, and we know their positions. 

**Assume a simple structure of 2 traits ✖ 2 SNPs** for n individuals

$(\vec y^{(1)},\vec y^{(2)})_{n\times 2}=X_{n\times 2}B_{2\times 2}+(\vec\epsilon^{(1)},\vec \epsilon^{(2)})_{n \times 2}$

For individual i,

$(y^{(1)},y^{(2)})^T=B_{2\times 2}(x_1,x_2)^T+(\epsilon^{(1)},\epsilon^{(2)})^T$

**Assuming perfect co-causality,** i.e. SNP 1 and 2 are causal/non-causal SNPs simultaneously.

**Assuming perfect genetic overlapping**, i.e. trait 1 and trait 2 have completely overlapping causal SNPs, i.e. SNP 1 and SNP2.

**Assuming** **random effect model**, i.e. $B=\begin{pmatrix}  \beta_1^{(1)} & \beta_1^{(2)} \\  \beta_2^{(1)}& \beta_2^{(2)} \end{pmatrix}=(\vec \beta^{(1)}, \vec\beta^{(2)})_{trait}=(\vec \beta_1,\vec \beta_2)^T_{SNP}$ is a random matrix: each column is a random vector for each trait; each row is also a random vector for a SNP. 

Marginally, for SNP-wise there is genetic effect correlation $\rho_s$, for trait wise, there is genetic correlation, $\rho_g$, i.e. **assuming marginally bivariate normal**,

trait-wise: **assuming two traits are equally heritable, i.e. identical $\sigma^2_g$**

$(\beta_j^{(1)},\beta_j^{(2)})^T \sim N(\vec0,\begin{pmatrix} \sigma^2_g & \rho_g \\  \rho_g& \sigma^2_g \end{pmatrix}),j=1,2$ ($\Sigma_g$)

SNP-wise

$(\beta_1^{(k)},\beta_2^{(k)})^T \sim N(\vec0,\begin{pmatrix} \sigma^2_g & \rho_s \\  \rho_s& \sigma^2_g \end{pmatrix}),j=1,2$ ($\Sigma_s$)

The complete but complicated genetic architecture, i.e. the joint distribution for random matrix B can be modelled by a copula function $C: [0,1]^2\rightarrow [0,1]$ s.t. coupling equation

$F(B)=C_1(F(\vec \beta_1),F(\vec \beta_2))=C_2(F(\vec \beta^{(1)}),F(\vec \beta^{(2)}))$



To further simplify, **assuming joint normal: different trait and different SNP, then no correlation in effect size**

unfinished

$(\beta_i^{(k)},\beta_j^{(l)})^T \sim N\{\vec0,\begin{pmatrix} \frac{h^2_k}{m_k} & \rho_{g_{kl}}\rho_{s_{ij}} (\frac{h^2_k}{m_k} \frac{h^2_l}{m_l})^{\frac{1}{2}} \\  \rho_{g_{kl}}\rho_{s_{ij}} (\frac{h^2_k}{m_k} \frac{h^2_l}{m_l})^{\frac{1}{2}} & \frac{h^2_l}{m_l} \end{pmatrix}\}$ 

$vec(\bold{B_{11:22}})=\begin{pmatrix}  \beta_1^{(1)} \\ \beta_1^{(2)} \\  \beta_2^{(1)}\\ \beta_2^{(2)} \end{pmatrix}\sim N(\vec 0,\Sigma), \Sigma=\begin{pmatrix} \frac{h^2_1}{m_1} & \rho_{g_{12}}(\frac{h^2_1}{m_1}\frac{h^2_2}{m_2})^{\frac{1}{2}}&\rho_{s_{12} } & 0 \\  \rho_{g_{12} }(\frac{h^2_1}{m_1}\frac{h^2_2}{m_2})^{\frac{1}{2}}& \frac{h^2_2}{m_2} & 0 &\rho_{s_{12} }\\ \rho_{s_{12} } & 0 &\frac{h^2_1}{m_1} & \rho_{g_{12} }\\ 0 &\rho_{s_{12} } & \rho_{g_{12} }&\frac{h^2_2}{m_2}\end{pmatrix}$

i.e.  for two effect size, $\Sigma_{..}=Cov(\beta_{j_1}^{(k_1)},\beta_{j_2}^{(k_2)})$, j for SNP index, k for trait index

$\Sigma_{..}=\sigma^2_g$ if $j_1=j_2,k_1=k_2$  (same SNP, same traits)

$\Sigma_{..}=\rho_s$ if $j_1\neq j_2,k_1=k_2$ (different SNPs, same traits)

$\Sigma_{..}=\rho_s$ if $j_1 = j_2,k_1\neq k_2$ (same SNPs, different traits)

$\Sigma_{..}=0$ if $j_1 \neq j_2,k_1\neq k_2$ (different SNPs, different traits)



**Further assuming no G✖E interaction, and no cross trait/cross individual environmental correlation** (like in GREML/LDSC model)

###### Sample covariance for the same trait (trait 1): same set of related individuals

$Cov(\vec y^{(1)},\vec y^{(1)}|X)_{n\times n}=Cov(X\vec\beta^{(1)},X\vec\beta^{(1)}|X)+\sigma^2_eI$

$=XCov(\vec\beta^{(1)},\vec\beta^{(1)})X^T+\sigma^2_eI$

$=X\Sigma_sX^T+\sigma^2_eI$

denote $X\Sigma_sX^T$ as $A$

For individual 1 and 2,

$A_{12}=\sigma^2_g\sum_{j=k} x_{1j} x_{2k} +\sum_{j\neq k}\rho_s^{(jk)} x_{1j}x_{2k}=m\sigma^2_g\overline{ x_{1j}x_{2j}}+m(m-1)\overline{\rho_s^{(jk)}x_{1j}x_{2k}}_{j\neq k}$

$=m^2\overline{\rho_s^{(jk)}x_{1j}x_{2k}}$  if we define $\rho_s^{(jk)}=\sigma^2_g$ for all $j=k$. 

(here j and k are both SNP index)

here, for our two SNP scenario, 

 $A=\begin{bmatrix} \sigma^2_g(x_{11}^2+x_{12}^2)+\rho_s(2x_{11}x_{12})&\sigma^2_g(x_{11}x_{21}+x_{12}x_{22})+\rho_s(x_{11}x_{22}+x_{12}x_{21})\\\sigma^2_g(x_{11}x_{21}+x_{12}x_{22})+\rho_s(x_{11}x_{22}+x_{12}x_{21}) & \sigma^2_g(x_{21}^2+x_{22}^2)+\rho_s(2x_{21}x_{22})\end{bmatrix}$

under conventional GRM model, $\rho_s^{(jk)}\equiv 0$

and genotypes have been SNP-wise normalized i.e. $E(x_{ij})=0,Var(x_{ij})=1$, so

$A=\sigma^2_g\begin{bmatrix} x_{11}^2+x_{12}^2&x_{11}x_{21}+x_{12}x_{22}\\x_{11}x_{21}+x_{12}x_{22}& x_{21}^2+x_{22}^2\end{bmatrix}$

$=m\sigma^2_g\begin{bmatrix} \overline x_{1.}^2&\overline{ x_{1.}x_{2.}}\\\overline{ x_{1.}x_{2.}}& \overline x_{2.}^2\end{bmatrix}=m\sigma^2_g G$

$ \overline x_{i.}^2\approx1$ if m is large enough.

$m=2 $（SNP number）

$G$ is GRM, usually approximated by an $p$ SNP-based matrix  $\hat G$.

$\hat A=p\frac{\hat h^2}{p}\hat G$  (note that we have$h^2=m\sigma^2_g$)

This approximation is considered to be accurate when the trait of polygenic enough, i.e. when the number of causal SNPs $m$ is large enough, and dispersed across the whole genome.

###### Sample covariance for different trait: : same set of related individuals

$Cov(\vec y^{(1)},\vec y^{(2)}|X)_{n\times n}=Cov(X\vec\beta^{(1)},X\vec\beta^{(2)}|X)+\sigma^2_eI$

$=XCov(\vec\beta^{(1)},\vec\beta^{(2)})X^T+\sigma^2_eI$

$=XE(\vec\beta^{(1)}\vec\beta^{(2)T})X^T+\sigma^2_eI$

$=XE\begin{pmatrix}  \beta_1^{(1)}\beta_1^{(2)} & \beta_1^{(1)}\beta_2^{(2)} \\  \beta_2^{(1)}\beta_1^{(2)}& \beta_2^{(1)}\beta_2^{(2)} \end{pmatrix}X^T+\sigma^2_eI$

denote $X\Sigma_sX^T$ as $S$

$S=X\begin{pmatrix}\rho_g & 0\\  0& \rho_g \end{pmatrix}X^T=m\rho_gG$

$m=2$

###### **Suppose we know $\rho_s$, it will provide a more accurate estimate of A and then a more accurate estimate of $h^2$ and $\rho_g$, under REML**

###### Extensing two SNPs to multiple SNPs

We have

 $Cov(\vec y^{(1)},\vec y^{(1)}|X)_{n\times n}=Cov(X\vec\beta^{(1)},X\vec\beta^{(1)}|X)+\sigma^2_eI$

$=XCov(\vec\beta^{(1)},\vec\beta^{(1)})X^T+\sigma^2_eI$

$=X\Sigma_sX^T+\sigma^2_eI$

denote $X\Sigma_sX^T$ as $A$

For individual 1 and 2,

$A_{12}=\sigma^2_g\sum_{j=k} x_{1j} x_{2k} +\sum_{j\neq k}\rho_s^{(jk)} x_{1j}x_{2k}=m\sigma^2_g\overline{ x_{1j}x_{2j}}+m(m-1)\overline{\rho_s^{(jk)}x_{1j}x_{2k}}_{j\neq k}$

$=m^2\overline{\rho_s^{(jk)}x_{1j}x_{2k}}$  if we define $\rho_s^{(jk)}=\sigma^2_g$ for all $j=k$. 

But now m is much larger than 2.

$\overline{ x_{1j}x_{2j}}$ is still average IBD proportion of two genomes at causal SNPs, which can be approximated by SNP data.

In GCTA, suppose the we have genome-wide SNP genotype data for n individual as $X_{n\times p}$

Estimate GRM $G$ by $S_{n\times n}=XX^T/n$

where  $x_{ij}=(w_{ij}-2p_j)/\sqrt{2p_j(1-p_j)}$, the standardized genotype

$S=\sum_j S_j$, $S_j$ is a n by n matrix for SNP j.

Genotypes have been SNP-wise normalized i.e. given j, $E(x_{ij})=0,Var(x_{ij})=1$, 

Assuming polygenicity, then

$\overline{ x_{1j}x_{2j}}\approx s_{12}$, which is the 1-2 element in $S$.

The **key** is to estimate $\overline{\rho_s^{(jk)}x_{1j}x_{2k}}_{j\neq k}$ and  $\overline{\rho_s^{(jk)}x_{1j}x_{1k}}_{j\neq k}$.

There are in total m(m-1)/2 $\rho_s^{(jk)}$, so when m is large, it is close to infinity, we can assume it follows a continuous distribution:

**Assuming a prior $\rho_s^{(jk)}\sim N(\mu_s,\sigma^2_s) $ for different j and k**

**Further assuming this distribution is independent from that of the LD of genotypes (key assumption)**, 

$\rho_s^{(jk)}\perp r_{jk}$, so

$\overline{\rho_s^{(jk)}x_{1j}x_{1k}}_{j\neq k}\approx E(\rho_s^{(jk)})E(x_{1j}x_{1k})\approx\mu_s\mu_r$, where $\mu_r$ is the average pair-wise LD, which could be zero.

$\overline{\rho_s^{(jk)}x_{1j}x_{2k}}_{j\neq k}\approx E(\rho_s^{(jk)})E(x_{1j}x_{2k})\approx\mu_s\mu_rs_{12}$,

Thus the new covariance matrix

$A\approx[m\sigma^2_g +m(m-1)\mu_s\mu_r]S$

$S$ is still the SNP approximated GRM, and we can calculate $\mu_r$ from reference LD pa444nel, so the unknown parameters are

$\sigma^2_g$, $m$, $\mu_s$, which can be estimated under REML.

###### Linking genetic effect correlation to LD?

Is it possible to estimate $\rho_s^{(jk)}$ individually under a fixed model, rather than globally estimate $\hat\rho_s$?

It is only possible when we have multiple trait information, like GWAS for multiple correlated traits.

##### Modelling co-causality and genetic overlap

Denote the non-zero entries of $B$ as 1, then convert $B_{p\times r}$ into a 0-1 causality matrix $\Pi_{p\times r}$ , which should be very sparse. (p SNPs, r traits)

$Y_{n\times r}=X_{n\times p}B_{p \times r}+E$

$\pi_{jk}=I(\beta_{jk}\neq0)$

We define:

- SNP-wise **co-causality** $\rho_c$ as the row-wise correlation of $\Pi_{p\times r}$

we have $\pi_{jk}=\pi_{jk}^2$ as 1^2^=1,0^2^=0, for SNP 1 and SNP2

$\rho_c^{12}=\frac{P(\pi_1=1,\pi_2=1)-E\pi_1E\pi_2}{(E\pi_1-E^2\pi_1)(E\pi_2-E^2\pi_2)}=\frac{P(\pi_1=1,\pi_2=1)-E\pi_1E\pi_2}{\sqrt{E\pi_1E\pi_2(1-E\pi_1)(1-E\pi_2)}}$

 empirically,

  $\rho_c^{12}\approx\frac{\sum_{k=1}^r\pi_{1k}\pi_{2k}/r-\sum_{k=1}^r\pi_{1k}\sum_{k=1}^r\pi_{2k}/r^2}{\sqrt{\sum_{k=1}^r(\pi_{1k}-\overline \pi_{1k})^2/r\times\sum_{k=1}^r(\pi_{2k}-\overline \pi_{2k})/r}}$ 

It is very easy to calculate. Supposing there are 100 traits, SNP1 and SNP2 are both causal SNPs for the first 5 traits, and non-causal for the rest. Then

$\rho_c=(0.05-0.05*0.05)/(0.05^2*0.95^2)=1$, which is **perfect co-causality.**

**If two SNPs are functionally independent (in terms of causality), $P(\pi_1=1,\pi_2=1)=E\pi_1E\pi_2=P(\pi_1=1)P(\pi_2=1)$**

then $\rho_c^{12}=0$

- trait-wise **genetic overlapping** as the column-wise correlation of $\Pi_{p\times r}$

Similarly, we have 

$\rho_o^{(12)}=\frac{P(\pi_{(1)}=1,\pi_{(2)}=1)-E\pi_{(1)}E\pi_{(2)}}{(E\pi_{(1)}-E^2\pi_{(1)})(E\pi_{(2)}-E^2\pi_{(2)})}=\frac{P(\pi_{(1)}=1,\pi_{(2)}=1)-E\pi_{(1)}E\pi_{(2)}}{\sqrt{E\pi_{(1)}E\pi_{(2)}(1-E\pi_{(1)})(1-E\pi_{(2)})}}$

##### Linking things together



for a unspecified trait, we model SNP1 and SNP2's joint effect size as ($\beta_1$, $\beta_2$)

and their causality status $(\pi_1,\pi_2)$,also as a random vector.

Assuming marginal effect size distribution for one SNP is of point-normal shape regardless of the trait,

$E(\beta_1, \beta_2)^T=(0,0)^T$,  $E(\pi_1,\pi_2)^T=(p_1,p_2)$ which measures the **pleiotropy** of SNPs.

$Var\beta_1=E \beta_1^2=E(E(\beta_1^2|\pi_1))=E(\overline h^2\pi_1)=p_1\overline {h^2/m}$

$\overline {h^2/m}=average\ h^2_{SNP}\ for\ all\ traits$ denoted as $h^2_m$

For trait k, per causal SNP heritability $h^2_{SNP(k)}=h_{(k)}^2/m=h_{(k)}^2/(N_p*p_{(k)})$

**heritability** divided by **the number of causal SNPs (=No. of SNPs*polygenicity)**, which measure the degree of micro-effect, i.e.**"infinitesmality"**

$\beta_1\beta_2|(\pi_1\pi_2=0)\equiv0$

denote $E\beta_1\beta_2|(\pi_1\pi_2=1)=\rho_s$

then $E\beta_1\beta_2|\pi_1\pi_2=\rho_s\pi_1\pi_2$

$E\beta_1\beta_2=E(E\beta_1\beta_2|\pi_1\pi_2)=\rho_sE\pi_1\pi_2=\rho_sp_{12}$

Denote $p_{12}=P(\pi_1=1,\pi_2=1)=E\pi_1\pi_2$

**Suppose** we know from biology functional information that these two SNPs are functionally closely related, like within the same gene region, then we may determine $\rho_c$ to be close to 1. More generally, we may regress SNP1's functional annotation on SNP2's functional annotation, the standardized coefficient might be an estimate of $\rho_c$. 

$\hat \rho_c^{12}=cor(C_1,C_2)$, $C_1$ and $C_2$ are dummy variables for SNP1 and SNP2, encoding their functional annotations.

If we also know the pleiotropy level i.e. $p_1$ and $p_2$ for these two SNPs, then we can solve $p_{12}$ from the formula above.

$(\beta_1, \beta_2)^T\sim \{\vec0,\begin{pmatrix}p_1\overline {h^2_m} & p_{12}\rho_s \\ p_{12}\rho_s & p_2\overline{h^2_m}\end{pmatrix}\}$

#### Algorithmic realization: penalized regression

#####  notations

<img src="C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416120405197.png" alt="image-20200416120405197" style="zoom:67%;" />

##### **Options**

##### 1) OLS per trait

Problems

- ignoring other responses (potentially useful information)
- high LD (co-linearity)

<img src="C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416120646306.png" alt="image-20200416120646306" style="zoom:50%;" />

##### 2) Ordinary Lasso: per trait

>  Robert Tibshirani. Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society. Series B (Methodological), pages 267–288, 1996. 
>
>  Efron et al. (2004) proposed least angle regression selection (LARS) 

- penality:

  for trait l,

$H=\alpha \sum_{j=1}^p   |\beta_{jl}|$

do it trait by trait

- realization: R package *glmnet* or *lars*

##### 3) Multivariate Lasso for all traits

- penalty form

$H=\alpha \sum_{j=1}^p \sum_{k = 1}^r \xi_{jk} |\beta_{jk}|$

let $\xi_{jk}\equiv1$

$H=\alpha \| \boldsymbol{B} \|_1$

- realization

  R package *lsgl* : set $\alpha=1$

##### 4) Group Lasso per trait： grouping features

 Yuan & Lin (2007) http://www.columbia.edu/~my2550/papers/glasso.final.pdf

>  M. Yuan and Y. Lin. Model selection and estimation in regression with grouped variables. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1):49–67, 2006.

- general penalty: J groups of features

  ![image-20200416230339788](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200416230339788.png)

   X~j~ is an n×p~j~ matrix corresponding to the jth factor and β~j~ is a coefﬁcient vector of size p~j~, j=1,...,J. 

- This is intermediate between the l1-penalty that is used in the lasso and the l2-penalty that is used in ridge regression

- most common case

$ p_1=...=p_J =1$

it is reduced to normal Lasso

- most used norm (that is used in this paper): $K_j=p_jI_{p_j}$, a weighted L-2 norm for each predictor group

<img src="C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416132028394.png" alt="image-20200416132028394" style="zoom:50%;" />

$\beta_l$ is a vector for group $l$ of predictors

This procedure acts like the lasso at the group level: depending on λ, an entire group of predictors may drop out of the model. 

- realization: R package *gglasso* or *grplasso*

##### 5) Sparse Group Lasso per trait:  grouping features

https://statweb.stanford.edu/~tibs/ftp/sparse-grlasso.pdf

- penalty: L groups of features

<img src="C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416132142604.png" alt="image-20200416132142604" style="zoom:50%;" />

 β = (β~1~,β~2~,...β~l~) is the entire parameter vector. For notational simplicity we omit the weights $\sqrt p_l$ 

- realization: R package *SGL*

##### 6) Multivariate Group Lasso: selecting a group of features

> K. Lounici, M. Pontil, A. B. Tsybakov, and S. van de Geer. Taking advantage of sparsity in multi-task learning. arXiv preprint arXiv:0903.1468, 2009.

- <img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200416233954876.png" alt="image-20200416233954876" style="zoom:50%;" />
- penalty: M features, T tasks
  - Assuming the set of relevant predictor variables is the same across the different equations (*structured sparsity assumption*)
  - data matrix for a task can be different
  - don't define groups: $p_1=...=p_M=1$
  - general mixed L-2,p norm

<img src="C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416135135454.png" alt="image-20200416135135454" style="zoom:50%;" />

- used in the paper: mixed L-2,1 norm

![image-20200416135306505](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200416135306505-1587048586087.png)

- extension: **grouping features**

  there are M groups, $p_1+...+p_M=p$

  we still used L-2,1 norm, but $\beta_j$ will be a T by p~j~ matrix rather than a vector

- **It results in the selection of the same predictors across all responses.  If a feature or a group of features is selected for one response, then it is selected for all responses. **

- realization: R package *lsgl* 

with penalties of the form

$H=\alpha \sum_{j=1}^p \sum_{k = 1}^r \xi_{jk} |\beta_{jk}| + (1-\alpha) \sum_{j = 1}^p \gamma_{j} \| \boldsymbol{\beta}_j \|_2$

let $\xi_{jk}\equiv1$

By default, $\gamma_{j} = \sqrt{r}$

$H=\alpha \| \boldsymbol{B} \|_{1,1}+(1-\alpha) \sum_{j = 1}^p \sqrt r \| \boldsymbol{\beta}_j \|_2$

let $\alpha=0$.

##### 7) Multivariate Sparse Group Lasso: Selecting arbitrary group of parameters

> Y. Li, B. Nan and J. Zhu (2015) Multivariate sparse group lasso for the multivariate multiple linear regression with an arbitrary group structure. Biometrics. DOI: 10.1111/biom.12292

- <img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200416234515667.png" alt="image-20200416234515667" style="zoom:50%;" />

Y column groups, X row groups, intersection block groups

- penalty

<img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200416234850350.png" alt="image-20200416234850350" style="zoom: 50%;" />

 the union of all groups inG does not need to contain all the elements of B, 

<img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417000729296.png" alt="image-20200417000729296" style="zoom:50%;" />

 For notational simplicity, special case is λ~jk~ = λ for all j, k.

- note about matrix norm:

  1. mixed norm/ **L-p,q norm**

  ![image-20200417011929795](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417011929795.png)

  the most used one is L-2,1 norm

  ![image-20200417012048855](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417012048855.png)

  also the L-2,2 norm, i.e. Frobenius norm

  ![image-20200417012231523](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417012231523.png)

  2. **induced norm**

  ![image-20200417012342684](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417012342684.png)

  the most used is spectrum norm, p-2 induced norm,

![image-20200417012513642](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417012513642.png)

where $A^*=\overline A^T$, the conjugate matrix's transpose.

normally we have ![image-20200417012614818](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417012614818.png)	

   3. **"Entrywise" norms** 

      ![image-20200417012832227](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417012832227.png)

      **here the author used $\|B_g\|_2$ as the p-2 "entrywise" norm**, i.e. L-2,2 norm for the block $B_g$; but in general it is still kind of L-2,1 norm for $B$, because of the sum across all groups. ("L2 plus sum" in the vector scenario)

      

- realization: R package *MSGLasso*

  



##### 8) Sparse Overlapping Sets Lasso (SoSLasso): less rigid than MGL, grouping similar features

**Suppose that the available features can be organized into overlapping subsets according to a notion of similarity**, and that the features useful in one task are similar, but not necessarily identical, to those best suited for other tasks. **In other words, a feature that is useful for one task suggests that the subset it belongs to may contain the features useful in other tasks** 

- need to define M groups of features according to similarity

  with a data matrix Φ~t~ ∈ R^n×p^ for each task t ∈ {1,2,...,T}. Can be the same one

  X here is the effect size matrix

- Penality

![image-20200416130612890](C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416130612890.png)

- realization: Matlab https://github.com/UWKCL/soslasso

##### 9) Multivariate Cluster Elastic Net (MCEN): learning to cluster responses

- defining groups/cluster of responses and permuting, to find best cluster
- within a group of traits, effect sizes are shared

<img src="C:\Users\22367\Desktop\fine mapping simulation\Simulation study v2.assets\image-20200416123441883.png" alt="image-20200416123441883" style="zoom:50%;" />

 Q, the total number of clusters

Let D = (D~1~,...,D~Q~) be a partition of the set{1,...,r}. (r traits)



##### 9) Curds and Whey method

 Breiman and Friedman (1997) 



##### Three strategies:

1. Do not define group structure, but use data-driven ones (like in MCEN)
2. Define trait+SNP group structure, by clustering traits based on genetic correlation, and clustering SNPs based on similarity of functional annotations
3. Do not even use group structure



#### Defining group structure on causality matrix

##### **Perfect genetic overlapping+perfect co-causality**

- $\rho_o$=1,$|\rho_c|=1$ : all traits share the same set of causal SNPs

  **Actually here the |$\rho$| can not be calculated within groups, as there is no variance. There is just perfect covariance.**

  Fig.a

<img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417033120333.png" alt="image-20200417033120333" style="zoom:33%;" />

- special case: causal SNP sets are bilaterally exclusive, $|\rho_o|=1$

  <img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417063927325.png" alt="image-20200417063927325" style="zoom:33%;" />



##### within-group perfect g.o.+c.o., between group imperfect

- 3 by 3

<img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417062712287.png" alt="image-20200417062712287" style="zoom:33%;" />

- 4 by 4

  <img src="C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417064944601.png" alt="image-20200417064944601" style="zoom: 33%;" />

  **can be arbitrary G~s~ by G~T~  with "overlapping"**

  A group is a row-column block-shaped intersection

#### Multivariate Penalized Regression Simulation

##### Basic setting

###### Genotypes

sample size **n** = **1000**

SNP number **p** = **451** of chromosome one block one

using genotyped SNPs, some remaining holes imputed by median, then generate phenotype

###### Genetic architecture

**Assuming perfect genetic overlapping** like fig.a

causal SNP proportion $\pi_c$

causal SNP number $m=ceiling[451*\pi_c]$

heritability interval $h^2\in I=$[0,0.25) ; [0.25, 0.5); [0.5, 0.75); [0.75, 1)

**In each simulation, select a interval for heritability→uniform random number→heritability for different trait **

###### Generating phenotypes

We follow the I-IV scenarios to generate joint effect size matrix $B$, and then assume a linear model

$Y_{n\times r}=XB_{p\times r}+E_{n\times r}$

trait number **r=100**

$E$ is a random residual error matrix, for each column k corresponding to trait k,

$\epsilon_i^{(k)} i.i.d.  N(0, \frac{Var(X\beta)(1-h^{(k)^2})}{h^{(k)2}}),i\in\N^n$

$h^{(k)2}\sim Uniform(I)$

Note: For each k, $Var(X\beta)/Var(y)=1/(1+\frac{h^2}{(1-h^2)})=h^2$

###### Plot power spectrum against per causal SNP heritability for different J kinds of penalized regression methods

**changing $\pi_c$ in each simulation, so after an simulation, we get an empirical power for J methods**

For a fixed interval of heritability, change $\pi_c$ 10 times, do 10 simulations, get 10 points for one method.

**Using mean(h^2^)/$\pi_C$ as x-axis.**

**Plot J methods on one plot, with J ✖ 10 points and J curves**

**Make four plots, annotated by average h^2^**

###### Four methods

1. Standard Lasso per trait
2. Multivariate Lasso for all traits
3. Multivariate Group Lasso: traits as one group, each SNP as a group--grouping SNPs, i.e. removing or selecting a SNP for all traits at the same time
4. Sparse Multivariate Group Lasso: traits as one group, each SNP as a group--grouping SNPs, but adding some element-wise sparsity

##### Scenario I:independent trait and independent SNP effect size

###### Effect size generation

generating each element of  one by one, following normal distribution N(0,h^(k)2^/m)

###### high $h^2\in [0.75,0.999]$

1) Standard Lasso

![image-20200417153638735](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417153638735.png)

![image-20200417153726813](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417153726813.png)



4) MSGLasso

![image-20200417154002746](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417154002746.png)



![image-20200417153822759](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417153822759.png)

###### low $h^2 \in [0.001,0.25]$

1) Standard Lasso

![](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417154918796.png)

![image-20200417154834836](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417154834836.png)

2) MLasso

![image-20200417160847017](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417160847017.png)

![image-20200417160909337](C:\Users\22367\Documents\WeChat Files\lane13318032043\FileStorage\File\2020-04\v2(1)(1).assets\image-20200417160909337.png)







##### Scenario II: only trait-wise genetic correlation



##### Scenario III: only SNP-wise effect correlation

##### Scenario IV: both trait-wise genetic correlation and SNP-wise effect correlation

##### Scenario V: imperfect group-wise co-causality & genetic overlapping without "overlapping"

##### Scenario VI: arbitrary causality structure with "overlapping"



### Real data analysis v.1

#### Implications from simulation analysis

1. Going ahead to UKBB real data: need to analyze all individuals
2. Chromosome by chromosome is better than LD block by block
3. Possible big data processing methods:
   - big data module in R
   - C++ lang
   - recursive incremental Bayesian learning
   - recursive Lasso
   - parallel running on server
4. Multivariate penalized regression should be better than the univariate one, so the hyperparameter tuning is critical and challenging: CV, BIC, AIC, quite many criteria. After fine tuning, it should have better performance. Our goal is not to compare different methods under different simulation scenarios, but to determine whether a method works for certain scenario, and most importantly, works for the real data.
5. L2,0 norm may be used in replace of L2,1 norm in sparse group lasso.
6. It is very necessary to review the multivariate analysis in genetics (like PHEN, M-SKAT, ACAT) and the application of penalized regression in genetics.
7. We can specify SNP group by functional annotations clustering and traits groups by genetic correlation clustering. Or most extreme: each/all as a group 
8. Using different $\lambda$ for different traits indicates that they have different degrees of sparsity, i.e. polygenicity.
9. Computational intensity is a concern for Bayesian methods. So is it possible to formulate some Bayesian methods (like point-mixture normal) and then develop multivariate penalized regression optimization methods for it?
10. Cross-validation may cause problem in lasso solution, which may have good prediction result but bad variable selection result, deviated from the purpose of fine mapping. BIC etc. may be good alternatives.
11. The idea of a 2 stage procedure where stage 1 uses multivariate ridge regression to combine information across GWAS on multiple correlated phenotypes, and stage 2 used thresholding to introduce sparsity.
12. Maybe do multivariate Lasso in a pairwise/bivariate manner (iteratively) for a whole chromosome, rather than in a Block-wise manner? 
13. Comparison should be made based on the same baseline, like using the same FWER/FDR control criterion. Then it is meaningful to compare the power and precision of different methods.
14. A mimicking of Knockoff method to detect significant results: simulaye

#### Computation speed analysis

##### Real genotypes+ simulated phenotypes

- 4cpu+**100gb** running memory
- *lsgl* R package, no parallel
- hyperparameter per fitting 

$\gamma$ = 0 ---- lasso

0<$\gamma$ < 1 ---- sparse group lasso

$\gamma$ = 1 ---- group lasso

$\lambda$ shrinkage factor $\in (0,1)$ given standardized data

| $(\gamma,\lambda)$ | \# SNP         | \# individuals | \# traits | data loading time (min) | approximate running time  per fitting (min) |
| ------------------ | -------------- | -------------- | --------- | ----------------------- | ------------------------------------------- |
| 0.5,0.1            | 450 (block1)   | 200,000        | 10        | <1                      | ~15                                         |
| 0.5,0.1            | 450            | 200,000        | 50        | <1                      | NA                                          |
| 0.5,0.1            | 450            | 200,000        | 100       | <1                      | NA                                          |
| 0,0.1              | 12113  (chr22) | 200,000        | 1         | ~10                     | ~5                                          |
| 0.5,0.1            | 12113  (chr22) | 200,000        | 10        | ~10                     | NA                                          |

**From the preliminary analysis, it seems univariate Lasso can only handle ~10 traits simultaneously for several hundred of SNPs. Or, fewer traits, but more SNPs**. In terms of CV, it cannot even handle 200k individual + one chromosome using CV-tuned $\lambda$ due to memory issue. **So possibly $\lambda$ has to be pre-determined before fitting, or we have to use CV on a smaller subset of original dataset (N~ 200k).** 

**The computation difficulty of CV was also indicated in a published paper (Cule et al 2012).**

#### Real data analysis on variable selection selection

##### univariate case

######  **Real genotypes+ real phenotypes**

**Genotypes**: Chromosome 22 (12113 SNPs) for 200k individuals

**10 Phenotypes**：

- **Quantitative**:  height, BMI, adjusted systolic blood pressure, adjusted diastolic blood pressure,

- **Binary**: diabetes, heart attack, stroke, high BP,  angina, 

- **Ordinal** heart problems

**Some of them are highly correlated**

 Case ratio for the 5 binary traits: 

[1] 0.048

[2] 0.032

[3] 0.023

[4] 0.015

[5] 0.274



###### **1) How to do univariate Lasso for 20k individuals + Chr22?**

We cannot use CV directly. **We use CV on a smaller sample size like 20k, and then we get the optimal lambda.** Here we use lambda.1se (largest value of lambda such that error is within 1 standard error of the minimum MSE), or we can use lambda.min.   And then we fixed a Lasso on this single lambda, on a large sample size like 200k.

observation: lambda.min is always the largest one among all lambdas in CV

**The number of variables selected by Lasso depends on $\lambda$, i.e. the strength of penalty (enforced sparsity).**

**Treating binary variable as quantitative ones (using Gaussian instead of Poisson/binomial model in *glmnet*) renders computation plausible.**

###### **2) ** **Elastic net**:

**Also do CV on a smaller sample size of 20k**

ridge--- 0<$\alpha$<1 ---lasso

less penalty, more variables selected

###### **3) ridge regression:**

**We cannot even do CV on 20k individuals for ridge. At this stage, we just arbitrarily give one.**

How about effect size thresholding？

**a) thresholding on |$\beta$|>0.01 --- variance explained by this SNP at least 10^-4^=0.01^2^ ~~ 1% of 10M SNPs = 100k causal SNPs ----   ~ 1/100k = 10^-5^  each SNP's variance explained**

**b) thresholding on |$\beta$|>0.0033--- variance explained by this SNP ~ 10^-5^=0.0033^2^** 

**c) thresholding on |$\beta$|>0.001--- variance explained by this SNP ~ 10^-6^=0.001^2^** 

glmnet did not work for ridge in terms of CV on 20k samples

###### 4） How about Bayesian way?

<img src="C:\Users\22367\Documents\markdown\image-20200522140151220.png" alt="image-20200522140151220" style="zoom:50%;" />

##### Results: No. of variables selected

| **Traits**                        | per SNP GWAS (p value<5*10^-8^) without LD clumping | per SNP GWAS  with LD clumping (r^2^ = 0.2) | per SNP GWAS (p value<10^-7^) without LD clumping | Univariate Lasso (two-step) | Univariate Elastic Net ($\alpha=0.5$) | Univariate Ridge (with thresholding $\beta^2$>=e-4,e-5,e-6 ,$\lambda=0.1$) |
| --------------------------------- | --------------------------------------------------- | ------------------------------------------- | ------------------------------------------------- | --------------------------- | ------------------------------------- | ------------------------------------------------------------ |
| height                            | 2                                                   | 2                                           | 4                                                 | 13                          | 57                                    | 211<br/>250<br/> 308                                         |
| BMI                               | 0                                                   | 0                                           | 2                                                 | 15                          | 45                                    | 260<br/>325<br/>410                                          |
| adjusted systolic blood pressure  | 0                                                   | 0                                           | 0                                                 | 12                          | 58                                    | 186<br/>125<br/>260                                          |
| adjusted diastolic blood pressure | 2                                                   | **1**                                       | 3                                                 | 13                          | 47                                    | 220<br/>246<br/>276                                          |
| diabetes                          | 0                                                   | 0                                           | 0                                                 | 3                           | **170**                               | 0<br/>17<br/>256                                             |
| heart attack                      | 2                                                   | 2                                           | 3                                                 | 2                           | **176**                               | 0<br/>20<br/>333                                             |
| stroke                            | 1                                                   | 1                                           | 1                                                 | NA                          | NA                                    | NA                                                           |
| high BP                           | 5                                                   | 5                                           | 5                                                 | NA                          | NA                                    | NA                                                           |
| angina                            | 0                                                   | 0                                           | 0                                                 | NA                          | NA                                    | NA                                                           |
| heart problems                    | 0                                                   | 0                                           | 0                                                 | 2                           | **211**                               | 149<br/>383<br/>587                                          |

Note: NA due to some coding error?

##### Predictive accuracy

Training=160k, Testing = 40k, other settings are the same.

**$R^2(\hat y,y)$**

| **Traits**                       | Univariate Lasso (two-step) | Univariate Elastic Net ($\alpha=0.5$) | Univariate Ridge (with thresholding) |
| -------------------------------- | --------------------------- | ------------------------------------- | ------------------------------------ |
| height                           |                             |                                       |                                      |
| BMI                              |                             |                                       |                                      |
| adjusted systolic blood pressure |                             |                                       |                                      |



##### Do these selected SNPs overlapping with each other?

Let's look at diabetes and heart attack:

**1. Are the smaller subset of selected SNPs included in larger subset of selected SNPs?**

The 3 and 2 SNPs selected for diabetes and heart attack  by Lasso (CV) **are also included** in the 256 SNPs selected by Ridge ($\lambda =0.2,\beta^2\geq 10^{-5}$) and 70 SNPs by Elastic Net (CV).

Using $\lambda=0.2$, for diabetes, for three thresholds, N~selected~ = 0/3/253  (vs 0/17/256), now the 3 SNPs selected by Lasso are not overlapping with the "3", but are included in the "253". 

【0/3/253 corresponds to $\beta^2$ = 0.0001, 0.00001, 0.000001 (e-6)】

**2. However, the 3 Lasso SNPs are not even included in those selected by $\beta^2$>=2e-6**

**That means the Lasso SNPs only have a moderate level of effect size magnitude. **

**3. Are the Lasso 3 and Ridge 3 have very small p-value in marginal testing? i.e. **Are the GWAS significant SNPs also selected by penalized regressions?

 **No!** Their p-values are quite large. That means the GWAS hits are not selected by Lasso nor Ridge.

#### Implementation of univariate ridge regression

##### Testing and p-value in ridge

###### theory: 'non-exact' *t*-type test

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-372



![image-20200429140514309](C:\Users\22367\Documents\markdown\image-20200429140514309.png)

Distribution of it under homogenous normal error (GLM assumption):

<img src="C:\Users\22367\Documents\markdown\image-20200529124415712.png" alt="image-20200529124415712" style="zoom:67%;" />

Correspondingly, a T statistic

<img src="C:\Users\22367\Documents\markdown\image-20200529091140592.png" alt="image-20200529091140592" style="zoom:50%;" />

<img src="C:\Users\22367\Documents\markdown\image-20200529091302690.png" alt="image-20200529091302690" style="zoom:50%;" />

Estimate of variance of residual error:

<img src="C:\Users\22367\Documents\markdown\image-20200529091327855.png" alt="image-20200529091327855" style="zoom:50%;" />

$v$ is the residual effective degrees of freedom given in Hastie & Tibshirani [[1990 paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-372#ref-CR18)]

$v=n-tr(2\bold H-\bold{HH'})=n-tr(\bold H)$

where <img src="C:\Users\22367\Documents\markdown\image-20200529123949398.png" alt="image-20200529123949398" style="zoom:50%;" />

under $H_0, T_λ\sim t_{v} $, when **v** is large enough, this approximates $N(0,1)$.

when $\lambda=0$, $tr(\bold H)=p$, all is reduced to OLS.



<img src="C:\Users\22367\Documents\markdown\image-20200529091750420.png" alt="image-20200529091750420" style="zoom:50%;" />

###### Real data application of testing

**Main challenge: algebra of large matrix**

1. Background setting

- R package *lmRidge* for fitting

- n=10k, p=1000 (from chr22) , 965 remained after MAF filtering 0.01% (exclude SNPs that are with too small MAF), phenotype=height

- lambda=0.1 (arbitrarily given)

- LD is troublesome for this package. Need to do some kind of LD clumping~~

  After that for all remaining SNPs LD<=**threshold=0.9**, and 884 SNPs pass it.

2. **Write an R function to calculate p value, and as v is large, we approximately have T~N(0,1) under H~0~.**

3. **Results**: p-values for all SNPs that are fitted simultaneously by a Ridge model

  > > **Computation only takes several minutes.**
  > >
  > > min(pval)
  > > [1] 0.0001160818
  > >
  > > **summary(pval)**
  > >
  > >
  > > | Min.      | 1st Qu.       | Median        | Mean      | 3rd Qu.       | Max       |
  > > | --------- | ------------- | ------------- | --------- | ------------- | --------- |
  > > | 0.0001161 | **0.2383228** | **0.4787446** | 0.4910669 | **0.7479664** | 0.9977548 |
  > >
  > > **A little bit smaller than uniform distribution?**

![image-20200710064017723](C:\Users\22367\Documents\markdown\image-20200710064017723.png)

**Comment:**

The results show consistency as it is not far from expected null hypothesis. This means that the statistical modelling procedure is not wrong! Also QQ plot shows that it has some power. The next step is to ensure the testing is powerful enough. And I should further investigate the p-value distribution given different lambda and power/FPR analysis using different p-value threshold in a simulation study.

##### Variable selection strategies in ridge

1. simple **thresholding on p value by different cutoffs.** This can be studied in simulations.
2. p-value+stability selection: in each of the bootstrap-like sub-sample, use p-value thresholding to select variables. Find the most stable subset.
3. **various pre-screening approaches**
   - simulate null phenotypes capturing LD information
   - obtain null p-value and $\hat \beta$ distribution for all SNPs. Due to LD, some of them are correlated (similar)
   - thresholding on the difference between p value from real data and p value from null data

##### Choice of $\lambda$ in ridge:

1. **Hoerl, Kennard & Baldwin** [1975 paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-372#ref-CR20)

<img src="C:\Users\22367\Documents\markdown\image-20200529112224903.png" alt="image-20200529112224903" style="zoom:50%;" />

<img src="C:\Users\22367\Documents\markdown\image-20200529112310805.png" alt="image-20200529112310805" style="zoom:50%;" />

m = p = the number of predictors

This is not plausible due to ill-conditioning (high LD)---we cannot even get OLS. 

**rationale of this approach:**

<img src="C:\Users\22367\Documents\markdown\image-20200529125337618.png" alt="image-20200529125337618" style="zoom:50%;" />

(this can be shown by differentiation)

 **2. Heritability-based extension of this approach:**

Under certain conditions (orthonormal X and standardized Y):

**For genotype matrix X, this means standardization plus no LD!**

$\lambda_{min}=p(1-\sum h^2_{SNP_j})/\sum h^2_{SNP_j}$, where $\sum h^2_{SNP_j} $ is the **unknown** variance explained by this set of J SNPs, which can be estimated by OLS **under good condition (n>>p, low LD).**

- If we further assume (genetic architecture, i.e. true model as) complete polygenicity so that all these p SNPs included follow a effect size distribution $N(0,\sigma^2_g)$, where $\sigma^2_g=h^2/m$ by CLT when m is large enough. (under complete polygenicity,**No. of  causal SNPs used=p, total No. of causal SNPs=total No. of SNPs=m**)

$\lambda_{min}=p(1-p \sigma^2_g)/(p \sigma^2_g)=1/\sigma^2_g-p=m/h^2-p$

- Or if we assume a point-normal effect size distribution as true model: 

  $\beta_j \sim (1-\pi)\delta+\pi N(0,\sigma^2_g)$, where $h^2=\sigma^2_gm\pi$

$\lambda_{min}=p(1-p\pi \sigma^2_g)/(p \pi\sigma^2_g)=1/(\pi\sigma^2_g)-p=m/h^2-p$

**The above two are the same because of assuming that the heritability are distributed evenly across genome.**

If we assume that the proportion of causal SNPs in p SNPs are not global $\pi$, but local $\pi'$

$\lambda_{min}=p(1-p\pi '\sigma^2_g)/(p \pi'\sigma^2_g)=1/(\pi'\sigma^2_g)-p=\frac{m\pi}{h^2\pi'}-p$

The ratio $\pi/\pi'$ is a tuning parameter. Particularly, suppose all the causal SNP are exactly these p SNPs $\pi'=1$, then, $\lambda=m\pi/h^2-p$

**Does heritability-based choice of $\lambda$ work better than ridge with fixed $\lambda$?**

If we use $\lambda=m/h^2-p\sim 10M$, seems too large. The algorithm cannot work?

Compare to $\lambda_{BLUP}=p*(1-h^2_{SNP})/h^2_{SNP}=p*(1-p/m*h^2)/(p/m*h^2)=m/h^2-p$

Exactly the same one! BLUP is the natural extension of HKB!



3.  **Cross-validation** (on a smaller sample size):

It could be computationally implausible even on a sample size of 20000 when p~ 10000.

4. **Bayesian perspective:**  Lawless & Wang [1976 paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-372#ref-CR21)

$\lambda=1/\sigma^2_0$ in Bayesian formulation, modelling the variance of normal prior $N(0,\sigma^2)$ again as unknown hyperparameters, $\lambda$ following a inverse-gamma prior (with mean as pre-determined ;good' one such as from approach one) , for example. Then update the posterior distribution of $\lambda$ based on data.

5. **Cule et al 2011, 2012** [paper2012](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.21750) [paper2011](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-372)in *Ridge* R Package that is based on PCA of X, similar to HKB method but in an adaptive way so that it works for p>n.

   

   > We found in order to obtain the best predictive performance a relaxed penalty was necessary. **CV to choose the parameters for the elastic net showed that that more than half of the time the alpha parameter was set to zero, which makes EN equivalent to RR** [Friedman et al., [2010](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.21750#gepi21750-bib-0014)]. The CV choice of the parameter that controls the amount of shrinkage of the regression coefficients resulted in a large parameter being chosen, implying strong shrinkage. **Together, the CV choice of EN parameters results in the equivalent of an RR with the regression coefficients shrunk close to zero. This reflects the structure of the data with many predictors of small effect.** RR‐CV and RR‐k offer approximately equivalent predictive performance, but RR-k is computationally more straightforward.

6. **AIC, BIC**

#####  Relation to heritability?

http://www.stat.yale.edu/~mjk56/temp/bigmemory-vignette.pdf

<img src="C:\Users\22367\Documents\markdown\image-20200529135435237.png" alt="image-20200529135435237" style="zoom:50%;" />

<img src="C:\Users\22367\Documents\markdown\image-20200529135455427.png" alt="image-20200529135455427" style="zoom:50%;" />

What is the relationship between c and $\lambda$?

   In the orthonormal case, **X**′**X**=**I** 

<img src="C:\Users\22367\Documents\markdown\image-20200529135646784.png" alt="image-20200529135646784" style="zoom:50%;" />

Therefore Lagrangian multiplier:

<img src="C:\Users\22367\Documents\markdown\image-20200529135748240.png" alt="image-20200529135748240" style="zoom:50%;" />

KKT condition for it: 

<img src="C:\Users\22367\Documents\markdown\image-20200529135947185.png" alt="image-20200529135947185" style="zoom:50%;" />

c = s = constraint = sum of ridge beta estimate

   <img src="C:\Users\22367\Documents\markdown\image-20200529115221272.png" alt="image-20200529115221272" style="zoom:50%;" />

   where the $\sum\beta_{OLS}^2=\hat h^2 \approx h^2$ when n is large enough.

shinkage: $1/(1+\lambda)$

#####    BLUP and Ridge: relationship

MTBLUP   LDpred

> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5841449/

![image-20200909234409570](C:\Users\22367\Documents\markdown\image-20200909234409570.png)

#### FWER control and FDR control in Lasso

###### Knockoffs to control FDR

R package gives threshold on $\lambda$ in Lasso

###### to control FWER

###### Real data application of testing

**Main challenge: algebra of large matrix**

1. Background setting

- R package *lmRidge* for fitting

- n=10k, p=1000 (from chr22) , 965 remained after MAF filtering 0.01% (exclude SNPs that are with too small MAF), phenotype=height

- lambda=0.1 (arbitrarily given)

- LD is troublesome for this package. Need to do some kind of LD clumping~~

  After that for all remaining SNPs LD<=**threshold=0.9**, and 884 SNPs pass it.

2. **Write an R function to calculate p value, and as v is large, we approximately have T~N(0,1) under H~0~.**

3. **Results**: p-values for all SNPs that are fitted simultaneously by a Ridge model

  > > **Computation only takes several minutes.**
  > >
  > > min(pval)
  > > [1] 0.0001160818
  > >
  > > **summary(pval)**
  > >
  > >
  > > | Min.      | 1st Qu.       | Median        | Mean      | 3rd Qu.       | Max       |
  > > | --------- | ------------- | ------------- | --------- | ------------- | --------- |
  > > | 0.0001161 | **0.2383228** | **0.4787446** | 0.4910669 | **0.7479664** | 0.9977548 |
  > >
  > > **A little bit smaller than uniform distribution?**

![image-20200710064017723](C:\Users\22367\Documents\markdown\image-20200710064017723.png)

**Comment:**

The results show consistency as it is not far from expected null hypothesis. This means that the statistical modelling procedure is not wrong! Also QQ plot shows that it has some power. The next step is to ensure the testing is powerful enough. And I should further investigate the p-value distribution given different lambda and power/FPR analysis using different p-value threshold in a simulation study.

#### Simulation comparison

##### Basic setting

###### Genotypes

sample size **n** = **10k**

SNP number **p** = **1000** , the first 1000 SNPs of chromosome 22 

using genotyped SNPs, some remaining holes imputed by median

- 964 SNPs remained after filtering those SNPs without any variation for these 10k people (exclude SNPs that are with too small MAF--too rare)

- **Strong LD could be troublesome for ridge computation. Need to do some kind of LD clumping**~~

  After that for all remaining SNPs LD<=**threshold=0.9**, and 883 SNPs pass it.

###### Generating phenotype: first only generate one continuous phenotype

a linear model for this very small genome

$Y_{n\times r}=XB_{p\times r}+E_{n\times r}$

trait number **r=1**

$E$ is a random residual error matrix, for each column k corresponding to trait k,

$\epsilon_i^{(k)} i.i.d.  N(0, \frac{Var(X\beta)(1-h^{(k)^2})}{h^{(k)2}}),i\in\N^n$

Note: For each k, $Var(X\beta)/Var(y)=1/(1+\frac{h^2}{(1-h^2)})=h^2$

###### Genetic architecture

causal SNP proportion $\pi_c\in$ {0.01,,0.025,0.05,0.1,0.125, 0.25,0.5,0.75,0.9,0.95}

causal SNP number $m=ceiling[883*\pi_c]$

heritability  $h^2\in ${0.1, 0.25, 0.5, 0.75, 0.9}

**In total 25 simulations for a method**

**In each simulation, select a $h^2, \pi_C$→generate phenotype→fit model, calculate FPR, FDR, Fscore and power. **

For a fixed heritability, changing $\pi_c$ 10 times will get 10 points for one method.

**Using $h^2/\pi_C$ and $m_{causal\ snps}$ as x-axis.**

##### Ridge

 **penalty** lambda=0.1 (arbitrarily given)

**Parameters P=(h^2^, $\pi$, number of causal SNPs, h^2^ per SNP)**

###### $\pi$=0.01

![image-20200710112437610](C:\Users\22367\Documents\markdown\image-20200710112437610.png)

###### $\pi$=0.025

![image-20200710114321901](C:\Users\22367\Documents\markdown\image-20200710114321901.png)





###### $\pi$=0.05

![image-20200710114619676](C:\Users\22367\Documents\markdown\image-20200710114619676.png)

**QQ plot for E and A**

![	](C:\Users\22367\Documents\markdown\image-20200710151131730.png)

###### $\pi$=0.1

![image-20200710114729805](C:\Users\22367\Documents\markdown\image-20200710114729805.png)

###### $\pi$=0.125

![image-20200710114145180](C:\Users\22367\Documents\markdown\image-20200710114145180.png)



###### $\pi$=0.25

![image-20200710144427607](C:\Users\22367\Documents\markdown\image-20200710144427607.png)

###### $\pi$=0.5

![image-20200710144828851](C:\Users\22367\Documents\markdown\image-20200710144828851.png)

###### $\pi$=0.75

![image-20200710145129167](C:\Users\22367\Documents\markdown\image-20200710145129167.png)



###### $\pi$=0.95

![image-20200710145429723](C:\Users\22367\Documents\markdown\image-20200710145429723.png)

**QQ plot for E and A**: ridge regression

![image-20200710150035269](C:\Users\22367\Documents\markdown\image-20200710150035269.png)

**QQ plot for E and A**: **per SNP GWAS for comparison-**--less powered

![image-20200710154456801](C:\Users\22367\Documents\markdown\image-20200710154456801.png)

**Tentative conclusion**

1. In general, higher heritability and higher per SNP heritability, better performance by ridge regression. Compared to GWAS that finds nothing under Bonferroni correction, ridge regression is much better.
2. In terms of FPR, p-value thresholding after ridge regression seems to control it very well even at a non-stringent cutoff.
3. When per SNP heritability > 0.01 (1%), ridge regression has considerable power (>50%) , except when heritability is as low as 0.1. And there seems to be certain "good" cutoff that can control FDR pretty well whilst retain adequate power, which lies between 10^-3^~10^-4^. (this should be related to N~SNP~?) This is close to Bonferroni correction.
4. When per SNP heritability goes below 0.005 (0.5%), true signals may be messed up with noise, so we see that power drops quite rapidly, suggesting that for ridge regression, there might be some obvious "detectability bound" on the signal intensity (beta effect size). However, we see in some cases that when h^2^ is as high as 0.9, SNPs with per SNP h^2^ as small as 0.001~0.002 is still detectable.
5. If we consider whole genome, 0.5% per SNP heritability for a h^2^ =60% trait corresponds to 120 SNPs, then true $\pi$ =120/10^6^ ~10^-4^ ~ 0.01%, very small! That suggests polygenic disease would be challenging, from the aspect of  whole-genome/chromosome wide variable selection. However, if we know that on a small region there are enriched heritability (causal SNPs), then the 1% per SNP heritability may be satisfied. A special advantage of ridge compared to per SNP GWAS is that simultaneous fitting takes LD into consideration.

##### improvement

**BH procedure**: in simulation, compare empirical BH FDR vs true FDR

empirical ones vs theoretical one

replications

recheck GWAS results

illustrate point 5 

large p, at least similar n/p ratio to real data's one

calculate a single index to make the comparison easier? For example, the number of false positives when the significance p-value cutoff is such that half the true associations are significant.

scale up the simulation to real data scale (whole chromosome, n>100k)

compare lambda_formula given vs lambda CV



#### Simulation comparison v2: RMT-provided $\lambda$

##### modification: RMT-based

- **similar n/p ratio** to real data's one

in UKBB we have n=500k,p~630k, maybe let's assume p/n→1

Here we will use the first **1000 SNPs** of chromosome 22 ,  so now we first use n~1k

- using genotyped SNPs
  -  some remaining holes imputed by median
  -  filtering those SNPs without any variation for these 10k people(exclude SNPs that are with too small MAF--too rare) ---> p=964
  -  Strong LD could be troublesome for ridge computation. Need to do some kind of LD clumping~~After that for all remaining SNPs LD<=threshold=0.9 ----> p=883
- take n=883, so p/n=1

> Bingxin Zhao and Hongtu Zhu. Cross-trait prediction accuracy of high-dimensional ridge-type estimators in genome-wide association studies, 2019.
>
> [Dobriban and Wager, 2018]. 

Applying the result from this paper, we will have

$\lambda_{opt}=p/n*(1-h_{SNP}^2)/h_{SNP}^2=(1-h_{SNP}^2)/h_{SNP}^2$

vs.

$\lambda_{BLUP}=M*(1-h_{SNP}^2)/h_{SNP}^2$

##### results

![](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.01.png)

![ridge_cutoff0814_0.9_0.1](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.1.png)

![ridge_cutoff0814_0.9_0.05](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.05.png)

![ridge_cutoff0814_0.9_0.5](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.5.png)

![ridge_cutoff0814_0.9_0.9](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.9.png)

![ridge_cutoff0814_0.9_0.025](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.025.png)

![ridge_cutoff0814_0.9_0.25](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.25.png)

![ridge_cutoff0814_0.9_0.75](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.75.png)

![ridge_cutoff0814_0.9_0.95](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.95.png)

![ridge_cutoff0814_0.9_0.125](C:\Users\22367\Documents\markdown\ridge_cutoff0814_0.9_0.125.png)

##### Verification for optimality of given lambda

###### $h^2=0.5, \pi=10\%$

optimal $\lambda$ =1, panel E

![ridge_cutoff0814lambda_0.5](C:\Users\22367\Documents\markdown\ridge_cutoff0814lambda_0.5.png)

![](C:\Users\22367\Documents\markdown\ridge_0814_differlambda_0.5.png)

###### $h^2=0.9, \pi=10\%$

optimal $\lambda$ = 0.1111, panel E

![ridge_cutoff0814lambda_0.9](C:\Users\22367\Documents\markdown\ridge_cutoff0814lambda_0.9.png)

![ridge_0814_differlambda_0.9](C:\Users\22367\Documents\markdown\ridge_0814_differlambda_0.9.png)

##### Comparison given CV $\lambda$

###### $h^2=0.75, \pi=10\%$

compare lambda: from 0.01, to formula given lambda=1/3, to CV lambda:=12, to BLUP lambda:=1000*1/3=333

**black line**: formula given,

 **red line**: CV



![image-20200814195730635](C:\Users\22367\Documents\markdown\image-20200814195730635.png)

If we look at panel B, cutoff=0.01 for p_value, giving highest power when controlling FWER<5% (panel A red line FWER>5% ) (edited) 

![image-20200814200818287](C:\Users\22367\Documents\markdown\image-20200814200818287.png)

Formula-given lambda (black line,0.25/0.75=1/3) are approximately the smallest lambda that gives similar performance as CV, given the best cutoff. 

BLUP lambda will be ~ 1000*1/3=333.3

Review: the predictors we used include all the causal SNPs. So $h^2_{SNP}=h^2$. If our data only account for a small **proportion L** of the real genome, then $h^2_{SNP}\approx L*h^2$



