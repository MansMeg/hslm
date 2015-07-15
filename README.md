## Bayesian linear regression using the horseshoe prior

As a part of my research I've become interested in the bayesian horseshoe prior for linear regression. But I had a hard time to find the nitty-gritty details of the derivations of the MCMC sampler or a simple implementation of the full MCMC sampler. This repo contains both the nitty-gritty details in deriving the sampler and a naive (slow) implementation in R.

### Reference
Carvalho, Carlos M., Nicholas G. Polson, and James G. Scott. "Handling sparsity via the horseshoe." International Conference on Artificial Intelligence and Statistics. 2009.
It can be found [here](http://jmlr.org/proceedings/papers/v5/carvalho09a/carvalho09a.pdf).

### The nitty-gritty details of the sampler 
The derivation of the sampler can be found [here](https://github.com/MansMeg/hslm/blob/master/Derivations/hslm.pdf). 

### Implementation in R

A naive (slow) implementation in R as an R package has been put together in this repository. The main purposes is to have a reference implementation and code that can be shared and developed further for others interested in this sampler. A faster implementation can be found in the ```monomvn``` R package.  

#### Install the package

```
install.packages("devtools")
devtools::install_github("MansMeg/hslm", subdir = "RPackage")
```

#### Basic usage

Training and test data (based on the diabetes dataset in the ```lars``` R package) is included to test the sampler. 

```
data(diabetes_x_train)
data(diabetes_y_train)
hs_res <- hslm(diabetes_y_train, diabetes_x_train)
colMeans(hs_res$beta[-(1:1000]),] # Mean of beta parameters with burnin (1000) removed
```
