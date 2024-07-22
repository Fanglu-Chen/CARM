#'@import arrangements
#'@import dplyr
#'@importFrom MASS ginv
#'@importFrom stats median
#'@importFrom stats na.omit
#'@importFrom stats var
#'@export
#'@title  Adaptive Randomization via Mahalanobis distance for Multi-arm design
#'
#'@description Randomize patients into treatment groups
#'for multi-arm trials using ARMM proposed by Haoyu Yang, Yichen Qin,
#'Yang Li, Fan Wang, and Feifang Hu.(2022)
#'
#'@param covariate a data frame. A row of the dataframe
#'corresponds to the covariate profile of a patient.
#'@param assignment a vector. If partial patients had been allocated
#', please input their allocation. IF
#'all the patients are not be allocated, please input
#''assignment = NA' directly.
#'@param K an integer; number of arms of the trial.
#'@param q the biased coin probability.
#'\eqn{q} should be larger than 1/2 and less than 1, default = 0.75
#'@param method Methods for calculating Mahalanobis distance, input one of these texts:
#''mean', 'max' or 'median'.
#'
#'
#'@details
#'Suppose \eqn{n} units (participants) are to be assigned to \eqn{K}
#'treatment groups. For each unit \eqn{i, i = 1, ..., n} and
#'treatment \eqn{j, j = 1, ..., K}, define the assignment
#'matrix \eqn{[T_{ij}]^{n*K}}, where
#'\eqn{T_{ij}=1} indicates unit \eqn{i} receives treatment \eqn{j}.
#' Consider \eqn{p} continuous covariates, let \eqn{x_i =
#' (x_{i1},...,x_{in})^T}.
#'
#'
#'
#'
#'The proposed method, namely the adaptive randomization
#'via Mahalanobis distance for multi-arm design (ARMM),
#'is outlined below. The implement of ARMM is similar to ARM.
#'
#'First assume that \eqn{n} units are in a sequence
#'and then assign the first \eqn{K} units to \eqn{K} treatment
#'groups randomly
#'as the initialization. Then,
#'the following units are assigned in blocks of \eqn{K}
#'sequentially and
#'adaptively until all the units
#'are assigned. For \eqn{K} units are assigned to \eqn{K}
#'groups, there are in total \eqn{K!} possible allocations.
#'Calculate \eqn{K!} potential overall
#'covariate imbalance measurement
#'according to pairwise Mahalanobis
#'distance under the \eqn{K!} possible allocations.
#'Choose the allocation which corresponds to the smallest
#'Mahalanobis
#'distance with a probability of \eqn{q} across all potential allocations.
#'Repeat the process until all units are assigned.
#'
#
#'
#'For any pair of treatments \eqn{s} and \eqn{t} among the \eqn{K}
#'treatment groups, calculate the Mahalanobis distance by:
#'
#'{\deqn{M_{s,t}(n) = 2n/K/K(\hat{x}_1 -\hat{x}_2)^Tcov(x)^{-1}(\hat{x}_1 -\hat{x}_2)}}
#'
#'In total, there are \eqn{C_K^2} pairs of Mahalanobis
#'distances among \eqn{K} treatment groups.Finally, calculate
#'the mean, the median or the maximum to represent the total imbalance.
#'
#'See the reference for more details.
#'
#'
#'

#'@return
#'An object of class "ARMM" is a list containing the following components:
#'\item{assignment}{Allocation of patients.}
#'\item{sample_size}{The number of patients from treatment 1 to treatment \eqn{K} respectively.}
#'\item{Mahalanobis_Distance}{Mahalanobis distance among treatment groups .}
#'
#'@references Yang H, Qin Y, Wang F, et al. Balancing covariates in multi-arm trials via adaptive randomization. Computational Statistics & Data Analysis, 2023, 179: 107642. https://doi.org/10.1016/j.csda.2022.107642


#'@examples
#'library(MASS)
#'#simulate covariates of patients
#'p <- 6; n <- 30
#'sigma <- diag(p); mean <- c(rep(0,p))
#'data <- mvrnorm(n, mean, sigma)
#'covariate <- as.data.frame(data)
#'#IF all the patients are not be allocated
#'ARMM(covariate = covariate, assignment = NA, K = 3, q = 0.75, method = 'mean')
#'#IF you had allocated partial patients
#'ARMM(covariate = covariate, assignment = c(1,2), K=4, q=0.75, method = 'max')




ARMM<-function(covariate, assignment, K, q=0.75, method){

  n<-nrow(covariate)
  p<-ncol(covariate)
  if (is.na(assignment)[1]){
    assignment<-data.frame(assignment=rep(NA,n))
  }
  else{
    aln<-length(assignment)
    assignment<-data.frame(assignment=c(assignment,rep(NA,n-aln)))
  }

  assigndata<-cbind(data.frame(assignment),
                    data.frame(covariate))
  names(assigndata)[1]<-'assignment'
  assigndata<-dplyr::arrange(assigndata, is.na(assigndata$assignment))


  if (nrow(assigndata[is.na(assigndata$assignment),])==n)
  {
    assigndata$assignment<-c(seq(1,K,1),rep(NA,(n-K)))
    assigndata<-circle_random(assigndata,K,p,q,method,n)
  }
  else
  {
    noassign<-setdiff(seq(1,K,1),unique(na.omit(assigndata$assignment)))
    if (length(noassign)==0){
      assigndata<-circle_random(assigndata,K,p,q,method,n)
    }
    else
    {
      assigndata$assignment<-c(assigndata[!is.na(assigndata$assignment),'assignment'],noassign,
                               rep(NA,(n-length(assigndata[!is.na(assigndata$assignment),'assignment'])-
                                         length(noassign))))
      assigndata<-circle_random(assigndata,K,p,q,method,n)
    }
  }

  R = NULL
  R$assignment<-assigndata[,1]
  R$sample_size<-as.data.frame(assigndata%>%group_by(assignment)%>%count(assignment))
  R$Mahalanobis_Distance<-pairwise_dis(assigndata,p,K,method)
  return(R)
}
