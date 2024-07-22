#'@import arrangements
#'@import dplyr
#'@importFrom MASS ginv
#'@importFrom stats median
#'@importFrom stats na.omit
#'@importFrom stats var
#'@export
#'@title Adaptive Randomization via Mahalanobis Distance
#'
#'@description Allocates patients to one of two treatments using
#'Adaptive Randomization via Mahalanobis Distance proposed by
#'Yichen Qin,Yang Li, Wei Ma, Haoyu Yang, and Feifang Hu.(2022)
#'@param covariate a data frame. A row of the dataframe
#'corresponds to the covariate profile of a patient.

#'@param assignment a vector. If partial patients had been allocated
#', please input their allocation. IF
#'all the patients are not be allocated, please input
#''assignment = NA' directly.
#'@param q the biased coin probability.
#'\eqn{q} should be larger than 1/2 and less than 1, default = 0.75

#'@details
#'Suppose that \eqn{n} patients are to be assigned to two treatment groups.
#'Consider \eqn{p}  continuous covariates for each patient.
#'\eqn{T_i} is the assignment of the \eqn{i}th patient.
#'The proposed procedure to assign units to treatment groups, namely adaptive
#'randomization via Mahalanobis distance (ARM), is outlined below.
#'
#'(1) Arrange all \eqn{n} units randomly into a sequence
#'\eqn{x_1,...,x_n}.
#'
#'(2) Assign the first two units with \eqn{T_1=1} and \eqn{T_2=2}.
#'
#'(3) Suppose that \eqn{2i} units have been assigned to
#'treatment groups,
#'for the \eqn{2i+1}-th and \eqn{2i+2}-th units:
#'
#' (3a) If the \eqn{2i+1}-th unit is assigned to treatment 1 and
#' the \eqn{2i+2}-th
#' unit to treatment 2, then calculate the potential
#' Mahalanobis distance, between the updated treatment groups.
#' with \eqn{2i+2} units, \eqn{M_1(2i + 2)}.
#'
#' (3b) Similarly, if the \eqn{2i+1}-th unit is
#' assigned to treatment 2 and
#' the \eqn{2i+2}-th unit to treatment 1, then calculate the
#' other potential Mahalanobis distance, \eqn{M_2(2i + 2)}.
#'
#'
#' (4) Assign the \eqn{2i+1}-th unit to treatment groups
#' according to the
#' following probabilities:
#'
#'if \eqn{ M_1(2i + 2) < M_2(2i + 2)}, \eqn{P(T_{2i+1} = 1)= q};
#'
#'if \eqn{ M_1(2i + 2) > M_2(2i + 2)}, \eqn{P(T_{2i+1} = 1)= 1-q};
#'
#'if \eqn{ M_1(2i + 2) = M_2(2i + 2)}, \eqn{P(T_{2i+1} = 1)= 0.5}.
#'
#'
#' (5) Repeat the last two steps until all units are assigned. If n is odd,
#' assign the last unit to two treatments with equal probabilities.
#'
#'Mahalanobis distance \eqn{M(n)} between the sample means across
#'different treatment groups is:
#'
#'\deqn{M(n)= np(1-p)(\hat{x_1} - \hat{x_2})^Tcov(x)^{-1}(\hat{x_1} - \hat{x_2}}
#'
#'
#'See the reference for more details.


#'
#'
#'@return
#'An object of class "ARM" is a list containing the following components:
#'\item{assignment}{Allocation of patients.}
#'\item{sample_size}{The number of patients in treatment 1 and treatment 2 respectively.}
#'\item{Mahalanobis_Distance}{Mahalanobis distance between treatment groups 1 and 2.}
#'
#'
#'@references Qin, Y., Y. Li, W. Ma, H. Yang, and F. Hu (2022). Adaptive randomization via mahalanobis distance. Statistica Sinica.DOI:<10.5705/ss.202020.0440>.

#'@examples
#'library(MASS)
#'#simulate covariates of patients
#'p <- 6; n <- 30
#'sigma <- diag(p); mean <- c(rep(0,p))
#'data <- mvrnorm(n, mean, sigma)
#'covariate <- as.data.frame(data)
#'#IF all the patients are not be allocated
#'ARM(covariate = covariate, assignment = NA, q=0.75)
#'#IF you had allocated partial patients
#'ARM(covariate = covariate,assignment = c(1,2),q=0.75)




ARM<-function(covariate, assignment, q=0.75){
  K=2
  method='none'

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


