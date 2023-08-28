#' Sort a partition by first-appearance order
#' 
#' Sort a partition (cluster membership vector) by the first-appearance order
#' That is, let t_j:=min\{i: z_i =j\}, then j<k implies t_j<t_k.
#'
#' @param z *vector&lt;int&gt; (n)*, unsorted partition, represented as cluster membership vector
#'
#' @return sorted partition (cluster membership vector)
#' @export
#'
#' @examples
#'
#' # partition_sort(c(2,1,5,1,2))
#' # should be 1 2 3 2 1
#' 
partition_sort <- function(z){ 
  k = length(unique(z))
  temp = as.integer(factor(z, labels = 1:k))
  replace(temp, unique(temp), 1:k)[temp]
}


#' Check partition validity
#' 
#' If partition with k clusters is represented as cluster labels, it must have labels 1,2,...k. 
#' This does not require partition to be sorted in first-appearance order but must include labels 1 through k without missing labels between. 
#'  
#'
#' @param z *vector&lt;int&gt; (n)*, partition
#'
#' @return *logical*, FALSE if invalid 
#' @export
#'
#' @examples
#' partition_vaild(c(2,1,5,1,2)) # FALSE, since label 3 and 4 are missing
#' partition_vaild(c(2,1,3,1,2)) # TRUE
#' partition_vaild(partition_sort(c(2,1,3,1,2))) # always TRUE
#'
partition_vaild <- function(z){
  k = length(unique(z))
  all(sort(unique(z))==1:k)
}

