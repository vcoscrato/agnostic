#' Agnostic Student's t-Test
#'
#' Performs t-tests under the agnostic perspective.
#'
#' @param x A vector of data.
#' @param y A vector of data or a constant.
#' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param alpha Desired type I probability of error.
#' @param beta Desired tupe II probability of error.
#' @param paired Boolean indicating if x and y are paired data vectors.
#'
#' @return A list containing the decision, the limits of the critical region and the test statistic value for the test.
#' @export
#'
#' @examples #One sample test
#' agnostic.t.test(rnorm(10), 0)
#'
#' #Two sample paired test
#' agnostic.t.test(rnorm(20, mean = 1), rnorm(20), paired = TRUE)
#'
#' #Simulation example (this should take some seconds to run)
#' library(tidyr)
#' library(ggplot2)
#'
#' n=10
#' mu.grid=seq(-2,2,length.out = 50)
#' B=5000
#' decision=matrix(NA,length(mu.grid),B)
#' for(i in 1:length(mu.grid))
#' {
#'   #print(i/length(mu.grid))
#'   for(b in 1:B)
#'   {
#'     decision[i,b]=agnostic.t.test(x=rnorm(n,mean = mu.grid[i]), alternative = "two.sided")$decision
#'   }
#' }
#'
#' decisions.names=c("Agnostic","Reject H0","Accept H0")
#' probability.decisions=t(apply(decision,1,function(x)
#' {
#'   return(apply(as.matrix(decisions.names),1,function(y)
#'   {
#'     mean(x==y)
#'   }))
#' }))
#' colnames(probability.decisions)=decisions.names
#' probability.decisions=cbind(mu.grid,as.data.frame(probability.decisions))
#' probability.decisions=gather(probability.decisions,"Decision","Probability",Agnostic:`Accept H0`)
#'
#' theme = theme_set(theme_minimal(base_size = 26))
#' theme = theme_update(legend.position="top", legend.title=element_blank(),panel.grid.major.x=element_blank())
#'
#' ggplot()+geom_line(data=probability.decisions,aes(x=mu.grid,y=Probability,linetype=Decision,color=Decision),size=2)+
#' geom_hline(yintercept = 0.05)+xlab(expression(mu))+geom_vline(xintercept = 0)+ylab("Probability of each decision")
agnostic.t.test <- function(x, y = 0, alternative = "two.sided", alpha = 0.05, beta = 0.05, paired = FALSE) {
  if(paired == TRUE) {
    if(length(x) != length(y)) {
      stop("Paired test need 2 vector of same length")
    } else {
      x <- x - y
      y = 0
    }
  }
  if(length(y) > 1) {
    t = (length(x)+length(y)-2)^(1/2)*((mean(x) - mean(y))/(((1/length(x)) + (1/length(y)))*(sum((x-mean(x))^2)+sum((y-mean(y))^2)))^(1/2))
    if(alternative == "two.sided") {
      c0 <- qt(alpha/2, df = length(x)+length(y)-2, lower.tail = FALSE)
      c1 <- qt((1-beta)/2, df = length(x)+length(y)-2, lower.tail = FALSE)
      if(abs(t) >= c0)
        decision <- "Reject H0"
      else if(abs(t) <= c1)
        decision <- "Accept H0"
      else
        decision <- "Agnostic"
      return(list(decision = decision, limits = c(c1,c0), abs.test.statistic = abs(t)))
    } else if(alternative == "less") {
      c0 <- qt(alpha, df = length(x)+length(y)-2, lower.tail = FALSE)
      c1 <- qt(1-beta, df = length(x)+length(y)-2, lower.tail = FALSE)
      if(t >= c0)
        decision <- "Reject H0"
      else if(t <= c1)
        decision <- "Accept H0"
      else
        decision <- "Agnostic"
      return(list(decision = decision, limits = c(c1,c0), test.statistic = t))
    } else if(alternative == "greater") {
      c0 <- qt(alpha, df = length(x)+length(y)-2, lower.tail = FALSE)
      c1 <- qt(1-beta, df = length(x)+length(y)-2, lower.tail = FALSE)
      if(t >= c0)
        decision <- "Accept H0"
      else if(t <= c1)
        decision <- "Reject H0"
      else
        decision <- "Agnostic"
      return(list(decision = decision, limits = c(c1,c0), test.statistic = t))
    } else
      stop("invalid alternative parameter, use one of 'two.sided','less' or 'greater'")
  } else {
    t = ((mean(x) - y)/(sd(x)/sqrt(length(x))))
    if(alternative == "two.sided") {
      c0 <- qt(alpha/2, df = length(x)-1, lower.tail = FALSE)
      c1 <- qt((1-beta)/2, df = length(x)-1, lower.tail = FALSE)
      if(abs(t) >= c0)
        decision <- "Reject H0"
      else if(abs(t) <= c1)
        decision <- "Accept H0"
      else
        decision <- "Agnostic"
      return(list(decision = decision, limits = c(c1,c0), abs.test.statistic = abs(t)))
    } else if(alternative == "less") {
      c0 <- qt(alpha, df = length(x)-1, lower.tail = FALSE)
      c1 <- qt(1-beta, df = length(x)-1, lower.tail = FALSE)
      if(t >= c0)
        decision <- "Reject H0"
      else if(t <= c1)
        decision <- "Accept H0"
      else
        decision <- "Agnostic"
      return(list(decision = decision, limits = c(c1,c0), test.statistic = t))
    } else if(alternative == "greater") {
      c0 <- qt(alpha, df = length(x)-1, lower.tail = FALSE)
      c1 <- qt(1-beta, df = length(x)-1, lower.tail = FALSE)
      if(t >= c0)
        decision <- "Accept H0"
      else if(t <= c1)
        decision <- "Reject H0"
      else
        decision <- "Agnostic"
      return(list(decision = decision, limits = c(c1,c0), test.statistic = t))
    } else
      stop("invalid alternative parameter, use one of 'two.sided','less' or 'greater'")
  }
}
