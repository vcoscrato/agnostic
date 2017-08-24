#' Summarizing Agnostic Linear Model Fits
#'
#' @param object an object of class "agnostic.lm", usually, a result of a call to agnostic.lm.
#' @param alpha probability of type I error.
#' @param beta probability of type II error.
#' @param d a vector indicating the desired Cohen's effect size where the probability of type II error is beta for each parameter. 0 by default.
#' @param plot.power logical ; if TRUE, draws the power function for the tests.
#' @param correlation logical; if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor logical; if TRUE, print the correlations in a symbolic form (see symnum) rather than as numbers.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list countaining various metrics from a agnostic linear model fit, there is a print method availiable for this.
#' @export
#'
#' @examples
#' mod1 <- agnostic.lm(rnorm(100) ~ rexp(100))
#' summary(mod1)
summary.agnostic.lm <- function(object, alpha = 0.05, beta = 0.05, d = NULL,
                               plot.power = FALSE,
                               correlation = FALSE, symbolic.cor = FALSE, ...)
{
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate",
                                               "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms))
    stop("invalid 'lm' object:  no 'terms' component")
  if (is.null(r <- object$qr))
    stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  Qr <- r
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual)
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept"))
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) *
      1e-30)
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  est <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  pvalue <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  decision <- rep("Agnostic", length(pvalue))
  decision[pvalue < alpha] <- "Accept H0"
  decision[pvalue > 1-beta] <- "Reject H0"
  if(!is.null(d)) {
    find.c.0=function(d,df,element.inverse,beta){
      c.0=seq(0,10,0.005)
      beta.true=rep(NA,length(c.0))
      for(i in 1:length(c.0))
      {
        beta.true[i]=suppressWarnings(pt(c.0[i],df,(1/sqrt(element.inverse))*d))-
          suppressWarnings(pt(-c.0[i],df,(1/sqrt(element.inverse))*d))
      }
      return(c.0[FNN::get.knnx(beta.true,beta,k=1)$nn.index])
    }
    elements.inverse=diag(R)
    if(length(d)==1)
      d=rep(d,length(elements.inverse))
    c.0=rep(NA,length(d))
    for(i in 1:length(d))
      c.0[i]=find.c.0(d[i],object$df.residual,elements.inverse[i],beta)
    c.1=qt(1-alpha/2,object$df.residual)
    if(any(c.0>c.1))
      stop("c.0>c.1: either decrease beta, decrease d, or decrease alpha")
    decision <- rep("Accept H0", length(tval))
    decision[tval > c.0] <- "Agnostic"
    decision[tval > c.1] <- "Reject H0"
    if(plot.power)
    {
      theme = theme_set(theme_minimal(base_size = 20))
      theme = theme_update(legend.position="top", legend.title=element_blank(), panel.grid.major.x=element_blank(),
                           panel.grid.major.y=element_blank(),
                           panel.grid.minor.y=element_blank(),
                           panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                           axis.text.y = element_text(colour="black",size=13), axis.text.x = element_text(size=13),axis.ticks.y= element_line(colour="black"))+ theme_update(axis.ticks.length=unit(.15, "cm"),panel.spacing.y = unit(1.5, "lines"))
      d.grid=seq(0,1,length.out = 1000)
      prob.accept=prob.reject=prob.agnostic=matrix(NA,length(d.grid),length(elements.inverse))
      colnames(prob.accept)=names(elements.inverse)
      colnames(prob.reject)=names(elements.inverse)
      colnames(prob.agnostic)=names(elements.inverse)
      for(i in 1:length(elements.inverse))
      {
        for(j in 1:length(d.grid))
        {
          prob.accept[j,i]=suppressWarnings(pt(c.0[i],fit$df.residual,
                                               (1/sqrt(elements.inverse[i]))*d.grid[j]))-
            suppressWarnings(pt(-c.0[i],fit$df.residual,
                                (1/sqrt(elements.inverse[i]))*d.grid[j]))
          prob.reject[j,i]=suppressWarnings(pt(c.1,fit$df.residual,
                                               (1/sqrt(elements.inverse[i]))*d.grid[j],lower.tail = FALSE))+
            suppressWarnings(pt(-c.1,fit$df.residual,
                                (1/sqrt(elements.inverse[i]))*d.grid[j]))
          prob.agnostic[j,i]=1-prob.reject[j,i]-prob.accept[j,i]
        }
      }
      prob.accept=cbind("Accept",d.grid,gather(as.data.frame(prob.accept),coefficient,prob,1:ncol(prob.accept)))
      prob.reject=cbind("Reject",d.grid,gather(as.data.frame(prob.reject),coefficient,prob,1:ncol(prob.reject)))
      prob.agnostic=cbind("Agnostic",d.grid,gather(as.data.frame(prob.agnostic),coefficient,prob,1:ncol(prob.agnostic)))
      names(prob.accept)[1]=names(prob.reject)[1]=names(prob.agnostic)[1]="Decision"
      probs=rbind(prob.accept,prob.reject,prob.agnostic)
      g=ggplot(data=probs ,aes(x=d.grid, y=prob, group=Decision)) +
        geom_line(aes(linetype=Decision, color=Decision),size=1.5)+
        facet_wrap( ~ coefficient, ncol=round(sqrt(length(elements.inverse))))+
        ylab("Probability of the Decision")+xlab("Cohen's d effect size")+ 
        geom_hline(yintercept=beta,color="black")+ 
        annotate("text", min(prob.accept$d.grid), beta+0.05,size=7, label = "beta",parse=TRUE,color="black")+
        geom_hline(yintercept=alpha,color="black")+ 
        annotate("text", min(prob.accept$d.grid), alpha+0.05, size=7,label = "alpha",parse=TRUE,color="black")
      print(g)
    }
  }
  ans$residuals <- r
  ans$coefficients <- cbind(est, se, tval, pvalue, 1-pvalue, decision)
  dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]],
                                     c("Estimate", "Std. Error", "t value", "P.H0", "P.H1", "Decision"))
  ans$aliased <- is.na(coef(object))
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept"))
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
                                                       df.int)/rdf)
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
                        numdf = p - df.int, dendf = rdf)
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action))
    ans$na.action <- z$na.action
  class(ans) <- "summary.agnostic.lm"
  ans
}

print.summary.agnostic.lm <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
                                      signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat(if (!is.null(x$weights) && diff(range(x$weights)))
    "Weighted ", "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L)
      structure(apply(t(resid), 1L, quantile), dimnames = list(nam,
                                                               dimnames(resid)[[2L]]))
    else {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  }
  else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  }
  else {
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L])
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                                                              colnames(coefs)))
      coefs[!aliased, ] <- as.data.frame(x$coefficients)
    }
    coefs[,-ncol(coefs)] <- round(as.numeric(coefs[,-ncol(coefs)]), 2)
    print(coefs)
  }
  cat("\nResidual standard error:", format(signif(x$sigma,
                                                  digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared,
                                           digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L],
                                                                                       digits = digits), "on", x$fstatistic[2L], "and",
        x$fstatistic[3L], "DF,  p.H0:", format.pval(pf(x$fstatistic[1L],
                                                       x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE),
                                                    digits = digits), ", p.H1 =", format.pval(1-pf(x$fstatistic[1L],
                                                                                                   x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE),
                                                                                              digits = digits))
    cat("\n")
  }
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}
