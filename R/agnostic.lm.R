#' Agnostic Linear Models
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights weights (that is, minimizing sum(w*e^2)); otherwise ordinary least squares is used. See also 'Details',
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param method the method to be used; for fitting, currently only method = "qr" is supported; method = "model.frame" returns the model frame (the same as with model = TRUE, see below).
#' @param model,x,y,qr logicals. If TRUE the corresponding components of the fit (the model frame, the model matrix, the response, the QR decomposition) are returned.
#' @param singular.ok logical. If FALSE (the default in S but not in R) a singular fit is an error.
#' @param contrasts an optional list. See the contrasts.arg of model.matrix.default.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one are specified their sum is used. See model.offset.
#' @param ... additional arguments to be passed to the low level regression fitting functions (see below).
#'
#' @return lm returns an object of class "agnostic.lm" or for multiple responses of class c("agnostic.mlm", "agnostic.lm"). There are print and summary methods availiable for these.
#' @export agnostic.lm
#' @S3method print,agnostic.lm
#'
#' @examples
#' mod1 <- agnostic.lm(rnorm(100) ~ rexp(100))
#' summary(mod1)
agnostic.lm <- function(formula, data, subset, weights, na.action, method = "qr",
                        model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
                        contrasts = NULL, offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
                     method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0,
                                                      3) else numeric(), residuals = y, fitted.values = 0 *
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=
                                                                                0) else if (is.matrix(y)) nrow(y) else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w))
      lm.fit(x, y, offset = offset, singular.ok = singular.ok,
             ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                 ...)
  }
  class(z) <- c(if (is.matrix(y)) "agnostic.mlm", "agnostic.lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr)
    z$qr <- NULL
  z
}

print.agnostic.lm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}
