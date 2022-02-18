# This script is independent of other code part of this package.
# It is kept here as it may be useful for some extra analysis or future integration.

#' Decompose distribution into components
#' @param dat a distribution represented as a numeric vector to decompose.
#' @param dist 'pois' or 'norm' (default).
#' @param seed seed number.
#' @param min_comp minimal number of components to fit, default is 1.
#' @param max_comp maximal number of components to fit, default is 10.
#' @param min_prior the minimum relative size of components, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' @param threshold default is `0.5`. Sometimes, the result components
#' include adjacent distributions with similar mu
#' (two and more distribution are very close), we use this threshold
#' to obtain a more meaningful fit with less components.
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' @param niter the maximum number of iterations.
#' @importClassesFrom flexmix FLXcontrol
#' @importFrom flexmix flexmix logLik
#' @export
#' @examples
#' set.seed(2021)
#' x <- c(rnorm(10, 0), rnorm(50, 1), rnorm(20, 4), rnorm(5, 10))
#' y <- gcap.extractComponents(x, max_comp = 5)
#' y
#' @testexamples
#' expect_equal(nrow(y$params), 3)
gcap.extractComponents <-
  function(dat,
           dist = "norm",
           seed = 2021L,
           model_selection = "BIC",
           threshold = 0.4,
           min_prior = 0.001,
           niter = 1000,
           nrep = 5,
           min_comp = 1,
           max_comp = 10) {
    control <- new("FLXcontrol")
    control@minprior <- min_prior
    control@iter.max <- niter

    stopifnot(max_comp >= min_comp)
    set.seed(seed, kind = "L'Ecuyer-CMRG")

    if (dist == "norm") {
      fit <-
        stepFlexmix_v2(
          dat ~ 1,
          model = flexmix::FLXMCnorm1(),
          k = seq(min_comp, max_comp),
          nrep = nrep,
          control = control
        )
      if (inherits(fit, "stepFlexmix")) {
        fit <- recur_fit_component(
          fit = fit, dist = dist,
          threshold = threshold,
          control = control,
          model_selection = model_selection
        )
      }
    } else if (dist == "pois") {
      fit <-
        stepFlexmix_v2(
          dat ~ 1,
          model = flexmix::FLXMCmvpois(),
          k = seq(min_comp, max_comp),
          nrep = nrep,
          control = control
        )
      if (inherits(fit, "stepFlexmix")) {
        fit <- recur_fit_component(
          fit = fit, dist = dist,
          threshold = threshold,
          control = control,
          model_selection = model_selection
        )
      }
    }

    list(
      fit = fit,
      params = get_component_parameter(fit)
    )
  }


recur_fit_component <- function(fit, dist, threshold, control, model_selection = "BIC") {
  fits <- fit
  fit <- flexmix::getModel(fits, which = model_selection)
  message("Select ", fit@k, " according to ", model_selection, ".")

  mu <- find_mu(fit)
  sub_len <- sum(diff(mu) < threshold)

  if (sub_len > 0) {
    message("The model does not pass the threshold for mu difference of adjacent distribution.")
    message("Trying to call the optimal model...")
  }

  while (sub_len > 0) {
    K <- fit@k
    fit <- tryCatch(flexmix::getModel(fits, which = as.character(K - sub_len)),
      error = function(e) {
        # If the model is not called by user
        # Call it
        if (dist == "norm") {
          flexmix::flexmix(
            dat ~ 1,
            model = flexmix::FLXMCnorm1(),
            k = K,
            control = control
          )
        } else {
          flexmix::flexmix(
            dat ~ 1,
            model = flexmix::FLXMCmvpois(),
            k = K,
            control = control
          )
        }
      }
    )

    mu <- find_mu(fit)
    sub_len <- sum(diff(mu) < threshold)
  }

  message("Finally, select ", fit@k, " after passing threshold ", threshold, ".")
  fit
}

find_mu <- function(fit) {
  mu <- flexmix::parameters(fit)
  if (is.matrix(mu)) {
    mu <- mu[1, ]
  }
  mu <- sort(mu, decreasing = FALSE)
  mu
}

get_component_parameter <- function(x) {
  paras <- flexmix::parameters(x)
  # weight is how many events assigned
  # to a cluster (component)
  # i.e. number of observations
  #
  # Info from package author:
  # the cluster sizes indicate the number
  # of observations assigned to each of the
  # clusters according to the a-posteriori probabilities.
  .get_weight <- function(mean, x) {
    wt_tb <- flexmix::clusters(x) %>%
      table()
    wt <- as.numeric(wt_tb)
    if (length(wt) == length(mean)) {
      return(wt)
    } else {
      names(wt) <- names(wt_tb)
      all_names <- seq_along(mean) %>%
        as.character()
      wt[setdiff(all_names, names(wt))] <- 0
      wt[sort(names(wt))] %>% as.numeric()
    }
  }

  if (is.null(dim(paras))) {
    # Assume it is pois distribution
    z <- data.table::data.table(
      mean = as.numeric(paras)
    )
    z$sd <- sqrt(z$mean)
    z$n <- .get_weight(z$mean, x)
  } else {
    # Assume it is normal distribution
    z <- data.table::data.table(
      mean = as.numeric(paras[1, ]),
      sd = as.numeric(paras[2, ])
    )
    z$n <- .get_weight(z$mean, x)
  }
  z
}


# stepFlex v2
stepFlexmix_v2 <- function(..., k = NULL, nrep = 3, verbose = TRUE, drop = TRUE, unique = FALSE) {
  MYCALL <- match.call()
  MYCALL1 <- MYCALL
  bestFlexmix <- function(...) {
    z <- new("flexmix", logLik = -Inf)
    logLiks <- rep(NA, length.out = nrep)
    for (m in seq_len(nrep)) {
      if (verbose) {
        cat(" *")
      }
      x <- try(flexmix(...), silent = TRUE)
      if (!is(x, "try-error")) {
        logLiks[m] <- logLik(x)
        if (logLik(x) > logLik(z)) {
          z <- x
        }
      }
    }
    return(list(z = z, logLiks = logLiks))
  }
  z <- list()
  if (is.null(k)) {
    RET <- bestFlexmix(...)
    z[[1]] <- RET$z
    logLiks <- as.matrix(RET$logLiks)
    z[[1]]@call <- MYCALL
    z[[1]]@control@nrep <- nrep
    names(z) <- as.character(z[[1]]@k)
    if (verbose) {
      cat("\n")
    }
  } else {
    k <- as.integer(k)
    logLiks <- matrix(nrow = length(k), ncol = nrep)
    for (n in seq_along(k)) {
      ns <- as.character(k[n])
      if (verbose) {
        cat(k[n], ":")
      }
      RET <- bestFlexmix(..., k = k[n])
      z[[ns]] <- RET$z
      logLiks[n, ] <- RET$logLiks
      MYCALL1[["k"]] <- as.numeric(k[n])
      z[[ns]]@call <- MYCALL1
      z[[ns]]@control@nrep <- nrep
      if (verbose) {
        cat("\n")
      }
    }
  }
  logLiks <- logLiks[is.finite(sapply(z, logLik)), , drop = FALSE]
  z <- z[is.finite(sapply(z, logLik))]
  rownames(logLiks) <- names(z)
  if (!length(z)) {
    stop("no convergence to a suitable mixture")
  }
  if (drop & (length(z) == 1)) {
    return(z[[1]])
  } else {
    z <- return(new("stepFlexmix",
      models = z, k = as.integer(names(z)),
      nrep = as.integer(nrep), logLiks = logLiks, call = MYCALL
    ))
    if (unique) {
      z <- unique(z)
    }
    return(z)
  }
}
