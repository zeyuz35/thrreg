#' Threshold Estimation
#'
#' Computes estimates, standard errors, and confidence intervals of parameters
#' in the threshold regression model.
#'
#' @param df Data frame.
#' @param yi Integer or character; index or column name of dependent (y)
#'     variable in \code{df}.
#' @param xi Integer or character vector; indexes or column names of
#'     independent (x) variables in \code{df}.
#' @param qi Integer or character; index or column name of threshold (q)
#'     variable in \code{df}.
#' @param h Integer; heteroskedasticity indicator.
#'     Set \code{h = 0} to impose homoskedasticity assumption;
#'     set \code{h = 1} to use White-correction for heteroskedasticity.
#' @param test.pvalue Numeric; p-value of the threshold test returned by
#'     \code{\link{thr_test_hom}} or \code{\link{thr_test_hom}}.
#' @param var.names Character vector; variable names with
#'     \code{length(var.names) == ncol(df)} corresponding to columns in
#'     \code{df} to be used in threshold regression table.
#'     Default is \code{colnames(df)}.
#' @param conf2 Numeric; confidence level for first step of two-step
#'     confidence regions for regression parameters. Default is \code{conf2 = .8}.
#' @param nonpar Integer; indicator for non-parametric method used to estimate
#'     nuisance scale in the presence of heteroskedasticity
#'     (only relevant if \code{h = 1}).
#'     Set \code{nonpar = 1} to estimate regressions using a quadratic;
#'     set \code{nonpar= 2} (default) to estimate regressions using an
#'     Epanechnikov kernel with automatic bandwidth.
#' @param graph Logical; graph indicator.
#'     Set \code{TRUE} (default) to view the graph of the concentrated
#'     likelihood in gamma;
#'     set \code{FALSE} otherwise.
#' @param signif.level Character; indicator for notation of statistical
#'     significance levels.
#'     Set \code{signif.level = "stars"} (default) to use stars:
#'     \eqn{*p < 0.1}, \eqn{**p < 0.05}, \eqn{***p < 0.01}.
#'     Set \code{signif.level = "colors"} to use red and blue tones for
#'     statistically significant positive and negative estimates, respectively.
#'     (see \strong{Note} for LaTeX specification of red and blue tones.)
#' @param inf.crit Logical; if \code{TRUE}, information criteria (AIC, BIC, HQC)
#'     for each regime are shown in threshold regression table. Default is
#'     \code{FALSE}.
#' @param digits Integer; number of decimal places to be used for estimated
#'     coefficients in threshold regression table. Default is \code{digits = 3}.
#'     (Will be used in \code{format(round(x, digits = digits), nsmall = digits)},
#'     see \code{\link{format}}, \code{\link{round}}.)
#' @param integer.digits Integer; number of integer digits (i.e. digits before
#'     the decimal point) to be used for estimated coefficients in threshold
#'     regression table. If \code{NULL} (default), maximum value will be inquired.
#' @param digits.thr Integer; number of decimal places to be used for threshold
#'     estimate in threshold regression table. Default is
#'     \code{digits.thr = digits}. (See \code{digits} for usage.)
#' @param header Character; header to be used in threshold regression table.
#' @param output.short Logical; if \code{FALSE} (default), full output is printed.
#'     If \code{TRUE}, only threshold regression table is printed.
#' @param signif.legend Logical; if \code{TRUE} (default), legend of
#'     significance levels is printed below threshold regression table.
#'
#' @details Do not include a constant in the independent variables;
#'         the function automatically adds an intercept to the regression.
#'
#' @return Least squares estimate of threshold parameter. Output is printed to
#'     console unless redirected to file with \code{\link{sink}} (see examples).
#'
#' @note
#'     \enumerate{
#'     \item Load these packages in LaTeX file:
#'     \code{
#'     \\usepackage{booktabs} \\usepackage[table]{xcolor} \\usepackage{siunitx}}.
#'     \item If \code{signif.level = "colors"}, add these commands in LaTeX file
#'           to obtain red and blue tones:
#'     \code{
#'     \\newcommand{\\RedA}{red} \\newcommand{\\RedB}{red!60}
#'     \\newcommand{\\RedC}{red!30} \\newcommand{\\BlueA}{blue!75}
#'     \\newcommand{\\BlueB}{blue!55} \\newcommand{\\BlueC}{blue!30}}.
#'     }
#'
#' @inherit thrreg author references
#'
# #' @family threshold regression functions
#'
#' @seealso
#'     \code{\link{thr_test_hom}} and \code{\link{thr_test_het}} for threshold
#'     tests under homoskedasticity and heteroskedasticity, respectively.
#'
#' @keywords htest models ts
#'
# #' @example R/thr_est_examples.R
# #' @inherit thr_est_examples examples
# #' @inheritSection thr_est_examples Packageexamples
#'
#' @examples
#' \donttest{
#' ## Performs part of the empirical work reported in Hansen (2000)
#' data <- dur_john
#' test <- thr_test_het(data, 1, 2:5, 6)
#'
#' qhat <- thr_est(data, 1, 2:5, 6, 1, test$p_value)
#' qhat
#' }
#'
#' @importFrom base graphics legend lines plot title
#' @importFrom stats pchisq pnorm
#'
#' @export
#'

# Note that this newer implementation still vomits out latex code to the console
# but this is safe to ignore
# and now have proper return values you can access
thr_est <- function(
  df, yi, xi, qi, h,
  test.pvalue = 0.05,
  var.names = colnames(df),
  conf2 = .95,
  nonpar = 2,
  graph = FALSE,
  latex = FALSE,
  signif.level = "stars",
  inf.crit = FALSE,
  digits = 3, integer.digits = NULL, digits.thr = digits,
  header = NULL, output.short = TRUE, signif.legend = TRUE) {
  # Warning Messages
  if ((h != 0)*(h != 1)) {
    cat("You have entered h = ", h, "\n",
        "This number must be either 0 (homoskedastic case) or 1 (heteoskedastic)",
        "\n", "The function will either crash or produce invalid results", "\n")
  }
  if ((nonpar != 1)*(nonpar != 2)*(h==1)){
    cat("You have entered nonpar = ", nonpar, "\n",
        "This number should be either 1 (quadratic regression)", "\n",
        "or 2 (kernel regression)", "\n",
        "The function will employ the quadratic regression method", "\n", "\n")
  }
  # Prepare data into matrix format -----------------------------
  dat <- as.matrix(df)
  if (is.character(yi)) yi <- which(colnames(df) == yi)
  if (is.character(xi)) xi <- which(colnames(df) %in% xi)
  if (is.character(qi)) qi <- which(colnames(df) == qi)

  n <- nrow(dat)
  q <- dat[,qi]
  qs <- order(q)
  q <- q[qs]
  y <- as.matrix(dat[qs,yi])
  x <- cbind(matrix(c(1),n,1),dat[qs,xi])
  k <- ncol(x)
  yname <- var.names[yi]
  qname <- var.names[qi]
  xname <- rbind("Const",as.matrix(var.names[xi]))
  # Begin actual optimization ---------------------------------
  mi <- solve(t(x)%*%x, tol = 1e-1000)
  beta <- mi%*%(t(x)%*%y)
  e <- y-x%*%beta
  ee <- t(e)%*%e
  sig <- ee/(n-k)
  xe <- x*(e %*% matrix(c(1),1,k))
  # variance covariance matrix
  if (h == 0) {
    varcov <- (mi)*sig
    se <- varcov |> diag() |> sqrt()
  }
  else{
    varcov <- (mi %*% t(xe) %*% xe %*% mi)
    se <- varcov |> diag() |> sqrt()
    }
  vy <- sum((y - mean(y))^2)
  r_2 <- 1-ee/vy

  # MK: Compute additional statistics (adjusted R^2, AIC, BIC, HQC) for linear
  # model (global OLS estimation, without threshold)
  r_2_adj <- 1 - (1-r_2)*(n-1)/(n-(k-1)-1)
  aic <- n*log(ee/n) + 2*k
  bic <- n*log(ee/n) + k*log(n)
  hqc <- n*log(ee/n) + 2*k*log(log(n))

  qs <- unique(q)
  qn <- length(qs)
  sn <- matrix(c(0),qn,1)

  irb <- matrix(c(0),n,1)
  mm <- matrix(c(0),k,k)
  sume <- matrix(c(0),k,1)
  ci <- 0

  r <- 1
  while (r <= qn){
    irf <- (q <= qs[r])
    ir <- irf - irb
    irb <- irf
    ci <- ci + sum(ir)
    xir <- as.matrix(x[ir%*%matrix(c(1),1,k)>0])
    xir <- matrix(xir,nrow=nrow(xir)/k,ncol=k)
    mm <- mm + t(xir)%*%xir
    xeir <- as.matrix(xe[ir%*%matrix(c(1),1,k)>0])
    xeir <- matrix(xeir,nrow=nrow(xeir)/k,ncol=k)
    sume <- sume + colSums(xeir)
    mmi <- mm - mm%*%mi%*%mm
    if ((ci > k+1)*(ci < (n-k-1))) {
      sn[r] <- ee - t(sume)%*%solve(mmi, tol = 1e-1000)%*%sume
    }
    else {
      sn[r] <- ee
      }
    r <- r+1
  }

  rmin <- which.min(sn)
  smin <- sn[rmin]
  qhat <- qs[rmin]
  sighat <- smin/n

  i1 <- (q <= qhat)
  i2 <- (1 - i1) > 0
  x1 <- as.matrix(x[i1%*%matrix(c(1),1,k)>0])
  x1 <- matrix(x1,nrow=nrow(x1)/k,ncol=k)
  y1 <- as.matrix(y[i1])
  x2 <- as.matrix(x[i2%*%matrix(c(1),1,k)>0])
  x2 <- matrix(x2,nrow=nrow(x2)/k,ncol=k)
  y2 <- as.matrix(y[i2])
  mi1 <- solve(t(x1)%*%x1, tol = 1e-1000)
  mi2 <- solve(t(x2)%*%x2, tol = 1e-1000)
  beta1 <- mi1%*%(t(x1)%*%y1)
  beta2 <- mi2%*%(t(x2)%*%y2)
  e1 <- y1 - x1%*%beta1
  e2 <- y2 - x2%*%beta2
  ej <- rbind(e1,e2)
  n1 <- nrow(y1)
  n2 <- nrow(y2)
  ee1 <- t(e1) %*% e1
  ee2 <- t(e2) %*% e2
  sig1 <- ee1/(n1-k)
  sig2 <- ee2/(n2-k)
  sig_jt <- (ee1 + ee2)/(n - k*2)

  # Variance covariance matrices
  if (h == 0) {
    # variance covariance matrix
    varcov1 <- mi1 * sig_jt
    varcov2 <- mi2 * sig_jt
    se1 <- varcov1 |> diag() |> sqrt()
    se2 <- varcov2 |> diag() |> sqrt()
  } else {
    xe1 <- x1*(e1 %*% matrix(c(1), 1, k))
    xe2 <- x2*(e2 %*% matrix(c(1), 1, k))
    varcov1 <- (mi1 %*% t(xe1) %*% xe1 %*% mi1)
    varcov2 <- (mi2 %*% t(xe2) %*% xe2 %*% mi2)
    se1 <- varcov1 |> diag() |> sqrt()
    se2 <- varcov2 |> diag() |> sqrt()
  }
  vy1 <- sum((y1 - mean(y1))^2)
  vy2 <- sum((y2 - mean(y2))^2)
  r2_1 <- 1 - ee1/vy1
  r2_2 <- 1 - ee2/vy2
  r2_joint <- 1 - (ee1+ee2)/vy

  # MK: Compute additional statistics (adjusted R^2, AIC, BIC, HQC) for full
  # threshold model
  # factor 2*(k-1)+1 in adjusted R^2: number of regressors (excluding the
  # constant) + threshold estimate
  r2_adj_joint <- 1 - (1-r2_joint)*(n-1)/(n-(2*(k-1)+1)-1)
  aic_joint <- n*log((ee1+ee2)/n) + 2*(2*k+1)
  bic_joint <- n*log((ee1+ee2)/n) + (2*k+1)*log(n)
  hqc_joint <- n*log((ee1+ee2)/n) + 2*(2*k+1)*log(log(n))

  # MK: Compute additional statistics (adjusted R^2, AIC, BIC, HQC) for each
  # regime
  r2_adj_1 <- 1 - (1-r2_1)*(n1-1)/(n1-(k-1)-1)
  r2_adj_2 <- 1 - (1-r2_2)*(n2-1)/(n2-(k-1)-1)
  aic_1 <- n1*log(ee1/n1) + 2*k
  aic_2 <- n2*log(ee2/n2) + 2*k
  bic_1 <- n1*log(ee1/n1) + k*log(n1)
  bic_2 <- n2*log(ee2/n2) + k*log(n2)
  hqc_1 <- n1*log(ee1/n1) + 2*k*log(log(n1))
  hqc_2 <- n2*log(ee2/n2) + 2*k*log(log(n2))


  if (h == 0) lr <- (sn-smin)/sighat
  if (h == 1){
    r1 <- (x %*% (beta1-beta2))^2
    r2 <- r1*(ej^2)
    qx <- cbind(q^0,q^1,q^2)
    qh <- cbind(qhat^0,qhat^1,qhat^2)
    m1 <- qr.solve(qx,r1)
    m2 <- qr.solve(qx,r2)
    g1 <- qh%*%m1
    g2 <- qh%*%m2
    if (nonpar==2){
      sigq <- sqrt(mean((q-mean(q))^2))
      hband <- 2.344*sigq/(n^(.2))
      u <- (qhat-q)/hband
      u2 <- u^2
      f <- mean((1-u2)*(u2<=1))*(.75/hband)
      df <- -mean(-u*(u2<=1))*(1.5/(hband^2))
      eps <- r1 - qx%*%m1
      sige <- (t(eps)%*%eps)/(n-3)
      hband <- sige/(4*f*((m1[3]+(m1[2]+2*m1[3]*qhat)*df/f)^2))
      u2 <- ((qhat-q)/hband)^2
      kh <- ((1-u2)*.75/hband)*(u2<=1)
      g1 <- mean(kh*r1)
      g2 <- mean(kh*r2)
    }
    eta2 <- g2/g1
    lr <- (sn-smin)/eta2
  }

  # MK
  tstars1 <- rep("", k)
  tstars2 <- rep("", k)
  # Intiialize empty lists
  qhat1_list <- list()
  qhat2_list <- list()
  beta1l_list <- list()
  beta1u_list <- list()
  beta2l_list <- list()
  beta2u_list <- list()

  # MK: Loop over different confidence levels ---------------------------
  # (I think this is really unnecessary)
  conf1_levels <- c(.9, .95, .99)

  for (i in 1:length(conf1_levels)) {

    conf1 <- conf1_levels[i]
    c1 <- -2*log(1 - sqrt(conf1))
    c2 <- -2*log(1 - sqrt(conf2))
    lr1 <- (lr >= c1)
    lr2 <- (lr >= c2)
    if (max(lr1) == 1) {
      qhat1 <- qs[which.min(lr1)]
      qhat2 <- qs[qn + 1 - which.min(rev(lr1))]
    }
    # qhat se?
    else {
      qhat1 <- qs[1]
      qhat2 <- qs[qn]
    }
    z <- which.max((pnorm(seq(.01, 3, by = .01))*2 - 1) >= conf1)/100;
    beta1l <- beta1 - se1*z
    beta1u <- beta1 + se1*z
    beta2l <- beta2 - se2*z
    beta2u <- beta2 + se2*z
    r <- 1
    while (r<=qn){
      if (lr2[r]==0){
        i1 <- (q <= qs[r])
        x1 <- as.matrix(x[i1%*%matrix(c(1),1,k)>0])
        x1 <- matrix(x1,nrow=nrow(x1)/k,ncol=k)
        y1 <- y[i1]
        if (qr(t(x1)%*%x1)$rank==ncol(t(x1)%*%x1)){
          mi1 <- solve(t(x1)%*%x1, tol = 1e-1000)
          b1 <- mi1%*%(t(x1)%*%y1)
          e1 <- y1 - x1%*%b1
          if (h==0){
            ser1 <- as.matrix(sqrt(diag(mi1)*(t(e1)%*%e1)/(nrow(y1)-k)))
          }else{
            xe1 <- x1*(e1%*%matrix(c(1),1,k))
            ser1 <- as.matrix(sqrt(diag(mi1%*%t(xe1)%*%xe1%*%mi1)))
          }
          beta1l <- apply((rbind(t(beta1l),t(b1 - ser1*z))),2,min)
          beta1u <- apply((rbind(t(beta1u),t(b1 + ser1*z))),2,max)
        }
        i2 <- (1-i1)>0
        x2 <- as.matrix(x[i2%*%matrix(c(1),1,k)>0])
        x2 <- matrix(x2,nrow=nrow(x2)/k,ncol=k)
        y2 <- y[i2]
        if (qr(t(x2)%*%x2)$rank==ncol(t(x2)%*%x2)){
          mi2 <- solve(t(x2)%*%x2, tol = 1e-1000)
          b2 <- mi2%*%(t(x2)%*%y2)
          e2 <- y2 - x2%*%b2
          if (h==0){
            ser2 <- as.matrix(sqrt(diag(mi2)*(t(e2)%*%e2)/(nrow(y2)-k)))
          }else{
            xe2 <- x2*(e2%*%matrix(c(1),1,k))
            ser2 <- as.matrix(sqrt(diag(mi2%*%t(xe2)%*%xe2%*%mi2)))
          }
          beta2l <- apply((rbind(t(beta2l),t(b2 - ser2*z))),2,min)
          beta2u <- apply((rbind(t(beta2u),t(b2 + ser2*z))),2,max)
        }
      }
      r <- r+1
    }

    # Define notation for significance levels for coefficient estimates
    if (signif.level == "stars") {
      # Define number of stars
      nstars <- (if (conf1 == 0.9) "$^\\ast$"
                 else if (conf1 == 0.95) "$^{\\ast\\ast}$"
                 else if (conf1 == 0.99) "$^{\\ast\\ast\\ast}$")

      tstars1[(0 <= beta1l | 0 >= beta1u)] <- nstars
      tstars2[(0 <= beta2l | 0 >= beta2u)] <- nstars

    } else if (signif.level == "colors") {
      # Define red/blue tones for positive/negative estimates
      colpos <- (if (conf1 == 0.9) "\\cellcolor{\\RedC}"
                 else if (conf1 == 0.95) "\\cellcolor{\\RedB}"
                 else if (conf1 == 0.99) "\\cellcolor{\\RedA}")
      colneg <- (if (conf1 == 0.9) "\\cellcolor{\\BlueC}"
                 else if (conf1 == 0.95) "\\cellcolor{\\BlueB}"
                 else if (conf1 == 0.99) "\\cellcolor{\\BlueA}")

      tstars1[((0 <= beta1l | 0 >= beta1u) & beta1 > 0)] <- colpos
      tstars1[((0 <= beta1l | 0 >= beta1u) & beta1 < 0)] <- colneg
      tstars2[((0 <= beta2l | 0 >= beta2u) & beta2 > 0)] <- colpos
      tstars2[((0 <= beta2l | 0 >= beta2u) & beta2 < 0)] <- colneg
    }

    qhat1_list[[i]] <- qhat1
    qhat2_list[[i]] <- qhat2
    beta1l_list[[i]] <- beta1l
    beta1u_list[[i]] <- beta1u
    beta2l_list[[i]] <- beta2l
    beta2u_list[[i]] <- beta2u

  }

  # Heteroskedasticity test
  het_test <- function(e, x) {
    e2 <- e^2
    x2 <- x^2
    v <- e2 - x2%*%qr.solve(x2,e2)
    e2 <- e2 - colMeans(e2)
    te <- nrow(e)%*%(1-(t(v)%*%v)/(t(e2)%*%e2))
    out <- 1-pchisq(te,ncol(x))
    out
  }

  # Latex output -----------------------------------------------------------
  if (latex == TRUE) {
    # Long Output -----------------------------------------------------------
    if (output.short == FALSE) {
      # Title ----------------------------------------------------------
      cat(paste0("\\subsection{Estimate Sample Split, Using ", qname, "}"), "\n\n")

      # Null Model (No Threshold) ---------------------------------------
      cat("\\subsubsection{Global OLS Estimation, Without Threshold}", "\n")
      cat("Dependent Variable:     ", yname, "\\\\\n")
      if (h==1) cat("Heteroskedasticity Correction Used", "\\\\\\\\")
      if (h==0) cat("OLS Standard Errors Reported", "\\\\\\\\")
      cat("\n")
      cat("\\begin{tabular}{l*{2}{r}}", "\\toprule", sep = "\n")
      cat("Variable ", "  &  ", "Estimate  ", "  &  ", "St Error", "\\\\\n")
      cat("\\midrule", "\n")
      tbeta <- format(beta, digits=4)
      tse <- format(se, digits=4)
      for (j in 1:k) {cat(xname[j], "  &  ", tbeta[j], "  &  ", tse[j], "\\\\\n")}
      cat("\\bottomrule", "\\end{tabular}", sep = "\n")
      cat("\\bigskip \\\\")
      cat("\n")
      cat("Observations:                      ", n, "\\\\\n")
      cat("Degrees of Freedom:                ", (n-k), "\\\\\n")
      cat("Sum of Squared Errors:             ", ee, "\\\\\n")
      cat("Residual Variance:                 ", sig, "\\\\\n")
      cat("R-squared:                         ", r_2, "\\\\\n")
      # MK
      cat("Adjusted R-squared:                ", r_2_adj, "\\\\\n")
      cat("AIC:                               ", aic, "\\\\\n")
      cat("BIC:                               ", bic, "\\\\\n")
      cat("HQC:                               ", hqc, "\\\\\n")
      # HTSK Test
      cat("Heteroskedasticity Test (P-Value): ", het_test(e, x), "\\\\\n")
      cat("\n")

      # Joint Model Output ---------------------------------------------
      # This is just a summary of joint model selection criteria,
      # not actual need for outputting joint model intercepts
      cat("\\subsubsection{Threshold Estimation}", "\n")
      cat("Threshold Variable:                ", qname, "\\\\\n")
      cat("Threshold Estimate:                ", qhat, "\\\\\n")
      tqhat1 <- format(qhat1_list, digits=4)
      tqhat2 <- format(qhat2_list, digits=4)
      for (i in 1:length(conf1_levels)) {
        tit <- paste(c("["),tqhat1[[i]],", ",tqhat2[[i]],c("]"),sep="")
        cat(conf1_levels[i], "Confidence Interval:          ", tit, "\\\\\n")
      }
      cat("Sum of Squared Errors:             ", (ee1+ee2), "\\\\\n")
      cat("Residual Variance:                 ", sig_jt, "\\\\\n")
      cat("Joint R-squared:                   ", r2_joint, "\\\\\n")
      # MK
      cat("Joint Adjusted R-squared:          ", r2_adj_joint, "\\\\\n")
      cat("Joint AIC:                         ", aic_joint, "\\\\\n")
      cat("Joint BIC:                         ", bic_joint, "\\\\\n")
      cat("Joint HQC:                         ", hqc_joint, "\\\\\n")
      # HTSK Test
      cat("Heteroskedasticity Test (P-Value): ", het_test(ej,x), "\\\\\n")
      cat("\n")

      # Regime 1 Output -----------------------------------------------
      tit <- paste(qname,"$\\leq$",format(qhat,digits=6),sep="")
      cat("\\subsubsection*{Regime 1:", tit, "}", "\n")
      cat("Parameter Estimates", "\\\\\n")
      cat("\\begin{tabular}{l*{2}{r}}", "\\toprule", sep = "\n")
      cat("Variable ", "  &  ", "Estimate  ", "  &  ", "St Error", "\\\\\n")
      cat("\\midrule", "\n")
      tbeta1 <- format(beta1, digits=4)
      tse1 <- format(se1, digits=4)
      for (j in 1:k) {cat(xname[j], "  &  ", tbeta1[j], "  &  ", tse1[j], "\\\\\n")}
      cat("\\bottomrule", "\\end{tabular}", sep = "\n")
      cat("\\bigskip \n")
      cat("\n")

      for (i in 1:length(conf1_levels)) {
        cat(conf1_levels[i], "Confidence Regions for Parameters", "\\\\\n")
        cat("\\begin{tabular}{l*{2}{r}}", "\\toprule", sep = "\n")
        cat("Variable ", " &   ", "Low         ", "  &  ", "High", "\\\\\n")
        cat("\\midrule", "\n")
        tbeta1l <- format(beta1l_list[[i]], digits=4)
        tbeta1u <- format(beta1u_list[[i]], digits=4)
        for (j in 1:k) {cat(xname[j], "  &  ", tbeta1l[j], "  &  ", tbeta1u[j],
                            "\\\\\n")}
        cat("\\bottomrule", "\\end{tabular}", sep = "\n")
        cat("\\bigskip \n")
        cat("\n")
      }

      cat("Observations:                      ", n1, "\\\\\n")
      cat("Degrees of Freedom:                ", (n1-k), "\\\\\n")
      cat("Sum of Squared Errors:             ", ee1, "\\\\\n")
      cat("Residual Variance:                 ", sig1, "\\\\\n")
      cat("R-squared:                         ", r2_1, "\\\\\n")
      # MK
      cat("Adjusted R-squared:                ", r2_adj_1, "\\\\\n")
      cat("AIC:                               ", aic_1, "\\\\\n")
      cat("BIC:                               ", bic_1, "\\\\\n")
      cat("HQC:                               ", hqc_1, "\\\\\n")
      #
      cat("\n")

      # Regime 2 Output --------------------------------------------------
      tit <- paste(qname,"$>$",format(qhat,digits=6),sep="")
      cat("\\subsubsection*{Regime 2:", tit, "}", "\n")
      cat("Parameter Estimates", "\\\\\n")
      cat("\\begin{tabular}{l*{2}{r}}", "\\toprule", sep = "\n")
      cat("Variable ", "  &  ", "Estimate  ", "  &  ", "St Error", "\\\\\n")
      cat("\\midrule", "\n")
      tbeta2 <- format(beta2, digits=4)
      tse2 <- format(se2, digits=4)
      for (j in 1:k) {cat(xname[j], "  &  ", tbeta2[j], "  &  ", tse2[j], "\\\\\n")}
      cat("\\bottomrule", "\\end{tabular}", sep = "\n")
      cat("\\bigskip \n")
      cat("\n")

      for (i in 1:length(conf1_levels)) {
        cat(conf1_levels[i], "Confidence Regions for Parameters", "\\\\\n")
        cat("\\begin{tabular}{l*{2}{r}}", "\\toprule", sep = "\n")
        cat("Variable ", "  &  ", "Low         ", "  &  ", "High", "\\\\\n")
        cat("\\midrule", "\n")
        tbeta2l <- format(beta2l_list[[i]], digits=4)
        tbeta2u <- format(beta2u_list[[i]], digits=4)
        for (j in 1:k) {cat(xname[j], "  &  ", tbeta2l[j], "  &  ", tbeta2u[j],
                            "\\\\\n")}
        cat("\\bottomrule", "\\end{tabular}", sep = "\n")
        cat("\\bigskip \n")
        cat("\n")
      }

      cat("Observations:                      ", n2, "\\\\\n")
      cat("Degrees of Freedom:                ", (n2-k), "\\\\\n")
      cat("Sum of Squared Errors:             ", ee2, "\\\\\n")
      cat("Residual Variance:                 ", sig2, "\\\\\n")
      cat("R-squared:                         ", r2_2, "\\\\\n")
      # MK
      cat("Adjusted R-squared:                ", r2_adj_2, "\\\\\n")
      cat("AIC:                               ", aic_2, "\\\\\n")
      cat("BIC:                               ", bic_2, "\\\\\n")
      cat("HQC:                               ", hqc_2, "\\\\\n")
      #
      cat ("\n")

    }
  }

  # Produce a graph of the criterion function and related test
  if (graph == TRUE) {
    xxlim <- range(qs)
    yylim <- range(rbind(lr,c1))
    clr <- matrix(c(1),qn,1)*c1
    plot(
      qs, lr, lty = 1, col = 1,
      xlim = xxlim, ylim = yylim, type = "l", ann = 0)
    lines(qs, clr, lty = 2, col = 2)
    xxlab <- paste(c("Threshold Variable: "),qname, sep = "")
    title(
      main = "Confidence Interval Construction for Threshold",
      xlab = xxlab, ylab = "Likelihood Ratio Sequence in gamma")
    tit <- paste(conf1*100,c("% Critical"), sep = "")
    legend(
      "bottomright", c("LRn(gamma)", tit),
      lty = c(1, 2), col = c(1, 2))
  }

  # Define notation for significance levels for threshold estimate
  test_pval <- test.pvalue

  if (signif.level == "stars") {
    # Define number of stars
    tstars_thr <- (if (test_pval <= 0.01) "$^{\\ast\\ast\\ast}$"
                   else if (test_pval <= 0.05) "$^{\\ast\\ast}$"
                   else if (test_pval <= 0.1) "$^\\ast$")

  } else if (signif.level == "colors") {
    # Define red/blue tones for positive/negative estimate
    tstars_thr <- if (qhat > 0) {
      (if (test_pval <= 0.01) "\\cellcolor{\\RedA}"
       else if (test_pval <= 0.05) "\\cellcolor{\\RedB}"
       else if (test_pval <= 0.1) "\\cellcolor{\\RedC}")
    } else if (qhat < 0) {
      (if (test_pval <= 0.01) "\\cellcolor{\\BlueA}"
       else if (test_pval <= 0.05) "\\cellcolor{\\BlueB}"
       else if (test_pval <= 0.1) "\\cellcolor{\\BlueC}")
    }
  }

  if (is.null(integer.digits)) {
    integer.digits <- max(nchar(n1), nchar(n2))
  }

  # More latex output ------------------------------------------------
  if (latex == TRUE) {
    # MK: Summarizing threshold regression table with Regime 1 (left) and
    # Regime 2 (right)
    if (output.short == FALSE) {
      cat("\\subsubsection*{Threshold Regression Table}", "\n")
    }
    cat(paste0("\\begin{tabular}{l*{2}{S[table-format=", integer.digits, ".",
               digits,"]r}}"),
        "\\toprule", sep = "\n")
    if (!is.null(header)) {
      cat("&", paste0("\\multicolumn{4}{c}{", header, "}"), "\\\\\n")
      cat("\\cmidrule(lr){2-5}", "\n")
    }
    cat("&", "\\multicolumn{2}{c}{Regime 1}", "&", "\\multicolumn{2}{c}{Regime 2}",
        "\\\\\n")
    cat("\\cmidrule(lr){2-3}", "\\cmidrule(lr){4-5}", sep = "\n")
    cat("&", paste0("\\multicolumn{2}{c}{", qname, tstars_thr, " $\\leq$ ",
                    format(round(qhat, digits=digits.thr), nsmall=digits.thr), "}"),
        "&", paste0("\\multicolumn{2}{c}{", qname, tstars_thr, " $>$ ",
                    format(round(qhat, digits=digits.thr), nsmall=digits.thr), "}"),
        "\\\\\n")
    cat("\\cmidrule(lr){2-3}", "\\cmidrule(lr){4-5}", sep = "\n")
    cat("{Variable}", "&", "{Estimate}", "&", "{Std error}", "&", "{Estimate}", "&",
        "{Std error}", "\\\\\n")
    cat("\\midrule", "\n")
    tbeta1 <- format(round(beta1, digits=digits), nsmall=digits)
    tse1 <- format(round(se1, digits=digits), nsmall=digits)
    tbeta2 <- format(round(beta2, digits=digits), nsmall=digits)
    tse2 <- format(round(se2, digits=digits), nsmall=digits)
    for (j in 1:k) {
      cat(xname[j], "&", paste0(tbeta1[j], tstars1[j]), "&", paste0("(", tse1[j], ")"),
          "&", paste0(tbeta2[j], tstars2[j]), "&", paste0("(", tse2[j], ")"), "\\\\\n")
    }
    r2_1 <- format(round(r2_1, digits=digits), nsmall=digits)
    r2_2 <- format(round(r2_2, digits=digits), nsmall=digits)
    r2_adj_1 <- format(round(r2_adj_1, digits=digits), nsmall=digits)
    r2_adj_2 <- format(round(r2_adj_2, digits=digits), nsmall=digits)
    aic_1 <- round(aic_1)
    aic_2 <- round(aic_2)
    bic_1 <- round(bic_1)
    bic_2 <- round(bic_2)
    hqc_1 <- round(hqc_1)
    hqc_2 <- round(hqc_2)
    cat("\\midrule", "\n")
    cat("\\#Obs", "&", n1, "& &", n2, "&", "\\\\\n")
    # cat("Degrees of Freedom ", "  &  ", (n1-k), "  &  & ", (n2-k), "  &  ", "\\\\\n")
    # cat("Sum of Squared Errors ", "  &  ", ee1, "  &  & ", ee2, "  &  ", "\\\\\n")
    # cat("Residual Variance ", "  &  ", sig1, "  &  & ", sig2, "  &  ", "\\\\\n")
    # cat("$R^2$ ", "  &  ", r2_1, "  &  & ", r2_2, "  &  ", "\\\\\n")
    cat("$R^2_\\text{adj}$", "&", r2_adj_1, "& &", r2_adj_2, "&", "\\\\\n")
    if (inf.crit == TRUE) {
      cat("AIC", "&", aic_1, "& &", aic_2, "&", "\\\\\n")
      cat("BIC", "&", bic_1, "& &", bic_2, "&", "\\\\\n")
      cat("HQC", "&", hqc_1, "& &", hqc_2, "&", "\\\\\n")
    }
    cat("\\bottomrule", "\\end{tabular}", sep = "\n")

    if (signif.legend == TRUE) {
      cat("\\smallskip \\\\\n")
      if (signif.level == "stars") {
        cat("$^\\ast p<0.1$; $^{\\ast\\ast} p<0.05$; $^{\\ast\\ast\\ast} p<0.01$", "\n")
      } else if (signif.level == "colors") {
        cat("\\colorbox{\\RedC}{\\makebox(20,6){}}/\\colorbox{\\BlueC}{\\makebox(20,6){}} $p<0.1$;",
            "\\colorbox{\\RedB}{\\makebox(20,6){}}/\\colorbox{\\BlueB}{\\makebox(20,6){}} $p<0.05$;",
            "\\colorbox{\\RedA}{\\makebox(20,6){}}/\\colorbox{\\BlueA}{\\makebox(20,6){}} $p<0.01$",
            sep = "\n")
      }
    }
  }
  # Clean up names
  xnames <- beta |> rownames()
  xnames[1] <- "(Intercept)"

  rownames(beta1) <- xnames
  rownames(beta2) <- xnames

  rownames(varcov) <- xnames
  colnames(varcov) <- xnames
  rownames(varcov1) <- xnames
  colnames(varcov1) <- xnames
  rownames(varcov2) <- xnames
  colnames(varcov2) <- xnames

  names(se) <- xnames
  names(se1) <- xnames
  names(se2) <- xnames

  # No Threshold --------------------------------------------
  null_model <- list(
    # Coefficients
    beta = beta,
    # Standard Errors
    se = se,
    varcov = varcov,
    # Other model selection components
    r2 = r_2,
    r2_adj = r_2_adj,
    aic = aic,
    bic = bic,
    hqc = hqc,
    e = e,
    # Sum of squared errors
    sse = ee,
    sig = sig,
    # HTSK test
    htsk_test = het_test(e, x)
  )
  # Full Model ----------------------------------------------
  full_model <- list(
    # Normal Summary stats
    n = n,
    dof = n - k,
    # Model Selection
    r2 = r2_joint,
    r2_adj = r2_adj_joint,
    aic = aic_joint,
    bic = bic_joint,
    hqc = hqc_joint,
    e = ej,
    sse = ee1 + ee2,
    htsk_test = het_test(ej, x)
  )
  # Regime 1 - this should be z < threshold -------------------
  regime_1 <- list(
    beta = beta1,
    # standard errors
    varcov = varcov1,
    se = se1,
    # Summary Stats
    n = n1,
    dof = n1 - k,
    # Other model selection components
    r2 = r2_1,
    r2_adj = r2_adj_1,
    aic = aic_1,
    bic = bic_1,
    hqc = hqc_1,
    e = e1,
    sse = ee1,
    sig = sig1
  )

  # Regime 2 - this should be z > threshold ------------------------
  regime_2 <- list(
    beta = beta2,
    # Standard errors,
    se = se2,
    varcov = varcov2,
    n = n2,
    dof = n2 - k,
    # Model Selection Criteria
    r2 = r2_2,
    r2_adj = r2_adj_2,
    aic = aic_2,
    bic = bic_2,
    hqc = hqc_2,
    e = e2,
    sse = ee2,
    sig = sig2
  )

  # Return --------------------------------------------------
  ret_list <- list(
    # The actual estimated threshold
    qhat = qhat,
    qhat_ci = c(qhat1_list[[2]], qhat2_list[[2]]),
    # global null model,
    null_model = null_model,
    # joint model
    full_model = full_model,
    # Individual regimes
    regime_1 = regime_1,
    regime_2 = regime_2
  )
  # Set class in case we want some methods for this later
  class(ret_list) <- "thrreg"
  return(ret_list)
}

# print method for thrreg ----------------------------------------------

print.thrreg <- function(obj, ...) {
  # This is a placeholder...
  # Print threshold first

  # Print Null Model

  # Print Joint Model

  # Print Regime 1

  # Print Regime 2
}

# coef() method -------------------------------------------------------

# This is useful for directly extracting the coefficients, just to make it
# a little easier to interface with

coef.thrreg <- function(
    obj, model = "all"
) {
  if (model == "all") {
    coef_list <- list(
      null_model = obj$null_model$beta,
      regime_1 = obj$regime_1$beta,
      regime_2 = obj$regime_2$beta
    )
  }
  else if (model == "null") {
    beta_list <- obj$null_model$beta
  }
  else if (model == "1") {
    beta_list <- obj$regime_1$beta
  }
  else if (model == "2") {
    beta_list <- obj$regime_2$beta
  }
  # return
  return(coef_list)
}

# vcov method ---------------------------------------------------------

# Similar to above, dedicated vcov method to make it easier to get standard
# errors, etc

# By default, this returns the vanilla variance covariance matrix of the
# coefficients

vcov.thrreg <- function(
  obj, type = "varcov", model = "all"
  ) {
  if (model == "all") {
    var_list <- list(
      null_model = obj$null_model$varcov,
      regime_1 = obj$regime_1$varcov,
      regime_2 = obj$regime_2$varcov
    )
    se_list <- list(
      null_model = obj$null_model$se,
      regime_1 = obj$regime_1$se,
      regime_2 = obj$regime_2$se
    )
  }
  else if (model == "null") {
    var_list <- obj$null_model$varcov
    se_list <- obj$null_model$se
  }
  else if (model == "1") {
    var_list <- obj$regime_1$varcov
    se_list <- obj$regime_1$se
  }
  else if (model == "2") {
    var_list <- obj$regime_2$varcov
    se_list <- obj$regime_2$se
  }
  # return
  if (type == "varcov") {
    return(var_list)
  }
  else if (type == "se") {
    return(se_list)
  }
}

# summary() --------------------------------------------------------
#
# Manually using coef() and varcov is still tedious
# so build in a summary function so that it is easier
# summary() will rely on coef and vcov

# model = "all" to show all models,
# can also be "null", "1", "2"
summary.thrreg <- function(obj, model = "all") {
  # Prepare coef list
  coef_list <- coef(obj)
  # Prepare se list
  se_list <- vcov(obj, type = "se")
  # Staple the two together
  ret_list <- foreach(i = 1:length(coef_list)) %do% {
    data.frame(
      coef = coef_list[[i]],
      se = se_list[[i]])
  }
  names(ret_list) <- c("null", "regime_1", "regime_2")
  if (model == "null") {
    ret_list <- ret_list$null
  }
  else if (model == "1") {
    ret_list <- ret_list$regime_1
  }
  else if (model == "2") {
    ret_list <- ret_list$regime_2
  }
  # return
  cat("Threshold Estimate:", obj$qhat)
  cat("\n")
  cat("Threshold Estimate CIs (95%):", obj$qhat_ci)
  cat("\n")
  print(ret_list)
}
