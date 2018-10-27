# require(pracma, quietly=TRUE)
require(lbfgs, quietly=TRUE)
require(lbfgsb3, quietly=TRUE)

auglag_lbfgs <- function(par, fn, gr=NULL,
    hin, hin.jac=NULL, heq, heq.jac=NULL,
    bfgs_method=c("lbfgs", "lbfgsb3"),
    control.outer=list(), control.lbfgs=list()) {

    # if (is.null() && is.null()) stop("Problem is unconstrained")

    control.outer.default <- 
        list(lam0 = 10, sig0 = 100, eps = 1e-07, itmax = 50, 
             method = "BFGS", trace = FALSE, NMinit = FALSE,
             ilack.max = 6, i.scale = 1, e.scale = 1, kkt2.check = TRUE)
    control.outer <- modifyList(control.outer.default, control.outer)

    control.lbfgs.default <- list()

    if (is.null(hin.jac))
        hin.jac <- function(x) pracma::jacobian(hin, x)

    if (is.null(heq.jac))
        heq.jac <- function(x) pracma::jacobian(heq, x)

    ans <- auglag3(par, fn, gr,
        hin, hin.jac, heq, heq.jac,
        bfgs_method=bfgs_method,
        control.outer, control.lbfgs)

    return(ans)
}


auglag3 <- function(par, fn, gr,
    hin, hin.jac, heq, heq.jac,
    bfgs_method=bfgs_method,
    control.outer=list(), control.lbfgs=list()) {

    sig <- control.outer$sig0
    lam0 <- control.outer$lam0
    trace <- control.outer$trace
    eps <- control.outer$eps
    itmax <- control.outer$itmax
    ilack.max <- control.outer$ilack.max
    method <- control.outer$method
    NMinit <- control.outer$NMinit
    kkt2.check <- control.outer$kkt2.check
    pfact <- 1

    control.optim <- list(reltol= 1e-06)

    fun <- function(par) {
        h0 <- hin(par)
        i0 <- heq(par)
        d0 <- c(h0, i0)
        active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
        inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
        d0[active] <- h0[active]
        d0[inactive] <- lam[inactive]/sig
        fn(par) - pfact * sum(lam * d0) + pfact * sig/2 * sum(d0 * d0)
    }
    
    gradient <- function(par) {
        h0 <- hin(par)
        i0 <- heq(par)
        d0 <- c(h0, i0)
        active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
        inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
        ij <- rbind(hin.jac(par)[active, , drop = FALSE], heq.jac(par
            ))
        if (length(inactive) > 0) 
            gr(par) - pfact * colSums(lam[-inactive] * ij) + pfact * sig * 
                drop(t(ij) %*% d0[-inactive]) else gr(par) - pfact * colSums(lam * ij) + pfact * sig * drop(t(ij) %*% 
            d0)
    }

    h0 <- hin(par)
    i0 <- heq(par)
    d0 <- c(h0, i0)
    
    lam <- rep(lam0, length(d0))
    active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
    inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
    d0[active] <- h0[active]
    d0[inactive] <- lam[inactive]/sig
    dmax <- max(abs(d0))
    
    obj <- fn(par)
    r <- obj
    feval <- 0
    geval <- 0
    ilack <- 0
    Kprev <- dmax
    sig0 <- sig/Kprev
    if (is.infinite(sig0)) 
        sig0 <- 1
    sig <- sig0

    K <- Inf
    if (trace) 
        cat("Min(hin): ", min(h0), "Max(abs(heq)): ", max(abs(i0)), "\n")
    for (i in 1:itmax) {
        if (trace) {
            cat("Outer iteration: ", i, "\n")
            cat("Min(hin): ", min(h0), "Max(abs(heq)): ", max(abs(i0)), "\n")
            cat("par: ", signif(par, 6), "\n")
            cat("fval =  ", signif(obj, 4), "\n \n")
        }
        par.old <- par
        obj.old <- obj
        r.old <- r
        if (NMinit && (i == 1)) {
            a <- optim(par = par, fn = fun, control = control.optim,
                       method = "Nelder-Mead")
        } else {
            if (sig > 1e+05) 
                control.optim$reltol <- 1e-10
            if (bfgs_method=="lbfgs") {
                a <- lbfgs::lbfgs(fun, gradient, par, invisible=1)
                r <- a$value
            } else if (bfgs_method=="lbfgsb3") {
                a <- lbfgsb3::lbfgsb3(par, fun, gradient,
                                      lower=rep(-10,52), upper=rep(10,52))
                r <- a$prm
            } else {
                stop("Unknown method ...")
            }
            # feval <- feval + a$counts[1]
            # geval <- geval + a$counts[2]
        }
        par <- a$par
        h0 <- hin(par)
        i0 <- heq(par)
        d0 <- c(h0, i0)
        active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
        inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
        if (length(active) > 0) 
            d0[active] <- h0[active]
        if (length(inactive) > 0) 
            d0[inactive] <- lam[inactive]/sig
        K <- max(abs(d0))
        if (K <= Kprev/4) {
            lam <- lam - d0 * sig
            Kprev <- K
        } else sig <- 10 * sig
        obj <- fn(par)
        
        pconv <- max(abs(par - par.old))
        if (pconv < eps) {
            ilack <- ilack + 1
        } else ilack <- 0
        
        if ((is.finite(r) && is.finite(r.old) && abs(r - r.old) < eps && K < 
            eps) | ilack >= ilack.max) 
            break
    }

    if (i == itmax) {
        a$convergence <- 7
        a$message <- "ALABaMA ran out of iterations and did not converge"
    } else if (K > eps) {
        a$convergence <- 9
        a$message <- "Convergence due to lack of progress in parameter updates"
    }
    return(a)
}