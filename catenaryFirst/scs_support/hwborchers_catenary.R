require(Matrix, quietly=TRUE)
require(ECOSolveR, quietly=TRUE)

N <- 1001                 # 2N + 1 variables
L <- 1; h <- 2/(N-1)

c <- c(rep(0,N), rep(1,N), 0)

A <- Matrix(0, nrow=5, ncol=2*N+1, sparse=TRUE)
A[1, 2*N+1] <- 1                # x[2*N+1] = 1
A[2, 1] <- 1; A[3, N] <- 1      # x[1] = 0; x[N] = 1
A[4, N+1] <- 1; A[5, 2*N] <- 1  # y[1] = 1; y[N] = 1

b = c(h, 0, 1, 1, 1)

G <- Matrix(0, nrow=3*(N-1), ncol=2*N+1, sparse=TRUE)

for (i in 1:(N-1)) {
    j <- 3*(i-1) + 1
    G[j, 2*N+1] <- -1
    G[j+1, i] <- -1; G[j+1, i+1] <- 1
    G[j+2, N+i] <- -1; G[j+2, N+i+1] <- 1
}

H <- rep(0, 3*(N-1))

quad <- as.integer(rep(3, N-1))

ecos_b <- c(b, H)
sol <- ECOS_csolve(c, G, H, dims=list(q=quad), A, b)
str(sol$x)

## using scs
library(scs)
sol_scs <- scs(A=rbind(A, G), b=c(b, H), obj=c, cone=list(f=NROW(A), q=quad))
str(sol_scs)
max(abs(sol$x - sol_scs$x))

## using ROI
library(slam)
library(ROI)

m <- NROW(A) + NROW(G)
soc_indices <- split(seq(NROW(A)+1, NROW(A)+NROW(G)), cumsum((seq_len(NROW(G))-1) %% 3 == 0))
AA <- as.simple_triplet_matrix(as.matrix(rbind(A, G)))
op <- OP( L_objective(c), 
          L_constraint(AA, dir=rep.int("==", m), rhs=c(b, H)),
          bounds = c(V_bound(li=seq_along(c), lb=rep.int(-Inf, length(c))),
                     as.C_bound(list(free=seq_len(NROW(A)), soc=soc_indices))) )

roi_ecos <- ROI_solve(op, solver="ecos")
roi_ecos$status
roi_ecos$solution


roi_scs <- ROI_solve(op, solver="scs", tol=1e-11)
roi_scs$status   ## status is UNBOUNDED_INACCURATE therfore ROI will return NA
roi_scs$solution 
roi_scs$message$x ## solution can still be obtained from the original message
str(roi_scs)


dims_ecos <- list(q=quad)
str(dims_ecos)
dims_scs <- c(list(f=3), dims_ecos)
str(dims_scs)
