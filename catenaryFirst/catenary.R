library(pracma)
# library()

N <- 51         # no. of points
L <- 2          # total length of chain
h <- L / (N-1)  # maximal length of each chain link

fnobj <- function(p) {
    sum(p[(N+1):(2*N)])  # sum(y)
}
grobj <- function(p) {
    c(rep(0, N), rep(1, N))
}

heq <- function(p) {
    c(p[1], p[N]-1, p[N+1], p[2*N])
}
heq.jac <- function(x) pracma::jacobian(heq, x)

hin <- function(p) {
    x <- p[1:N]; y <- p[(N+1):(2*N)]
    h^2 - diff(x)^2 - diff(y)^2
}
hin.jac <- function(x) pracma::jacobian(hin, x)

x0 <- seq(0, 1, length.out=N)
p0 <- c(x0, rep(0, N))

system.time(
  sol <- 0
)

x <- sol$par[1:N]; y <- sol$par[(N+1):(2*N)]

plot(c(0,1), c(-1,0), type='n')
lines(x, y, col="blue", lwd=2)
grid()

curve(0.22964*cosh((x-0.5)/0.22964)-1.02603, 0, 1,
      col="red", add=TRUE)

