# This explores the behaviour of the different methods for non-orthogonal batch effects.

source("functions.R")

#########################################################
# Creating L-shapes in two dimensions.

set.seed(0)
ncells <- 200
all.L <- list()

for (x in 1:2) {    
    hor <- cbind(runif(ncells, min = 0, max = 2),
                 rnorm(ncells, mean = 0, sd = 0.1))
    ver <- cbind(rnorm(ncells, mean = 0, sd = 0.1), 
                 runif(ncells, min = 0, max = 2))
    all.L[[x]] <- rbind(hor, ver)
}

L1 <- all.L[[1]]
L2 <- all.L[[2]]

#########################################################
# Simulating a batch vector with orthogonal and non-orthogonal components.

shift <- matrix(c(0.3,-0.3), ncells*2, 2, byrow=TRUE)
L2b <- L2 + shift

plot(L1, col="red", pch=rep(1:2, each=ncells),
    xlim=range(c(L1[,1], L2b[,1])),
    ylim=range(c(L1[,2], L2b[,2])))
points(L2b, pch=rep(1:2, each=ncells))

ngenes <- 2000
proj <- matrix(rnorm(ngenes*2), nrow=ngenes, ncol=2)
L1p <- tcrossprod(proj, L1)
L2p <- tcrossprod(proj, L2b) + rnorm(ngenes)

out <- runAllMethods(L1p, L2p)
cluster.id <- rep(rep(1:2, each=ncells), 2)

pdf("nonorth_nonorth.pdf")
for (x in names(out$mat)) {
    plotResults(out$mat[[x]], cluster.id, out$batch, main=x,
            pch.choices=c(16, 2))
}
dev.off()

#########################################################
# Simulating a batch vector with orthogonal and rotation components.

theta <- 20*pi/180
rotM <- rbind(c(cos(theta),-sin(theta)),
              c(sin(theta),cos(theta)))

L2r <- tcrossprod(L2, rotM)
plot(L1, col="red", pch=rep(1:2, each=ncells),
    xlim=range(c(L1[,1], L2r[,1])),
    ylim=range(c(L1[,2], L2r[,2])))
points(L2r, pch=rep(1:2, each=ncells))

L1p <- tcrossprod(proj, L1)
L2p <- tcrossprod(proj, L2r) + rnorm(ngenes)

out <- runAllMethods(L1p, L2p)
cluster.id <- rep(rep(1:2, each=ncells), 2)

pdf("nonorth_rotated.pdf")
for (x in names(out$mat)) {
    plotResults(out$mat[[x]], cluster.id, out$batch, main=x,
            pch.choices=c(16, 2))
}
dev.off()

