# plotting for timing experiment
pbmc <- read.table("~/CI_filesystem//mnt/scratcha/jmlab/morgan02/10X/MNN/data/PBMC/PBMC_norm.tsv.gz",
                   h=T, sep="\t", stringsAsFactors=FALSE)
rownames(pbmc) <- pbmc$gene_id

pbmc.samples <- seq(1:10)/10
n.pbmc <- length(colnames(pbmc)) - 1

n.samples <- floor(pbmc.samples * n.pbmc)

timings <- read.table("~/CI_filesystem/mnt/scratcha/jmlab/morgan02/10X/MNN/results/68K_PBMC_timing.tsv",
                      h=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(timings) <- c("user", "system", "total", "cumulative", "child", "run")
timings$NCells <- n.samples

# plot the points and a linear fit
lm.time <- lm(user ~ NCells, data=timings)

png("Droplet/PBMC68k_timings.png",
    height=3.75, width=5.25, res=300, units="in")
par(mar=c(5.1, 4.1, 2.1, 4.1))
plot(x=timings$NCells, y=lm.time$fitted.values/60,
     lty='dashed', col='red', type="n",
     ylab="CPU time (mins)", xlab="N Cells",
     xlim=c(1, 70000), ylim=c(0, 1200/60),
     yaxt="n", bty="n")
lines(x=timings$NCells, y=lm.time$fitted.values/60,
      lty='dashed', col='red')
lines(x=timings$NCells, y=timings$user/60,
      lty=1, col='grey')
points(x=timings$NCells, y=timings$user/60,
       pch=21, bg='grey', col='black', cex=1.5)
axis(side=2,
     at=c(0, 200, 400, 600, 800, 1000, 1200)/60,
     label=floor(c(0, 200, 400, 600, 800, 1000, 1200)/60))
dev.off()

