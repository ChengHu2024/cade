.libPaths(c(file.path("..", ".Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(cade); library(limma) })
set.seed(1)
G <- 500
mu <- matrix(runif(3*G, 1, 200), nrow=3,
             dimnames=list(c("A","B","C"), paste0("g", seq_len(G))))
mu["A","g1"] <- 500; mu["B","g1"] <- 5; mu["C","g1"] <- 2

sim <- simulate_pseudobulk(mu, c(A=.6,B=.3,C=.1), c(A=.2,B=.5,C=.3),
                           25, 25, intrinsic=c(g1=-1, g2=1), noise=0.05, seed=1)

frac <- rbind(matrix(rep(c(.6,.3,.1),25), ncol=3, byrow=TRUE),
              matrix(rep(c(.2,.5,.3),25), ncol=3, byrow=TRUE))
colnames(frac) <- c("A","B","C")

E <- t(log2(sim$bulk + 1))
tt <- topTable(eBayes(lmFit(E, model.matrix(~ as.integer(sim$is_case)))),
               coef=2, number=Inf, sort.by="none")
obs <- setNames(tt$logFC, rownames(tt))

d <- decompose_composition(obs, frac, mu, sim$is_case)

cat("=== 逐基因核对（g1, g2 是注入了改变的）===\n")
for (g in c("g1","g2","g3")) {
  i <- which(d$gene == g)
  t <- sim$truth[sim$truth$gene == g, ]
  cat(sprintf("%-3s  obs=%+7.3f  pred=%+7.3f  resid=%+7.3f | truth: comp=%+7.3f intr=%+5.2f total=%+7.3f\n",
      g, d$obs[i], d$pred[i], d$resid[i], t$composition, t$intrinsic, t$total))
}

cat("\n=== merge 之后的样子 ===\n")
m <- merge(d, sim$truth, by="gene")
cat("列名:", paste(names(m), collapse=", "), "\n")
print(head(m[, c("gene","obs","pred","resid","composition","intrinsic")], 3))

cat("\n=== 手算 recovery ===\n")
ii <- m$gene %in% c("g1","g2")
print(m[ii, c("gene","resid","intrinsic")])
cat("resid/intrinsic:", round(m$resid[ii] / m$intrinsic[ii], 4), "\n")
cat("mean =", round(mean(m$resid[ii] / m$intrinsic[ii]), 4), "\n")