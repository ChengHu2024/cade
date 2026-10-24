.libPaths(c(file.path("..", ".Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(cade); library(limma) })

cat("=== 已安装包冒烟测试（500 基因）===\n\n")
set.seed(1)
G <- 500
## 3 类细胞的每细胞均值，线性尺度
mu <- matrix(runif(3 * G, 1, 200), nrow = 3,
             dimnames = list(c("A","B","C"), paste0("g", seq_len(G))))
## 让 g1 在 A 类里很高、在病例中大幅减少
mu["A","g1"] <- 500; mu["B","g1"] <- 5; mu["C","g1"] <- 2

sim <- simulate_pseudobulk(mu,
        frac_ctrl = c(A=.6, B=.3, C=.1),
        frac_case = c(A=.2, B=.5, C=.3),
        n_ctrl = 25, n_case = 25,
        intrinsic = c(g1 = -1, g2 = 1), noise = 0.05, seed = 1)
cat("1. 伪 bulk:", nrow(sim$bulk), "x", ncol(sim$bulk), "\n")

frac <- rbind(matrix(rep(c(.6,.3,.1), 25), ncol=3, byrow=TRUE),
              matrix(rep(c(.2,.5,.3), 25), ncol=3, byrow=TRUE))
colnames(frac) <- c("A","B","C")

E <- t(log2(sim$bulk + 1))
tt <- topTable(eBayes(lmFit(E, model.matrix(~ as.integer(sim$is_case)))),
               coef = 2, number = Inf, sort.by = "none")
obs <- setNames(tt$logFC, rownames(tt))

d <- decompose_composition(obs, frac, mu, sim$is_case)
cat("2. 分解:", nrow(d), "基因\n\n")

sc <- score_simulation(d, sim$truth, c("g1","g2"))
cat("3. 恢复度\n")
cat("   composition 相关  :", round(sc$cor_pred, 4), "\n")
cat("   composition 误差  :", round(sc$median_error_pred, 4), "\n")
cat("   intrinsic 恢复率  :", round(sc$recovery_intrinsic, 3), "\n")
cat("   假阳性率(<0.3)    :", round(sc$false_positive_rate, 4), "\n\n")

cal <- check_calibration(d)
cat("4. 校准\n")
print(round(cal[, c("level","n","slope","r2")], 3))
cat("   判定:", calibration_verdict(cal), "\n\n")

attr <- attribute_genes(d, min_obs = 0.5)
cat("5. 归类（|obs|>0.5 才判）\n")
print(table(attr$class))
cat("\n")

gen <- c("g1","g2")
cat("6. g1/g2 的分解\n")
sub <- d[d$gene %in% gen, c("gene","obs","pred","resid")]
for (k in seq_len(nrow(sub)))
  cat(sprintf("   %-3s obs=%+7.3f pred=%+7.3f resid=%+7.3f\n",
      sub$gene[k], sub$obs[k], sub$pred[k], sub$resid[k]))
cat("   真值 intrinsic: g1 = -1, g2 = +1\n\n")

r <- check_reference(rep(c("A","B","C"), c(500,100,5)), min_cells = 50)
cat("7. 参考检查: ok =", r$ok, "\n   ", paste(r$flags, collapse="; "), "\n\n")

## 错误路径检查
cat("8. 错误路径（必须报错）\n")
e1 <- tryCatch({ compose_prediction(frac[, "A", drop=FALSE], mu, sim$is_case); "无报错" },
               error = function(e) "正确报错")
cat("   缺细胞类型:", e1, "\n")
e2 <- tryCatch({ decompose_composition(unname(obs), frac, mu, sim$is_case); "无报错" },
               error = function(e) "正确报错")
cat("   obs 无名:  ", e2, "\n\n")

cat("=== 冒烟测试完成 ===\n")