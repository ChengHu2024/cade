# ============================================================================
# cade 通用性压力测试
#
# 目的：不是证明"能跑"，而是找出"什么条件下会失效"。
# 每种场景都有已知答案，用 score_simulation() 量化恢复度。
# ============================================================================
.libPaths(c("C:/文献/主动脉瘤_从头重算_20261022/packages/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(cade); library(limma) })

OUT <- "C:/文献/主动脉瘤_从头重算_20261022/results/package_validation"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------- 公共设置 ----------
set.seed(20261023)
G <- 800
mu <- matrix(runif(3 * G, 2, 300), nrow = 3,
             dimnames = list(c("A", "B", "C"), paste0("g", seq_len(G))))
## 让三个基因的表达有细胞类型特异性，这样组成变化才有信号
mu["A", "g1"] <- 900; mu["B", "g1"] <- 3;   mu["C", "g1"] <- 2
mu["B", "g2"] <- 900; mu["A", "g2"] <- 3;   mu["C", "g2"] <- 2
mu["C", "g3"] <- 900; mu["A", "g3"] <- 3;   mu["B", "g3"] <- 2

obs_from_bulk <- function(bulk, is_case) {
  E <- t(log2(bulk + 1))
  tt <- topTable(eBayes(lmFit(E, model.matrix(~ as.integer(is_case)))),
                 coef = 2, number = Inf, sort.by = "none")
  setNames(tt$logFC, rownames(tt))
}

frac_of <- function(f, n) {
  m <- matrix(rep(f, n), ncol = 3, byrow = TRUE)
  colnames(m) <- c("A", "B", "C"); m
}

run_case <- function(tag, fc, fa, n_ctrl, n_case, intrinsic, noise,
                     use_frac = NULL, ref_mu = NULL, note = "") {
  sim <- simulate_pseudobulk(mu, fc, fa, n_ctrl, n_case,
                             intrinsic = intrinsic, noise = noise, seed = 7)
  obs <- obs_from_bulk(sim$bulk, sim$is_case)
  fr  <- if (is.null(use_frac)) frac_of(fa, n_case) else use_frac
  fr  <- rbind(frac_of(fc, n_ctrl), fr)
  mu_use <- if (is.null(ref_mu)) mu else ref_mu

  d <- decompose_composition(obs, fr, mu_use, sim$is_case)
  ig <- names(intrinsic)
  sc <- score_simulation(d, sim$truth, ig)

  ## 组成变化的幅度
  comp_span <- max(abs(sim$truth$composition), na.rm = TRUE)
  data.frame(
    scenario = tag,
    n_ctrl = n_ctrl, n_case = n_case,
    noise = noise,
    n_intrinsic = length(intrinsic),
    intrinsic_total_abs = if (is.null(intrinsic)) 0 else sum(abs(intrinsic)),
    max_comp_change = round(comp_span, 3),
    cor_pred = round(sc$cor_pred, 3),
    err_pred = round(sc$median_error_pred, 3),
    recovery = round(sc$recovery_intrinsic, 3),
    fpr = round(sc$false_positive_rate, 3),
    note = note, stringsAsFactors = FALSE
  )
}

cat("=============================================================\n")
cat("cade 压力测试\n")
cat("=============================================================\n\n")

BIG_C   <- c(A = .75, B = .15, C = .10)
BIG_A   <- c(A = .15, B = .70, C = .15)
SMALL_C <- c(A = .55, B = .30, C = .15)
SMALL_A <- c(A = .45, B = .38, C = .17)
SAME    <- c(A = .50, B = .35, C = .15)

RES <- list()
add <- function(...) RES[[length(RES) + 1]] <<- run_case(...)

## ---- 1. 纯组成变化，无细胞内改变 ----
add("1 纯组成变化", BIG_C, BIG_A, 25, 25, NULL, 0.05)
cat("[1] 纯组成变化        完成\n")

## ---- 2. 纯细胞内改变，无组成变化 ----
add("2 纯细胞内改变", SAME, SAME, 25, 25, c(g1 = 1.5), 0.05)
cat("[2] 纯细胞内改变      完成\n")

## ---- 3. 两者都有 ----
add("3 组成+细胞内", BIG_C, BIG_A, 25, 25, c(g3 = -1, g5 = 1), 0.05)
cat("[3] 组成+细胞内       完成\n")

## ---- 4. 组成变化很小 ----
add("4 组成变化很小", SMALL_C, SMALL_A, 25, 25, c(g1 = 1), 0.05)
cat("[4] 组成变化很小      完成\n")

## ---- 5. 高噪声 ----
add("5 高噪声(sd=0.3)", BIG_C, BIG_A, 25, 25, c(g1 = 1), 0.30)
cat("[5] 高噪声            完成\n")

## ---- 6. 样本量很小 ----
add("6 样本量小(6 vs 6)", BIG_C, BIG_A, 6, 6, c(g1 = 1), 0.05)
cat("[6] 样本量小          完成\n")

## ---- 7. 无噪声（代数检验）----
add("7 无噪声(代数检验)", BIG_C, BIG_A, 20, 20, c(g1 = 1, g2 = -.5), 0)
cat("[7] 无噪声            完成\n")

## ---- 8. 参考缺少一个细胞类型 ----
## 把 C 类细胞的表达抹成 0，模拟"参考里没有这个细胞类型的数据"
mu_bad <- mu; mu_bad["C", ] <- 0
r8 <- tryCatch(
  run_case("8 参考缺一类细胞", BIG_C, BIG_A, 25, 25, c(g1 = 1), 0.05,
           ref_mu = mu_bad, note = "C 类被清零"),
  error = function(e) data.frame(scenario = "8 参考缺一类细胞",
                                 note = paste("报错:", conditionMessage(e))))
RES[[length(RES) + 1]] <- r8
cat("[8] 参考缺一类细胞    完成\n")

## ---- 9. 疾病态参考（关键失败模式）----
## 用病例组的组成估一个"疾病态"参考：把 mu 按病例组成加权后重新分配到三类
pseudo_dis <- mu
for (k in rownames(mu)) {
  other <- setdiff(rownames(mu), k)
  pseudo_dis[k, ] <- 0.7 * mu[k, ] + 0.3 * colMeans(mu[other, , drop = FALSE])
}
add("9 疾病态参考(混杂)", BIG_C, BIG_A, 25, 25, c(g1 = 1), 0.05,
    ref_mu = pseudo_dis, note = "参考被其他类型污染 30%")
cat("[9] 疾病态参考        完成\n")

R <- do.call(rbind, RES)
rownames(R) <- NULL
write.csv(R, file.path(OUT, "stress_test.csv"), row.names = FALSE)

cat("\n=============================================================\n")
cat("结果\n")
cat("=============================================================\n")
op <- options(width = 200)
print(R[, c("scenario","cor_pred","err_pred","recovery","fpr","note")],
      row.names = FALSE)
options(op)

cat("\n=============================================================\n")
cat("判读\n")
cat("=============================================================\n")
for (i in seq_len(nrow(R))) {
  cat(sprintf("\n%s\n", R$scenario[i]))
  cat(sprintf("  组成预测相关 %.3f | 误差 %.3f | 细胞内恢复 %.3f | 假阳性 %.3f\n",
              R$cor_pred[i], R$err_pred[i], R$recovery[i], R$fpr[i]))
  if (!is.na(R$recovery[i])) {
    if (R$recovery[i] > 0.85) cat("  → 恢复良好\n")
    else if (R$recovery[i] > 0.5) cat("  → 部分恢复，幅度被低估\n")
    else cat("  → ⚠ 恢复失败，残差不代表真实细胞内改变\n")
  }
}
cat("\n输出:", OUT, "\n")
