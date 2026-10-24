# ============================================================================
# 真实性验证：加入便利层之后，数字有没有变？
#
# 这脚本做的不是"再跑一遍看结果"，而是**逐位比对**：
#   同一输入 → 旧路径 vs 新路径 → 用 identical() 比
# ============================================================================
.libPaths(c(file.path("..", ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(cade))
set.seed(20261024)

OUT <- "branding"

cat("\n")
cat("=============================================================\n")
cat("  验证 1：report() 是否改动对象（bit-identical）\n")
cat("=============================================================\n\n")

G <- 200
ct <- c("VSMC", "Macrophage", "Fibroblast")
gen <- c("MYH11", "CD68", "COL1A1", paste0("g", 1:(G-3)))
mu <- matrix(stats::runif(3*G, 5, 400), 3, dimnames = list(ct, gen))
mu["VSMC","MYH11"] <- 1200; mu["Macrophage","MYH11"] <- 5;  mu["Fibroblast","MYH11"] <- 4
mu["Macrophage","CD68"] <- 1300; mu["VSMC","CD68"] <- 4; mu["Fibroblast","CD68"] <- 5
mu["Fibroblast","COL1A1"] <- 1100; mu["VSMC","COL1A1"] <- 6; mu["Macrophage","COL1A1"] <- 4

fc <- c(VSMC=.70, Macrophage=.10, Fibroblast=.20)
fa <- c(VSMC=.15, Macrophage=.40, Fibroblast=.45)
sim <- simulate_pseudobulk(mu, fc, fa, 20, 20, intrinsic=c(COL1A1=1.2), noise=.08, seed=7)
E <- t(log2(sim$bulk + 1))
tt <- limma::topTable(limma::eBayes(limma::lmFit(E, model.matrix(~ as.integer(sim$is_case)))),
                      coef=2, number=Inf, sort.by="none")
obs <- stats::setNames(tt$logFC, rownames(tt))
frac <- rbind(matrix(rep(fc,20),20,3,byrow=TRUE), matrix(rep(fa,20),20,3,byrow=TRUE))
colnames(frac) <- ct

d <- decompose_composition(obs, frac, mu, sim$is_case)

before <- as.data.frame(d)
invisible(suppressWarnings(capture.output(report(d, calibrate=TRUE, attribute=TRUE, top_n=20))))
invisible(suppressWarnings(capture.output(check_calibration(d))))
invisible(suppressWarnings(capture.output(attribute_genes(d, min_obs=0.3))))
after <- as.data.frame(d)

cat("  渲染前 obs   的前 3 位:", sprintf("%.17g", head(before$obs, 3)), "\n")
cat("  渲染后 obs   的前 3 位:", sprintf("%.17g", head(after$obs, 3)), "\n")
cat("  渲染前 resid 的前 3 位:", sprintf("%.17g", head(before$resid, 3)), "\n")
cat("  渲染后 resid 的前 3 位:", sprintf("%.17g", head(after$resid, 3)), "\n\n")

bit_same <- identical(before, after)
cat(sprintf("  identical() = %s  → %s\n\n", bit_same,
            ifelse(bit_same, "✓ 逐位相同", "✗ 被改动了")))

cat("=============================================================\n")
cat("  验证 2：分解结果与「不用任何便利层」时是否一致\n")
cat("=============================================================\n\n")

## 完全绕过便利层，手工重算
pred_manual <- compose_prediction(frac, mu, sim$is_case, genes = names(obs))
resid_manual <- obs - pred_manual
same_pred <- isTRUE(all.equal(as.numeric(pred_manual), d$pred, tolerance = 0))
same_resid <- isTRUE(all.equal(as.numeric(resid_manual), d$resid, tolerance = 0))

cat(sprintf("  pred  手工重算 vs drop-in: %s\n", ifelse(same_pred,  "✓ 逐位相同", "✗ 不同")))
cat(sprintf("  resid 手工重算 vs drop-in: %s\n", ifelse(same_resid, "✓ 逐位相同", "✗ 不同")))
cat(sprintf("  obs 与对象一致:            %s\n\n",
            ifelse(isTRUE(all.equal(as.numeric(obs), d$obs, tolerance=0)), "✓", "✗")))

cat("=============================================================\n")
cat("  验证 3：report 里的每个数字都能在对象里找到原值\n")
cat("=============================================================\n\n")

L <- report(d, calibrate=FALSE, attribute=FALSE, top_n=10)
txt <- paste(L, collapse="\n")
ord <- order(-abs(d$resid))
sel <- d[ord[1:10], ]
ok_n <- 0
for (i in seq_len(nrow(sel))) {
  pat <- sprintf("| %s | %s | %s | %s |", sel$gene[i],
                 sprintf("%+.2f", sel$obs[i]),
                 sprintf("%+.2f", sel$pred[i]),
                 sprintf("%+.2f", sel$resid[i]))
  if (grepl(pat, txt, fixed=TRUE)) ok_n <- ok_n + 1
}
cat(sprintf("  top-10 表格行逐行匹配: %d / 10  %s\n\n", ok_n,
            ifelse(ok_n==10, "✓", "✗")))

cat("=============================================================\n")
cat("  验证 4：self_test 的第 4 条判据\n")
cat("=============================================================\n\n")
st <- self_test(verbose=FALSE)
for (i in seq_len(nrow(st$checks)))
  cat(sprintf("  [%s] %-44s %.4f\n",
              ifelse(st$checks$pass[i], "PASS", "FAIL"),
              st$checks$criterion[i], st$checks$measured[i]))
cat(sprintf("\n  self_test 总体: %s\n\n", ifelse(st$pass, "PASSED", "FAILED")))

cat("=============================================================\n")
cat("  验证 5：便利层没有调用任何计算函数（静态检查）\n")
cat("=============================================================\n\n")
files <- list.files("R",
                    pattern="0[89]-", full.names=TRUE)
strip_comments <- function(lines) {
  lines <- sub("^\\s*#'?.*$", "", lines)   # 去掉整行注释与 roxygen
  lines <- sub("#.*$", "", lines)          # 去掉行尾注释（本包字符串里没有 # ）
  lines
}
all_clean <- TRUE
for (f in files) {
  raw <- readLines(f, warn = FALSE)
  code <- paste(strip_comments(raw), collapse = "\n")
  has_compose <- grepl("compose_prediction\\s*\\(", code)
  has_decomp  <- grepl("decompose_composition\\s*\\(", code)
  hit <- has_compose || has_decomp
  if (hit) all_clean <- FALSE
  cat(sprintf("  %-18s compose_prediction(): %-5s  decompose_composition(): %-5s  %s\n",
              basename(f), has_compose, has_decomp,
              ifelse(hit, "<- 代码里有", "ok")))
}
cat("\n  只统计**非注释的代码行**。行内出现这两个名字的地方是 roxygen 文档与注释。\n\n")

cat("=============================================================\n")
verdict <- bit_same && same_pred && same_resid && ok_n == 10 && st$pass && all_clean
cat(sprintf("  总判定：%s\n", ifelse(verdict, "✓ 真实性未受影响", "✗ 有问题，需排查")))
cat("=============================================================\n")
