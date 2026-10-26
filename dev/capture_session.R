# ============================================================================
# 捕获 CADE 的真实控制台会话（不手写任何输出，全部来自实际调用）
# ============================================================================
.libPaths(c("C:/文献/主动脉瘤_从头重算_20261022/packages/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(cade); library(limma) })
set.seed(42)

OUT <- "C:/文献/主动脉瘤_从头重算_20261022/packages/cade/branding"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "console_session.txt")

## 把一次调用连同它产生的输出一起记下来
rec <- function(cmd, expr) {
  out <- capture.output({
    res <- tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)),
                                                       class = "cap_err"))
    if (inherits(res, "cap_err")) cat("Error: ", res$msg, "\n", sep = "")
  })
  list(cmd = cmd, out = out)
}

steps <- list()

## ---------- 1 ----------
steps[[1]] <- list(hdr = "1. 载入", items = list(
  rec("library(cade)", invisible(library(cade)))
))

## ---------- 2 ----------
steps[[2]] <- list(hdr = "2. 先问她能不能用", items = list(
  rec("self_test()", self_test(verbose = TRUE))
))

## ---------- 3 ----------
G  <- 300
ct <- c("VSMC", "Macrophage", "Fibroblast")
KEY <- c("MYH11", "CD68", "COL1A1")
gen <- c(KEY, paste0("g", seq_len(G - length(KEY))))
mu <- matrix(stats::runif(3*G, 5, 400), 3, dimnames = list(ct, gen))
mu["VSMC","MYH11"] <- 1200
mu["Macrophage","CD68"] <- 1300
mu["Fibroblast","COL1A1"] <- 1100
fc <- c(VSMC = .70, Macrophage = .10, Fibroblast = .20)
fa <- c(VSMC = .15, Macrophage = .40, Fibroblast = .45)
sim <- simulate_pseudobulk(mu, fc, fa, n_ctrl = 20, n_case = 20,
                           intrinsic = c(COL1A1 = 1.2), noise = 0.08, seed = 7)
E <- t(log2(sim$bulk + 1))
tt <- limma::topTable(limma::eBayes(limma::lmFit(E, model.matrix(~ as.integer(sim$is_case)))),
                      coef = 2, number = Inf, sort.by = "none")
obs <- stats::setNames(tt$logFC, rownames(tt))
frac <- rbind(matrix(rep(fc, 20), 20, 3, byrow = TRUE),
              matrix(rep(fa, 20), 20, 3, byrow = TRUE))
colnames(frac) <- ct

steps[[3]] <- list(hdr = "3. 造一份有标准答案的数据", items = list(
  rec("mu <- matrix(runif(3*300, 5, 400), 3, dimnames = list(ct, gen))",
      cat(sprintf("%d 细胞类型 x %d 基因\n", nrow(mu), ncol(mu)))),
  rec("sim <- simulate_pseudobulk(mu, frac_ctrl = c(VSMC=.70, Macrophage=.10, Fibroblast=.20),
                   frac_case = c(VSMC=.15, Macrophage=.40, Fibroblast=.45),
                   n_ctrl = 20, n_case = 20, intrinsic = c(COL1A1 = 1.2), seed = 7)",
      cat(sprintf("bulk: %d 样本 x %d 基因\n", nrow(sim$bulk), ncol(sim$bulk)))),
  rec("sim$truth[sim$truth$gene %in% c(\"MYH11\",\"CD68\",\"COL1A1\"), ]",
      print(sim$truth[sim$truth$gene %in% c("MYH11","CD68","COL1A1"), ], row.names = FALSE))
))

## ---------- 4 ----------
d <- decompose_composition(obs, frac, mu, sim$is_case)
obs_txt <- sprintf("obs: %d 基因（由 limma 从 bulk 算出）", length(obs))
frac_txt <- sprintf("frac: %d 样本 x %d 细胞类型（由反卷积算出）", nrow(frac), ncol(frac))
steps[[4]] <- list(hdr = "4. 拆", items = list(
  rec("length(obs); dim(frac)", { cat(obs_txt, "\n"); cat(frac_txt, "\n") }),
  rec("d <- decompose_composition(obs, frac, mu, sim$is_case)", invisible(NULL)),
  rec("d", print(d)),
  rec("as.data.frame(d)[order(-abs(d$resid)), ][1:4, c(\"gene\",\"obs\",\"pred\",\"resid\")]",
      print(as.data.frame(d)[order(-abs(d$resid)), ][1:4, c("gene","obs","pred","resid")], row.names = FALSE))
))

## ---------- 5 ----------
cal <- check_calibration(d, levels = c(0, 0.5, 0.75))
steps[[5]] <- list(hdr = "5. 判断她的话能不能信", items = list(
  rec("check_calibration(d)", print(cal, row.names = FALSE)),
  rec("calibration_verdict(cal)", cat(calibration_verdict(cal), "\n"))
))

## ---------- 6 ----------
sc <- score_simulation(d, sim$truth, "COL1A1")
steps[[6]] <- list(hdr = "6. 对照真值打分", items = list(
  rec("score_simulation(d, sim$truth, \"COL1A1\")", print(sc, row.names = FALSE))
))

## ---------- 7 ----------
a <- attribute_genes(d, min_obs = 0.3)
steps[[7]] <- list(hdr = "7. 归类", items = list(
  rec("a <- attribute_genes(d, min_obs = 0.3)", invisible(NULL)),
  rec("table(a$class)", print(table(a$class)))
))

## ---------- 8 ----------
steps[[8]] <- list(hdr = "8. 给她错的输入，她拒绝", items = list(
  rec("decompose_composition(obs, frac, log2(mu + 1), sim$is_case)",
      decompose_composition(obs, frac, log2(mu + 1), sim$is_case))
))

frac_bad <- frac
colnames(frac_bad) <- c("VSMC", "VSMC", "Fibroblast")
steps[[8]]$items[[2]] <- rec(
  "colnames(frac_bad) <- c(\"VSMC\",\"VSMC\",\"Fibroblast\")  # 重复列名",
  decompose_composition(obs, frac_bad, mu, sim$is_case))

steps[[9]] <- list(hdr = "9. 参考合不合格，也先查", items = list(
  rec("check_reference(rep(ct, c(600, 30, 300)), min_cells = 100)",
      print(check_reference(rep(ct, c(600, 30, 300)), min_cells = 100)))
))

## ---------- 写出 ----------
lines <- character(0)
for (s in steps) {
  lines <- c(lines, "",
             strrep("=", 62),
             paste0("# ", s$hdr),
             strrep("=", 62))
  for (it in s$items) {
    lines <- c(lines, "", paste0("> ", it$cmd))
    if (length(it$out)) lines <- c(lines, it$out)
  }
}
writeLines(lines, LOG, useBytes = TRUE)
cat("已写出:", LOG, "\n")
cat("行数:", length(lines), "\n")
