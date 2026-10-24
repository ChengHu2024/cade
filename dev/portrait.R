# ============================================================================
# CADE 肖像 v2
# ============================================================================
.libPaths(c(file.path("..", ".Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(ggplot2) })
OUT <- "branding"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FAM   <- "Microsoft YaHei"
INK   <- "#16232e"
BLUE  <- "#2c6e9e"
SAND  <- "#cbb894"
RED   <- "#a83a42"
GREEN <- "#37684e"
PALE  <- "#f7f5f0"
GREY  <- "#6b7885"

p <- ggplot() + theme_void(base_family = FAM) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  theme(plot.background = element_rect(fill = PALE, colour = NA),
        plot.margin = margin(0, 0, 0, 0))

## ---------- 顶带 ----------
p <- p +
  annotate("rect", xmin = 0, xmax = 100, ymin = 87, ymax = 100,
           fill = BLUE, colour = NA) +
  annotate("text", x = 4, y = 96.4, label = "CADE", hjust = 0,
           family = FAM, size = 20, fontface = 2, colour = "white") +
  annotate("text", x = 4, y = 91.2,
           label = "Composition Attribution of Differential Expression",
           hjust = 0, family = FAM, size = 5.0, colour = "#cfe3f0") +
  annotate("text", x = 96, y = 95.6, label = "差异表达的组成归因",
           hjust = 1, family = FAM, size = 5.6, colour = "white") +
  annotate("text", x = 96, y = 91.0, label = "v0.9.0   ·   2026-10-23   ·   MIT",
           hjust = 1, family = FAM, size = 3.3, colour = "#a8c6dd")

## ---------- 核心等式 ----------
p <- p +
  annotate("rect", xmin = 4, xmax = 96, ymin = 72, ymax = 84,
           fill = "white", colour = SAND, linewidth = 0.8) +
  annotate("text", x = 50, y = 81.4, label = "她只做一件事",
           family = FAM, size = 3.9, colour = BLUE, fontface = 2) +
  annotate("text", x = 50, y = 77.0,
           label = "log2FC 观测   =   log2FC 组成   +   残差（细胞自己）",
           family = FAM, size = 5.6, colour = INK, fontface = 2) +
  annotate("text", x = 50, y = 74.0,
           label = "观测到的变化  =  细胞搬走了  +  细胞变了",
           family = FAM, size = 3.4, colour = GREY)

## ---------- 三个模块 ----------
mx <- c(18, 50, 82)
mt <- c("拆开了什么", "能不能信", "没骗你")
ms <- c("分解", "诊断", "验证")
mn <- c("5 个函数", "7 个函数", "5 个函数")
mc <- c(BLUE, GREEN, RED)
md <- list(
  c("compose_prediction()", "decompose_composition()", "expression_filter()",
    "check_inputs()", "attribute_sensitivity()"),
  c("check_calibration()", "calibration_verdict()", "check_reference()",
    "compare_references()", "attribute_genes()", "cohens_d()", "pool_cohens_d()"),
  c("simulate_pseudobulk()", "score_simulation()", "self_test()",
    "plot_calibration()", "plot_decomposition()")
)

for (i in 1:3) {
  p <- p +
    annotate("rect", xmin = mx[i]-14, xmax = mx[i]+14, ymin = 35, ymax = 68,
             fill = "white", colour = mc[i], linewidth = 1.0) +
    annotate("rect", xmin = mx[i]-14, xmax = mx[i]+14, ymin = 62.5, ymax = 68,
             fill = mc[i], colour = NA) +
    annotate("text", x = mx[i], y = 65.2, label = mt[i],
             family = FAM, size = 4.4, fontface = 2, colour = "white") +
    annotate("text", x = mx[i], y = 59.6, label = ms[i],
             family = FAM, size = 3.6, colour = mc[i], fontface = 2) +
    annotate("text", x = mx[i], y = 57.0, label = mn[i],
             family = FAM, size = 3.0, colour = GREY)
  yy <- 53.8
  for (fn in md[[i]]) {
    p <- p + annotate("text", x = mx[i], y = yy, label = fn,
                      family = FAM, size = 2.75, colour = INK)
    yy <- yy - 2.75
  }
}

## ---------- 压力测试 ----------
p <- p +
  annotate("text", x = 4, y = 31.5, label = "她经得起什么", hjust = 0,
           family = FAM, size = 4.2, fontface = 2, colour = INK) +
  annotate("text", x = 4, y = 29.0,
           label = "107 个单元测试   ·   R CMD check: Status OK   ·   8 个已知答案场景",
           hjust = 0, family = FAM, size = 3.1, colour = GREY)

tt <- data.frame(
  x = c(12, 34, 56, 78),
  lab = c("纯组成变化\n→ 1.000", "纯细胞内改变\n→ 0.998",
          "高噪声 sd=0.3\n→ 0.950", "仅 6 vs 6 样本\n→ 0.961"),
  col = c(BLUE, BLUE, BLUE, BLUE)
)
p <- p +
  annotate("rect", xmin = 2, xmax = 22, ymin = 18, ymax = 28, fill = "white",
           colour = BLUE, linewidth = 0.55) +
  annotate("rect", xmin = 24, xmax = 44, ymin = 18, ymax = 28, fill = "white",
           colour = BLUE, linewidth = 0.55) +
  annotate("rect", xmin = 46, xmax = 66, ymin = 18, ymax = 28, fill = "white",
           colour = BLUE, linewidth = 0.55) +
  annotate("rect", xmin = 68, xmax = 88, ymin = 18, ymax = 28, fill = "white",
           colour = BLUE, linewidth = 0.55) +
  annotate("rect", xmin = 90, xmax = 98, ymin = 18, ymax = 28, fill = "#fdf3f3",
           colour = RED, linewidth = 0.55) +
  annotate("text", x = 12, y = 23, label = tt$lab[1], family = FAM,
           size = 3.0, colour = BLUE, lineheight = 1.35) +
  annotate("text", x = 34, y = 23, label = tt$lab[2], family = FAM,
           size = 3.0, colour = BLUE, lineheight = 1.35) +
  annotate("text", x = 56, y = 23, label = tt$lab[3], family = FAM,
           size = 3.0, colour = BLUE, lineheight = 1.35) +
  annotate("text", x = 78, y = 23, label = tt$lab[4], family = FAM,
           size = 3.0, colour = BLUE, lineheight = 1.35) +
  annotate("text", x = 94, y = 23, label = "参考缺一类\n假阳性 16%", family = FAM,
           size = 2.7, colour = RED, lineheight = 1.35)

## ---------- 拒绝 ----------
p <- p +
  annotate("text", x = 50, y = 12.6,
           label = "她最常做的一件事，是拒绝",
           family = FAM, size = 3.6, colour = RED, fontface = 2) +
  annotate("text", x = 50, y = 9.2,
           label = "composition prediction is uninformative: slope 0.03, R2 0.017",
           family = FAM, size = 2.9, colour = "#7a4a52", fontface = 3) +
  annotate("text", x = 50, y = 6.6,
           label = "—  do not interpret the decomposition",
           family = FAM, size = 2.9, colour = "#7a4a52", fontface = 3) +
  annotate("text", x = 50, y = 2.4,
           label = "给胡骋　·　2026 年 10 月 24 日 00:20",
           family = FAM, size = 3.0, colour = "#9aa5ae")

ggsave(file.path(OUT, "CADE_portrait.png"), p, width = 10, height = 7.2,
       dpi = 300, bg = PALE, type = "cairo")
cat("saved png\n")
