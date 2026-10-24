# ============================================================================
# 把真实控制台会话渲染成图（两栏：正常工作 / 拒绝）
# ============================================================================
.libPaths(c(file.path("..", ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

OUT <- "branding"
TXT <- file.path(OUT, "console_session.txt")

raw <- readLines(TXT, encoding = "UTF-8", warn = FALSE)

## 太长的行按 92 字符折行（控制台里用缩进续行表示）
wrap_line <- function(s, width = 92) {
  if (nchar(s, type = "width") <= width) return(s)
  pieces <- character(0)
  rest <- s
  first <- TRUE
  while (nchar(rest, type = "width") > width) {
    cut <- width
    ## 尽量在空格处断开
    seg <- substr(rest, 1, cut)
    sp <- max(gregexpr(" ", seg, fixed = TRUE)[[1]])
    if (sp > 40) cut <- sp
    pieces <- c(pieces, substr(rest, 1, cut))
    rest <- substring(rest, cut + 1)
    first <- FALSE
  }
  pieces <- c(pieces, rest)
  ## 续行加缩进，视觉上属于上一条命令
  c(pieces[1], paste0("    ", pieces[-1]))
}
raw <- unlist(lapply(raw, wrap_line), use.names = FALSE)

## ---------- 切分成两段 ----------
i8 <- grep("^# 8\\.", raw)[1]
A <- raw[1:(i8 - 3)]
B <- raw[(i8 - 2):length(raw)]

## ---------- 着色规则 ----------
classify <- function(x) {
  ifelse(grepl("^>", x),                                  "cmd",
  ifelse(grepl("^Error", x),                              "err",
  ifelse(grepl("^\\s*\\*\\s", x),                          "err2",
  ifelse(grepl("PASS|PASSED|informative", x),             "ok",
  ifelse(grepl("^#|^=+$", x),                             "hdr",
  ifelse(grepl("^(\\s*)$", x),                            "blank", "out"))))))
}

COLS <- c(cmd = "#7ec8e3", err = "#ff7b72", err2 = "#ffa198",
          ok  = "#7ee787", out = "#c9d1d9", hdr = "#d2a8ff", blank = "#0d1117")

BG   <- "#0d1117"
FONT <- "Consolas"

draw_panel <- function(lines, title, subtitle) {
  d <- data.frame(y = seq_along(lines), txt = lines,
                  cls = classify(lines), stringsAsFactors = FALSE)
  d$y <- rev(d$y)
  n <- nrow(d)
  ggplot(d, aes(x = 0, y = y)) +
    geom_text(aes(label = txt, colour = cls), family = FONT, size = 2.65,
              hjust = 0, vjust = 0.5, lineheight = 1.0) +
    scale_colour_manual(values = COLS, guide = "none") +
    annotate("text", x = 0, y = n + 3.2, label = title, hjust = 0,
             family = FONT, size = 4.2, fontface = 2, colour = "#ffffff") +
    annotate("text", x = 0, y = n + 1.6, label = subtitle, hjust = 0,
             family = FONT, size = 2.7, colour = "#8b949e") +
    coord_cartesian(xlim = c(0, 96), ylim = c(0, n + 4.6), expand = FALSE) +
    theme_void() +
    theme(plot.background = element_rect(fill = BG, colour = NA),
          plot.margin = margin(14, 14, 14, 14))
}

pA <- draw_panel(A, "A · 正常工作",
                 "cade 0.9.0   ·   以下每一行都是真实控制台输出，没有一行是手写的")
pB <- draw_panel(B, "B · 她会拒绝",
                 "同一批数据，只把输入改坏两个地方")

hA <- max(4, length(A)) * 0.106 + 1.4
hB <- max(4, length(B)) * 0.106 + 1.4

ggsave(file.path(OUT, "CADE_session_A.png"), pA,
       width = 9.2, height = hA, dpi = 260, bg = BG)
ggsave(file.path(OUT, "CADE_session_B.png"), pB,
       width = 9.2, height = hB, dpi = 260, bg = BG)
cat("A:", length(A), "行  → CADE_session_A.png\n")
cat("B:", length(B), "行  → CADE_session_B.png\n")
