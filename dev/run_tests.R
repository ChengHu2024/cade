suppressPackageStartupMessages({ library(testthat) })
pkg <- "C:/文献/主动脉瘤_从头重算_20261022/packages/cade"
cat("=== load_all ===\n")
ok <- tryCatch({ pkgload::load_all(pkg, quiet = TRUE); TRUE },
               error = function(e) { cat("load failed:", conditionMessage(e), "\n"); FALSE })
if (!ok) quit(status = 1)
cat("loaded\n\n=== tests ===\n")
res <- testthat::test_dir(file.path(pkg, "tests", "testthat"),
                          reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nPASS %d | FAIL %d | WARN %d | SKIP %d\n",
            sum(df$passed), sum(df$failed), sum(df$warning), sum(df$skipped)))