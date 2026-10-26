# CADE 中文速查表

**版本**：0.9.0（2026-10-24）
**对应英文包**：`cade_0.9.0.tar.gz`，`R CMD check`: Status OK，143 个测试

> 包内的错误信息与修复建议是英文的（R 包必须保持 ASCII 可移植性）。
> 这份文件是它们的中文对照。**内容与英文版逐条对应，不是另写一套。**

---

# 一、她做事的顺序

```
1. library(cade)                     载入
2. self_test()                       先证明自己算得对（4 条判据）
3. simulate_pseudobulk()             造一份有标准答案的数据（可选，但建议）
4. decompose_composition()           拆
5. check_calibration()               判断自己说的话能不能信
6. score_simulation()                对照真值打分（若有标准答案）
7. attribute_genes()                 归类
8. report()                          出报告
```

**第 2 步和第 5 步是"证明"和"判断"，只有第 4 步在真正干活。**

---

# 二、输入检查：她说什么，你该改什么

## 2.1 参考矩阵 μ

| 她说 | 中文意思 | 你该怎么做 |
|---|---|---|
| `the largest per-cell mean in mu is X` | μ 的最大值太小，**很可能已经是 log 尺度** | 还原成线性：`log2(CP10k+1)` → `mu <- 2^L - 1`；`log2(CP10k)` → `mu <- 2^L`；`ln(CP10k+1)` → `mu <- expm1(L)`。**还原后最大值应达到数百以上**（线粒体/核糖体基因总是很高）|
| `looks like raw counts` | 这是原始计数，不是每细胞均值 | 先按文库大小归一化：`cp10k <- sweep(counts, 2, colSums(counts), '/') * 1e4`；再按类型平均 |
| `has no cell-type row names` | μ 的行没有细胞类型名 | `rownames(mu) <- 细胞类型标签` |
| `has duplicated cell-type row names` | 行名重复 | 合并：`mu <- rowsum(mu, group = rownames(mu)) / ...` |
| `all-zero expression in mu` | 某类细胞完全没有表达数据 | 要么从 frac 和 μ 里一起删掉，要么换参考 |
| `contains non-finite` | 有 NA / Inf | 先处理掉，不要留着 |

## 2.2 细胞比例 frac

| 她说 | 你该怎么做 |
|---|---|
| `has no cell-type column names` | `colnames(frac) <- 反卷积输出的细胞类型名`（不能是 V1/V2）|
| `has duplicated cell-type names` | 先判断是真重复还是命名冲突。命名冲突用 `make.unique()`；同一群体就合并成一行 |
| `row(s) of frac sum to zero` | 删掉这些样本，**同时同步删 obs 与 is_case** |
| `contains negative` | NNLS/CIBERSORT 的负值先截断到 0 再归一化 |

## 2.3 观测值 obs

| 她说 | 你该怎么做 |
|---|---|
| `is not numeric` | `as.numeric(obs)`，但**先确认它真的是 log2 倍数变化，不是线性比值** |
| `has no names` | `names(obs) <- 基因名`，否则无法与 μ 对齐 |
| `duplicated gene names` | 先用 `limma::avereps()` 或 `rowsum()` 把同名探针折叠成一个基因 |
| `contains non-finite` | 删掉 NA 行 |

## 2.4 分组 is_case

| 她说 | 中文意思 | 你该怎么做 |
|---|---|---|
| **`contains NA`** | **最危险的一条**：NA 会让下游**全部**变成 NA，而且**不报错**——所有基因都会被标成"无变化" | `is_case <- !is.na(is_case) & is_case`，或直接删掉这些样本 |
| `is not logical` | 不是逻辑值 | `is_case <- as.logical(is_case)`；若是 0/1，用 `is_case == 1` |
| `has length X but frac has Y` | 长度不一致 | 检查是否某步删了样本但没同步 |
| `no TRUE` / `no FALSE` | 只有一组 | 只有病例的数据做不了组成归因；若是方向反了，试 `is_case <- !is_case` |
| `smaller group` | 较小一组样本太少 | 每组至少 5 例；达不到就只报方向、不报幅度 |

## 2.5 对齐

| 她说 | 你该怎么做 |
|---|---|
| `present in mu but not in frac` | μ 有而反卷积没产出这类细胞 → 换反卷积设置或删掉这类 |
| `present in frac but not in mu` | 反卷积产出了 μ 里没有的类型 → 参考缺这一类 |
| `differ in case, space or encoding` | 只是命名不一致，用 `trimws()` + `tolower()` 统一 |

---

# 三、校准怎么读

## 3.1 三档斜率

```
全部基因      理想 → 0.3–0.5 属正常（低表达基因噪声大）
前 50% 表达   理想 → 0.7 以上
前 25% 表达   理想 → 0.8–1.4
```

**判据是"前 25%"这一档，不是全部基因。**

## 3.2 三种判读

| 判读 | 中文 | 你能说什么 |
|---|---|---|
| `informative` | 组成项可用 | 可以解释残差 |
| `weakly calibrated` | 弱校准 | **只报方向，不报幅度** |
| `uninformative` | 组成项不可用 | **不要解释残差** |

## 3.3 绝对不能做的事

**不要用 `|logFC| > 0.5` 去筛选基因再报校准度。**

那是循环论证——用结果筛结果。**只能用表达量百分位筛。**

---

# 四、归类怎么读

## 4.1 五类

| 类 | 含义 |
|---|---|
| `composition-driven` | 变化主要由细胞组成解释 |
| `mixed` | 两者都有 |
| `cell-intrinsic` | 变化主要发生在细胞内 |
| **`over-predicted`** | **残差与观测反号**——组成项预测得比实际更大，此时有符号的 pct 不是有意义的份额 |
| `unchanged` | 观测接近 0 |

## 4.2 阈值敏感性

`report()` 会自动算：在 6 个阈值下标签是否稳定。

**如果某个基因"随阈值变类"，就不能把它当成稳健发现写进论文。**

---

# 五、双参考怎么选

| 情形 | 她给的建议 |
|---|---|
| 只有一套参考有正常组织 | 用它；另一套只作敏感性分析 |
| 两套都有正常组织，但结论不一致 | **不一致本身就是发现**——两套都报，不要选 |
| 两套的出处都不知道 | 无法建议。**先去查参考是拿什么组织建的** |

## 最重要的一条

**不要因为某套参考的校准斜率更接近 1 就选它。**

**那个判据是循环的**——你用输出选输入，等于用结论选前提。

---

# 六、快速上手

```r
library(cade)

## 1. 先自检（约 1 秒）
self_test()

## 2. 准备三个东西
##    obs     带基因名的 log2 倍数变化
##    frac    样本 × 细胞类型 的比例矩阵
##    mu      细胞类型 × 基因 的每细胞均值（线性）
##    is_case 逻辑向量

## 3. 有错就先诊断
diagnose_inputs(obs, frac, mu, is_case)

## 4. 拆
d <- decompose_composition(obs, frac, mu, is_case)

## 5. 判断能不能信
check_calibration(d) |> calibration_verdict()

## 6. 出报告
report(d, file = "cade_report.md")
```

---

# 七、她的四条硬边界（代码修不了）

| 边界 | 后果 |
|---|---|
| **参考被污染 10–20% 时她察觉不到** | 校准仍说可用，但残差已不可信。**只能靠知道参考怎么来的** |
| **低表达基因的假阳性取决于 log 变换形式** | 用 `log2(bulk+1)` 会抬高假阳性；按表达量分层报告 |
| **细胞类型特异的改变被稀释** | 同一个改变，在何处高表达可差 11 倍 |
| **反卷积的系统性偏差被残差吸收** | 反卷积偏差与细胞内改变在本框架下不可区分 |

**这四条已写进论文 β 的 §5.3b，作为本文自己方法的定量边界公开。**
