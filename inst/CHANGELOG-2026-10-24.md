# CADE 0.9.0 便利层改动记录

**日期**：2026-10-24 01:00
**改动**：加入 3 个新函数（`diagnose_inputs`、`report`、`recommend_reference`）
**结论**：**真实性未受影响，有 5 项独立验证为证**

---

# 一、加了什么

| 函数 | 做什么 | 是否计算 |
|---|---|---|
| `diagnose_inputs()` | 把 `check_inputs()` 的每条错误/警告配上一句"怎么办" | ❌ 只做字符串匹配 |
| `report()` | 一次输出可直接贴进论文的 Markdown 报告（含校准判读、top 基因、归类、稳定性、已知答案打分） | ❌ 只转印已有对象 |
| `recommend_reference()` | 根据两套参考的**出处事实**（有无正常组织、细胞数是否充足）给建议 | ❌ 只读 `compare_references()` 的输出 |

**同时**：`self_test()` 从 3 条判据增加到 4 条，第 4 条是"渲染报告前后对象是否逐位相同"。

---

# 二、你问的那句："改动之后跑数据真实性还能得到保障吗？"

## 2.1 五重验证（可重跑）

脚本：`packages/cade/dev/verify_authenticity.R`

### 验证 1：`report()` 不改动对象

```
渲染前 obs 前 3 位: -2.2616096070181566  1.9619146669038268  2.3375418210240713
渲染后 obs 前 3 位: -2.2616096070181566  1.9619146669038268  2.3375418210240713

渲染前 resid 前 3 位: -0.067126029365297324  -0.0044043293726041544  1.1902032231062987
渲染后 resid 前 3 位: -0.067126029365297324  -0.0044043293726041544  1.1902032231062987

identical(before, after) = TRUE   →  逐位相同
```

### 验证 2：与"完全不用便利层"的结果一致

```
pred   手工重算 vs drop-in:  逐位相同
resid  手工重算 vs drop-in:  逐位相同
obs    与对象一致:           逐位相同
```

### 验证 3：报告里的每个数字都能在对象里找到原值

```
top-10 表格行逐行匹配:  10 / 10
```

### 验证 4：`self_test()` 的第 4 条判据

```
[PASS] composition term recovers the truth          1.0000
[PASS] injected within-cell change recovered        0.9984
[PASS] no false positives above tolerance           0.0000
[PASS] reporting does not alter the decomposition   1.0000
```

### 验证 5：静态检查——便利层不调用任何计算函数

```
08-remedies.R      compose_prediction(): FALSE  decompose_composition(): FALSE  ok
09-report.R        compose_prediction(): FALSE  decompose_composition(): FALSE  ok
```

**只统计非注释的代码行。** 行内出现这两个名字的地方是 roxygen 文档与注释。

**总判定：真实性未受影响。**

## 2.2 结构性保证（不靠自觉）

| 机制 | 作用 |
|---|---|
| `09-report.R` 里 `report()` 只读 `decomp` 的列，从不调用 `compose_prediction()` | 结构上不可能改数字 |
| 测试 `report numbers are transcription, not recomputation` | 一旦有人让 `report()` 重算，测试立刻失败 |
| `self_test()` 第 4 条 | 安装后一条命令即可自证 |
| `report()` 的文档明写"transcription, not computation" | 使用者能看到这个承诺 |

---

# 三、包的其他变化

## 3.1 修掉的三个问题

| 问题 | 修法 |
|---|---|
| 中文 remed 字符串破坏 ASCII 可移植性（R CMD check WARNING） | 包内全部改英文；中文版单独出一份 `inst/CADE-CN-cheatsheet.md` |
| `capture.output` 未导入（NOTE） | 加入 `importFrom(utils, capture.output)` |
| `branding/` 与 `dev/` 被打进 tarball（840 KB） | 加 `.Rbuildignore`，包降到 **52 KB** |

## 3.2 当前状态

```
R 源文件     10 个 / 1,782 行
man 文档     21 个
测试         143 个，全部通过
R CMD check  Status: OK（0 ERROR / 0 WARNING / 0 NOTE）
包大小       52 KB
```

---

# 四、便利层的中文对照

包内是英文（可移植性要求）。中文速查表在：

```
packages/cade/inst/CADE-CN-cheatsheet.md
packages/cade/branding/CADE_中文速查表.md
```

**逐条对应，不是另写一套。**

---

# 五、还有一件没做

**`recommend_reference()` 的"出处事实"需要使用者自己提供**（`a_has_normal`、`a_cells_per_type`）。

**这是有意的**：那些事实不能从数据里推断，只能从参考是怎么建的来知道。**如果让它从数据猜，就又变成循环论证了。**

---

# 六、给明天的你

```r
library(cade)

self_test()                                    # 4 条判据，1 秒
diagnose_inputs(obs, frac, mu, is_case)        # 有错先看这个
d <- decompose_composition(obs, frac, mu, is_case)
report(d, file = "cade_report.md")             # 一次出全
```

**跑真实数据时不需要做任何额外的事——新加的三个函数都是可选的。**
