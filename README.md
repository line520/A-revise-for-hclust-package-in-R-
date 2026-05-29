# A-revise-for-hclust-package-in-R-
# 手工层次聚类：修正 `hclust` 的可靠替代

本项目基于 **Lance–Williams 通用递推公式**，从零实现了层次聚类的六种常用方法（最短距离、最长距离、类平均、重心法、中间距离法、Ward 法）。  
与 R 自带的 `stats::hclust` 相比，本实现更加透明、可控，且修正了部分方法（如重心法和 Ward 法）在原始包中可能存在的计算偏差问题。

> **背景**：R 语言（≤ 4.6.0 及更早版本）中的 `hclust` 函数在实现某些方法（尤其是 `centroid` 和 `median`）时，**距离矩阵的平方/非平方转换逻辑不一致**，导致聚类结果与教科书公式不符。  
> 本仓库从头编写了 `hclust_manual` 函数，严格遵循 Lance–Williams 递推公式，并正确处理距离的平方转换，确保所有方法的输出完全符合理论预期。

---

## 📐 数学原理

### 1. Lance–Williams 递推公式

设类 $A$ 和 $B$ 合并为新类 $C = A \cup B$，$|A|, |B|$ 为类的大小（样本数）。对于任意另一个类 $X$，新类 $C$ 与 $X$ 之间的距离 $d(C, X)$ 可以由 $d(A, X)$、$d(B, X)$ 以及 $d(A, B)$ 递推得到：

$$
d(C, X) = \alpha_A d(A, X) + \alpha_B d(B, X) + \beta \, d(A, B) + \gamma |d(A, X) - d(B, X)|
$$

对于**非加权**的六种方法，$\gamma = 0$，系数 $(\alpha_A, \alpha_B, \beta)$ 如下表所示：

| 方法        | $\alpha_A$               | $\alpha_B$               | $\beta$                     | 备注                           |
|-------------|--------------------------|--------------------------|-----------------------------|--------------------------------|
| single      | $1/2$                    | $1/2$                    | $0$                         | 同时取 $\gamma = -1/2$ 可写为 $d(C,X)=\min(d(A,X),d(B,X))$ |
| complete    | $1/2$                    | $1/2$                    | $0$                         | 同时取 $\gamma = +1/2$ 可得最大值   |
| average     | $\frac{|A|}{|A|+|B|}$    | $\frac{|B|}{|A|+|B|}$    | $0$                         |                                   |
| centroid    | $\frac{|A|}{|A|+|B|}$    | $\frac{|B|}{|A|+|B|}$    | $-\frac{|A||B|}{(|A|+|B|)^2}$ | **必须使用平方距离**             |
| median      | $1/2$                    | $1/2$                    | $-1/4$                      | **必须使用平方距离**             |
| ward        | $\frac{|A|+|X|}{|A|+|B|+|X|}$ | $\frac{|B|+|X|}{|A|+|B|+|X|}$ | $-\frac{|X|}{|A|+|B|+|X|}$  | **必须使用平方距离**             |

> 对于 `centroid`, `median`, `ward` 三种方法，Lance–Williams 公式要求初始距离矩阵为**欧氏距离的平方**。  
> 本函数自动识别这些方法，在内部将输入距离平方后再递推，最后将合并高度（`height`）开方还原为原始距离，从而保证结果的正确性。

### 2. 算法流程

1. 输入距离矩阵 `D`（$n \times n$，对称，对角为 0）。
2. 若方法为 `centroid`, `median`, `ward`，则令 `D = D^2`。
3. 初始化：每个样本自成一类，类大小为 1。
4. 重复 $n-1$ 次：
   - 在当前的类间距离矩阵中找出最小距离及对应的两个类 $p, q$。
   - 记录合并信息（`merge`）和合并高度（`height`）。
   - 使用 Lance–Williams 公式计算新类与其余各类的距离。
   - 更新距离矩阵、类大小和成员列表。
5. 若第 2 步做了平方转换，则将 `height` 开方还原。
6. 返回一个与 `stats::hclust` 结构相似的列表。

---

## 🚀 使用方法

### 环境要求
- **R 版本**：4.0 及以上（推荐 4.6.0）
- **依赖包**：仅需基础 R（无需额外安装）

### 文件说明
- `my_hclust.R`：核心函数文件，包含 `hclust_manual` 和 `plot_dendrogram`。
- `ex6.1.R`：示例脚本 – 一维数据的六种方法对比。
- `ex6.2.R`：示例脚本 – 美国 10 个城市飞行距离的四种方法聚类。

### 快速开始

```r
# 1. 加载函数
source("my_hclust.R")

# 2. 准备距离矩阵（例如 5 个点的一维坐标）
x <- c(1, 2, 6, 8, 11)
D <- as.matrix(dist(x))

# 3. 进行聚类（以 Ward 法为例）
hc <- hclust_manual(D, method = "ward", labels = paste0("G", 1:5))

# 4. 查看合并过程
print(hc$merge)
print(hc$height)

# 5. 绘制树状图
plot_dendrogram(hc, main = "Ward 聚类结果")
```

### 函数说明

#### `hclust_manual(D, method, labels = NULL)`

| 参数     | 类型        | 说明                                                                 |
|----------|-------------|----------------------------------------------------------------------|
| `D`      | 距离矩阵    | $n \times n$ 对称矩阵，非负，对角为 0。                              |
| `method` | 字符串      | 可选：`"single"`, `"complete"`, `"average"`, `"centroid"`, `"median"`, `"ward"`。 |
| `labels` | 字符向量    | 样本标签，长度 = $n$。默认为 `as.character(1:n)`。                  |

**返回值**：一个列表，包含以下组件：
- `merge`：$ (n-1) \times 2$ 矩阵，每行表示合并的两个类编号（负数表示原始样本）。
- `height`：长度为 $n-1$ 的向量，每次合并时的距离（原始欧氏距离）。
- `labels`：样本标签。
- `method`：使用的聚类方法。
- `members`：列表，记录每个（中间）类包含的原始样本索引。

#### `plot_dendrogram(hc, main = "", cex.label = 0.8)`

- `hc`：由 `hclust_manual` 返回的对象。
- `main`：图形标题。
- `cex.label`：Y 轴标签字体大小。

---

## 📊 示例演示

### 示例 1：一维数据（ex6.1.R）

5 个点：$x = [1, 2, 6, 8, 11]$。  
比较六种方法的合并过程及树状图。运行后会在控制台打印每一步的合并详情，并逐一显示树状图。

```r
source("ex6.1.R")
```

### 示例 2：美国 10 个城市的飞行距离（ex6.2.R）

使用四种方法（最短距离、类平均、重心法、Ward 法）对城市进行聚类。  
距离矩阵为真实飞行里程（英里）。

```r
source("ex6.2.R")
```

---

## ⚠️ 注意事项

1. **距离矩阵必须满足三角不等式**（对于 centroid / median / Ward 方法，需为欧氏距离的平方）。本函数不检查输入是否为有效距离，请确保 `D` 来自 `dist()` 或合理构造。
2. **内存占用**：由于递推过程每次重构距离矩阵，当 $n > 1000$ 时可能较慢。建议 $n \leq 500$。
3. **与原版 `hclust` 的差异**：
   - 对于 `centroid` 和 `median` 方法，原始 `hclust` 有时不恰当地在平方距离和非平方距离之间切换。本实现严格按照 Lance–Williams 系数，并在最终输出时将高度转换回原始距离。
   - 原版 `hclust` 的 `ward.D` 方法对应本实现的 `method = "ward"`（平方距离递推）。注意不要混淆 `ward.D2`。
4. **绘图功能**：`plot_dendrogram` 为基础 R 图形，简单但够用。如需更美观的图形，可将 `hc` 对象转换为 `dendrogram` 后使用 `ggplot2` 或 `ggdendro`。

---

## 📝 参考文献

- Lance, G. N., & Williams, W. T. (1967). A general theory of classificatory sorting strategies. *The Computer Journal*, 9(4), 373–380.
- Murtagh, F. (1985). Multidimensional clustering algorithms. *Physica-Verlag*.
- R Core Team (2025). R: A language and environment for statistical computing. (Version 4.6.0)

---

## 🧪 测试与验证

我们使用 R 自带的 `USArrests` 数据集和模拟高斯混合数据，将本函数的聚类结果与 `stats::hclust` 进行了对比：
- 对于 `single`, `complete`, `average` 方法，两者结果完全一致。
- 对于 `centroid`, `median`, `ward` 方法，本函数修正了原始版本中的错误，结果与教科书公式及 `fastcluster` 包（正确实现）相同。

欢迎 fork 并提交 issue 或 pull request。

---

## 📄 许可证

MIT License – 可自由使用、修改和分发，请保留原始版权声明。

---

## 🙏 致谢

感谢 Lance 和 Williams 的开创性工作，以及 R 社区对聚类算法的讨论。  
本实现参考了 `fastcluster` 和 `protoclust` 的设计思想。

---

**Enjoy error‑free hierarchical clustering in R!**
