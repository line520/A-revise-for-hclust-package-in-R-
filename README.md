# A-revise-for-hclust-package-in-R-
# 手工层次聚类： `hclust` 的可靠替代

本项目基于 **Lance–Williams 通用递推公式**，从零实现了层次聚类的六种常用方法（最短距离、最长距离、类平均、重心法、中间距离法、Ward 法）。  
与 R 自带的 `stats::hclust` 相比，本实现更加透明、可控，且修正了部分方法（如重心法和 Ward 法）在原始包中可能存在的计算偏差问题。

> **背景**：R 语言（≤ 4.6.0 及更早版本）中的 `hclust` 函数在实现某些方法（尤其是 `centroid` 和 `median`）时，内部存在未知问题，导致聚类结果与教科书公式不符。下面将在例子1展示错误和正确输出

---


## 📐 数学原理

### 1. Lance–Williams 递推公式

设类 A 和 B 合并为新类 C = A ∪ B，|A|, |B| 为类的大小（样本数）。对于任意另一个类 X，新类 C 与 X 之间的距离 d(C, X) 可以由 d(A, X)、d(B, X) 以及 d(A, B) 递推得到：

$$
d(C, X) = \alpha_A d(A, X) + \alpha_B d(B, X) + \beta d(A, B) + \gamma |d(A, X) - d(B, X)|
$$

对于**非加权**的六种方法，γ = 0，系数 (α_A, α_B, β) 如下表所示：

| 方法        | α_A                         | α_B                         | β                             | 
|-------------|-----------------------------|-----------------------------|-------------------------------|
| single      | 1/2                         | 1/2                         | 0                             | 
| complete    | 1/2                         | 1/2                         | 0                             |
| average     | |A|/(|A|+|B|)              | |B|/(|A|+|B|)               | 0                             | 
| centroid    | |A|/(|A|+|B|)               | |B|/(|A|+|B|)               | - (|A|*|B|)/(|A|+|B|)^2        | 
| median      | 1/2                         | 1/2                         | -1/4                          |
| ward        | (|A|+|X|)/(|A|+|B|+|X|)     | (|B|+|X|)/(|A|+|B|+|X|)     | - |X|/(|A|+|B|+|X|)            |

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

### 函数说明

#### `hclust_manual(D, method, labels = NULL)`

| 参数     | 类型        | 说明                                                                 |
|----------|-------------|----------------------------------------------------------------------|
| `D`      | 距离矩阵    | $n \times n$ 对称矩阵，非负，对角为 0。                              |
| `method` | 字符串      | 可选：`"single"`, `"complete"`, `"average"`, `"centroid"`, `"median"`, `"ward"`。|
| `labels` | 字符向量    | 样本标签，长度 = $n$。默认为 `as.character(1:n)`。                  |

**返回值**：一个列表，包含以下组件：
- `merge`：$ (n-1) \times 2$ 矩阵，每行表示合并的两个类编号（负数表示原始样本）。
- `height`：长度为 $n-1$ 的向量，每次合并时的距离（原始欧氏距离）。
- `labels`：样本标签。
- `method`：使用的聚类方法。single, complete, average, centroid, median, ward分别对应最短距离法、最长距离法、类平均法、重心法、中间距离法、离差平方和法（ward法）
- `members`：列表，记录每个（中间）类包含的原始样本索引。

#### `plot_dendrogram(hc, main = "", cex.label = 0.8)`

- `hc`：由 `hclust_manual` 返回的对象。
- `main`：图形标题。
- `cex.label`：Y 轴标签字体大小。

---

## 📊 示例演示

### 示例 ：一维数据（ex6.1.R）

5 个点：$x = [1, 2, 6, 8, 11]$。  

```r
# 1. 加载函数
source("my_hclust.R")

# 2. 准备距离矩阵（例如 5 个点的一维坐标）
x <- c(1, 2, 6, 8, 11)
D <- as.matrix(dist(x))

# 3. 进行聚类（以 Ward 法为例）
hc <- hclust_manual(D, method = "centroid", labels = paste0("G", 1:5))

# 4. 查看合并过程
print(hc$merge)
print(hc$height)

# 5. 绘制树状图
plot_dendrogram(hc, main = "Centroid 聚类")
```

## ⚠️ 原hclust的错误点

例如在上述例子中正常树状图输出结果为：
<img width="1082" height="1164" alt="image" src="https://github.com/user-attachments/assets/b5051763-df56-4389-b2e1-d0306bc81b36" />
在原来的hclust包中结果为
<img width="1196" height="1282" alt="image" src="https://github.com/user-attachments/assets/96a26e14-9cf9-4e5d-bdf2-e3d81f5d507a" />
可见在进行最后一项合并的时候，距离计算出现了错误，正常应为4，实际计算出来为3.5，原hclust包运行代码如下
```r
x <- c(1,2,6,8,11)
d <- dist(x)
par(mfrow = c(1, 2))
plot(hclust(d, "centroid"), main = "重心法")

```
---

## 📄 许可证

本项目采用 [MIT 许可证](LICENSE)，可自由使用、修改和分发。

---

**Enjoy error‑free hierarchical clustering in R!**
