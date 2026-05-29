# 加载手工聚类函数库
source("my_hclust.R")

# 1. 准备数据
x <- c(1, 2, 6, 8, 11)
n <- length(x)
labels <- paste0("G", 1:n, "=", x)

# 2. 构造距离矩阵（一维欧氏距离即绝对值）
D <- as.matrix(dist(x, method = "euclidean"))

# 3. 分别用六种方法聚类
hc_single   <- hclust_manual(D, method = "single",   labels = labels)
hc_complete <- hclust_manual(D, method = "complete", labels = labels)
hc_average  <- hclust_manual(D, method = "average",  labels = labels)
hc_centroid <- hclust_manual(D, method = "centroid", labels = labels)
hc_median   <- hclust_manual(D, method = "median",   labels = labels)
hc_ward     <- hclust_manual(D, method = "ward",     labels = labels)

# 4. 输出合并过程（可选）
all_hc <- list(
  最短距离法 = hc_single,
  最长距离法 = hc_complete,
  类平均法   = hc_average,
  重心法     = hc_centroid,
  中间距离法 = hc_median,
  Ward法     = hc_ward
)

for (name in names(all_hc)) {
  cat("\n==========", name, "==========\n")
  hc <- all_hc[[name]]
  n <- length(hc$labels)
  
  # 初始化类名：前 n 个为原始城市，后面为合并产生的新类
  class_names <- character(2 * n - 1)
  class_names[1:n] <- hc$labels
  
  # 逐步输出合并过程
  for (i in 1:(n - 1)) {
    left_id  <- hc$merge[i, 1]
    right_id <- hc$merge[i, 2]
    left_name  <- class_names[left_id]
    right_name <- class_names[right_id]
    
    # 新类的 ID 一定为 n + i（与聚类函数内部的 new_id 顺序一致）
    new_id <- n + i
    new_name <- paste0("新类", i)
    class_names[new_id] <- new_name
    
    cat(sprintf("步骤%2d: 合并 (%s) 与 (%s) 形成 %s, 距离 = %.1f\n",
                i, left_name, right_name, new_name, hc$height[i]))
  }
  
  # 最终类成员（此时只剩下一个类，直接显示其包含的所有原始城市）
  final_members <- hc$members[[1]]  # 最后一个成员列表即为全部样本
  cat("最终类成员:\n")
  cat("  ", paste(hc$labels[final_members], collapse = ", "), "\n")
}

# 5. 绘制树状图
for (name in names(all_hc)) {
  plot_dendrogram(all_hc[[name]], main = name)
  cat("\n当前显示：", name, "\n")
  if (which(names(all_hc) == name) < length(all_hc)) {
    readline("按 Enter 键查看下一张图...")
  }
}