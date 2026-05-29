# 加载升级后的函数
source("my_hclust.R")

# ============== 输入飞行里程距离矩阵 ==============
# 10个城市
cities <- c("亚特兰大","芝加哥","丹佛","休斯顿","洛杉矶",
            "迈阿密","纽约","旧金山","西雅图","华盛顿")

# 距离矩阵（单位：英里）
D <- matrix(c(
  0,  587,1212, 701,1936, 604, 748,2139,2182, 543,
  587,  0, 920, 940,1745,1188, 713,1858,1737, 597,
  1212,920,  0, 879, 831,1726,1631, 949,1021,1494,
  701, 940,879,   0,1374, 968,1420,1645,1891,1220,
  1936,1745,831,1374,   0,2339,2451, 347, 959,2300,
  604,1188,1726,968,2339,   0,1092,2594,2734, 923,
  748, 713,1631,1420,2451,1092,   0,2571,2408, 205,
  2139,1858,949,1645, 347,2594,2571,   0, 678,2442,
  2182,1737,1021,1891, 959,2734,2408, 678,   0,2329,
  543, 597,1494,1220,2300, 923, 205,2442,2329,   0
), nrow = 10, byrow = TRUE)

rownames(D) <- cities
colnames(D) <- cities

# ============== 四种聚类方法 ==============
hc_single   <- hclust_manual(D, method = "single",   labels = cities)
hc_average  <- hclust_manual(D, method = "average",  labels = cities)
hc_centroid <- hclust_manual(D, method = "centroid", labels = cities)
hc_ward     <- hclust_manual(D, method = "ward",     labels = cities)

# ============== 输出合并过程 ==============
all_hc <- list(最短距离法 = hc_single, 
               类平均法 = hc_average,
               重心法 = hc_centroid, 
               Ward法 = hc_ward)

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

# ============== 绘制树状图 ==============
for (name in names(all_hc)) {
  plot_dendrogram(all_hc[[name]], main = name)
  cat("\n当前显示：", name, "\n")
  if (which(names(all_hc) == name) < length(all_hc)) {
    readline("按 Enter 键查看下一张图...")
  }
}