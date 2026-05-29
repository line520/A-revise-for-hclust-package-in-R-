# ============================================================
# 层次聚类核心函数（采用Lance-Williams 通用递推实现）
# 仅依赖距离矩阵，无需原始坐标
# 方法：single, complete, average, centroid, median, ward
# 
# 各个参数含义：
# single      最短距离法
# complete    最长距离法
# average     类平均法
# centroid    重心法
# median      中间距离法
# ward        离差平方和法
# ============================================================

# 工具函数：查找距离矩阵中最小非对角元素
find_min_dist <- function(D, active) {
  m <- length(active)
  min_val <- Inf
  min_i <- min_j <- 0
  for (ii in 1:(m - 1)) {
    for (jj in (ii + 1):m) {
      if (D[ii, jj] < min_val) {
        min_val <- D[ii, jj]
        min_i <- ii
        min_j <- jj
      }
    }
  }
  c(min_i, min_j, min_val)
}

# ============================================================
# 通用层次聚类
# D         : 距离矩阵（若 method 为 ward/centroid/median，自动内部平方）
# method    : 方法名
# labels    : 样品标签（字符向量）
# 返回列表  : merge, height, labels, method, members
hclust_manual <- function(D, method, labels = NULL) {
  n <- nrow(D)
  if (is.null(labels)) labels <- as.character(1:n)
  
  # 对需要平方距离的方法，内部转换为平方距离
  if (method %in% c("ward", "centroid", "median")) {
    D <- D^2          # 后续所有递推都在平方距离上进行
    squared <- TRUE
  } else {
    squared <- FALSE
  }
  
  active <- 1:n
  sizes <- rep(1, n)         # 类大小
  Dmat <- D
  merge <- matrix(0, n-1, 2)
  height <- numeric(n-1)
  members <- as.list(1:n)
  
  for (k in 1:(n-1)) {
    idx <- find_min_dist(Dmat, active)
    p <- idx[1]; q <- idx[2]; d_pq <- idx[3]
    
    real_p <- active[p]
    real_q <- active[q]
    merge[k, ] <- c(real_p, real_q)
    height[k] <- d_pq            # 此时是平方距离（若 squared=TRUE）或原始距离
    
    # 新类成员与大小
    new_members <- c(members[[p]], members[[q]])
    size_p <- sizes[p]; size_q <- sizes[q]
    size_new <- size_p + size_q
    
    # ----- Lance-Williams 递推 -----
    other <- setdiff(1:length(active), c(p, q))
    new_dist <- numeric(length(other))
    
    for (i in seq_along(other)) {
      o <- other[i]
      d_po <- Dmat[p, o]
      d_qo <- Dmat[q, o]
      
      if (method == "single") {
        new_dist[i] <- min(d_po, d_qo)
      } else if (method == "complete") {
        new_dist[i] <- max(d_po, d_qo)
      } else if (method == "average") {
        new_dist[i] <- (size_p * d_po + size_q * d_qo) / size_new
      } else if (method == "centroid") {
        # Lance-Williams 系数：重心法（平方距离）
        alpha_p <- size_p / size_new
        alpha_q <- size_q / size_new
        beta    <- - size_p * size_q / (size_new^2)
        new_dist[i] <- alpha_p * d_po + alpha_q * d_qo + beta * d_pq
      } else if (method == "median") {
        # 中间距离法（平方距离）
        alpha_p <- 0.5
        alpha_q <- 0.5
        beta    <- -0.25
        new_dist[i] <- alpha_p * d_po + alpha_q * d_qo + beta * d_pq
      } else if (method == "ward") {
        # Ward 法（平方距离，使用 Lance-Williams 系数）
        denom <- size_p + size_q + sizes[o]
        alpha_p <- (size_p + sizes[o]) / denom
        alpha_q <- (size_q + sizes[o]) / denom
        beta    <- - sizes[o] / denom
        new_dist[i] <- alpha_p * d_po + alpha_q * d_qo + beta * d_pq
      }
    }
    
    # 更新 active, sizes, members, Dmat
    new_active <- active[-c(p, q)]
    new_sizes <- sizes[-c(p, q)]
    new_members_list <- members[-c(p, q)]
    new_id <- if (length(new_active) > 0) max(active) + 1 else max(active) + 1
    new_active <- c(new_active, new_id)
    new_sizes <- c(new_sizes, size_new)
    new_members_list <- c(new_members_list, list(new_members))
    
    m <- length(new_active)
    Dnew <- matrix(0, m, m)
    for (i in seq_along(other)) {
      for (j in seq_along(other)) {
        Dnew[i, j] <- Dmat[other[i], other[j]]
      }
    }
    for (i in seq_along(other)) {
      Dnew[i, m] <- new_dist[i]
      Dnew[m, i] <- new_dist[i]
    }
    
    active <- new_active
    sizes <- new_sizes
    Dmat <- Dnew
    members <- new_members_list
  }
  
  # 若使用过平方距离，将高度转换回原始距离（便于绘图和解释）
  if (squared) {
    height <- sqrt(height)
  }
  
  list(merge = merge, height = height, labels = labels,
       method = method, members = members)
}


# ============================================================
# 树状图绘制函数
plot_dendrogram <- function(hc, main = "", cex.label = 0.8) {
  merge <- hc$merge
  height <- hc$height
  n <- nrow(merge) + 1
  labels <- hc$labels
  
  # 计算所有节点的 y 坐标（叶节点按顺序 1:n）
  y_pos <- numeric(2 * n - 1)
  x_pos <- numeric(2 * n - 1)
  y_pos[1:n] <- 1:n
  x_pos[1:n] <- 0
  
  # 按高度升序处理（若 height 非单调则先排序）
  ord <- order(height)
  for (i in ord) {
    left <- merge[i, 1]
    right <- merge[i, 2]
    new_node <- n + i
    y_left <- if (left <= n) y_pos[left] else y_pos[left]
    y_right <- if (right <= n) y_pos[right] else y_pos[right]
    y_pos[new_node] <- (y_left + y_right) / 2
    x_pos[new_node] <- height[i]
  }
  
  # 收集所有线段端点需要的 x 坐标范围
  all_x <- c(0, height)
  for (i in 1:(n-1)) {
    left <- merge[i, 1]; right <- merge[i, 2]
    all_x <- c(all_x, x_pos[left], x_pos[right], height[i])
  }
  x_min <- min(0, min(all_x, na.rm = TRUE))
  x_max <- max(all_x, na.rm = TRUE)
  x_range <- c(x_min - 0.05 * (x_max - x_min), x_max + 0.05 * (x_max - x_min))
  
  # 自动调整左边距以容纳长标签
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  max_str_width <- max(strwidth(labels, units = "inches", cex = cex.label))
  if (max_str_width > 0.5) {
    mar <- par("mar")
    mar[2] <- max(mar[2], max_str_width * 2)  # 增大左边距
    par(mar = mar)
  }
  
  # 绘图
  plot(0, 0, type = "n",
       xlim = x_range,
       ylim = c(0.5, n + 0.5),
       xlab = "距离", ylab = "", main = main, yaxt = "n")
  axis(2, at = 1:n, labels = labels, las = 1, cex.axis = cex.label)
  
  # 绘制树枝
  for (i in 1:(n - 1)) {
    left <- merge[i, 1]; right <- merge[i, 2]
    new_node <- n + i
    h <- height[i]
    y_left <- y_pos[left]
    y_right <- y_pos[right]
    y_new <- y_pos[new_node]
    x_left <- x_pos[left]
    x_right <- x_pos[right]
    
    # 水平连接线
    segments(h, y_left, h, y_right, lwd = 1.5)
    # 从子节点到水平线的垂直线
    segments(x_left, y_left, h, y_left, lwd = 1.5)
    segments(x_right, y_right, h, y_right, lwd = 1.5)
  }
}