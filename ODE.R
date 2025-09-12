# y: 0/1 true labels; s: numeric score (6-item sum)
ode_at <- function(y, s, t) {
  mean((s >= t) == y)
}

# Compute ODE for a vector of thresholds (defaults to all unique scores)
ode_sweep <- function(y, s, thresholds) {
  acc <- sapply(thresholds, function(t) ode_at(y, s, t))
  data.frame(threshold = thresholds, ODE = acc)
}

# Example usage:
y <- PCL_CAPS_ALL$caps_ptsd1
s6 <- rowSums(PCL_CAPS[, c("pcl_q10","pcl_q2","pcl_q7","pcl_q15","pcl_q4","pcl_q9")])
s20<-rowSums(PCL_CAPS[, c(3:23)])

# ODE at a specific threshold:
ode_at(y, s6, t = 10)

# Sweep all thresholds and pick the best:
tbl <- ode_sweep(y, s6, thresholds = c(1:16) )

tbl <- ode_sweep(y, s20, thresholds = c(19:38) )


#--------------- spearman correlation of EPDS vs. PCL6 -------------------
# Create the data frame of EPDS and PCL6 items
df <- read.table(header = TRUE, text = "EPDS PCL6 ...")

res <- cor.test(df$EPDS, df$PCL6, method = "spearman", exact = FALSE)
print(res)

