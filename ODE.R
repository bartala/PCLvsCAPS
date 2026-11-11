# Run data_analysis.R first (it creates PCL_CAPS)

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
item_names <- paste0("pcl_q", 1:20)
y   <- PCL_CAPS$caps_ptsd1
s6  <- rowSums(PCL_CAPS[, c("pcl_q10","pcl_q2","pcl_q7","pcl_q15","pcl_q4","pcl_q9")])
s20 <- rowSums(PCL_CAPS[, item_names])

# ODE at a specific threshold:
ode_at(y, s6, t = 10)


# Sweep all thresholds and pick the best:
tbl6  <- ode_sweep(y, s6,  thresholds = 1:16)
tbl20 <- ode_sweep(y, s20, thresholds = 19:38)

print(tbl6[which.max(tbl6$ODE), ])    # best 6-item cutoff
print(tbl20[which.max(tbl20$ODE), ])  # best 20-item cutoff

#--------------- spearman correlation of EPDS vs. PCL-5 6-items -------------------
# Create the data frame of EPDS and PCL6 items
df <- read.table(header = TRUE, text = "EPDS PCL6 ...")

res <- cor.test(df$EPDS, df$PCL6, method = "spearman", exact = FALSE)
print(res)

