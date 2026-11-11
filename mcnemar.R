library(stats)
 y_true = PCL_CAPS_ALL_xlsx_BIRTH$caps_ptsd1
 
item_names <- paste0("pcl_q", 1:20)
total20 <- rowSums(PCL_CAPS_ALL_xlsx_BIRTH[, item_names], na.rm = TRUE)
total20 <- factor(total20 >= 23, levels = c(FALSE, TRUE))

total6 <- rowSums(PCL_CAPS_ALL_xlsx_BIRTH[, c('pcl_q10','pcl_q2','pcl_q7','pcl_q15','pcl_q4','pcl_q9')], na.rm = TRUE)
total6 <- factor(total6 >= 7, levels = c(FALSE, TRUE))


tab <- table(total20, total6)
mcnemar.test(tab)  # overall comparison


tab_sens <- table(total20[y_true==1], total6[y_true==1])
mcnemar.test(tab_sens)

tab_spec <- table(total20[y_true==0], total6[y_true==0])
mcnemar.test(tab_spec)
