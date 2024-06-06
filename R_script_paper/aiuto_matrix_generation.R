# help build matrix by hand
a_beta <- 1
b_beta <- 10

H <- 2000
bins_ref <- c(0, 0.0002, 0.0007, 0.002, 0.007, 0.02, 0.07, 0.2, 0.5, 1)
#bins <- c(0, 0.0002, 0.002, 0.02, 0.07, 0.2, 0.5, 1)
counts <- matrix(NA, nrow = 100, length(bins_ref) -1)

for (iter in 1:100){
  pis <- rbeta(H, shape1 = a_beta, shape2 = b_beta)
  
  categorized_numbers <- cut(pis, breaks = bins_ref, right = TRUE)
  counts[iter, ] <- table(categorized_numbers)
}

print(colMeans(counts))

