elements <- c(0, 0.6, 0.32, 0, 0.14, 0.4, 0, 0.68, 0, 0.3, 0.6, 0.35, 0, 0.12, 0, 0, 0, 0 ,0, 0.56, 0, 0.05, 0, 0.88, 0)
P <- matrix(elements, nrow = 5, ncol = 5)

#verif
rowSums(P)

P2 = P%*%P
P2

P3 = P2%*%P
P3

P4 = P3%*%P
P4

P5 = P4%*%P
P5