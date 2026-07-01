png("permuted_correlation_distribution.png", width = 5, height = 5, units = "in", res = 300)
par(mar = c(3, 2, 2, 2))
x <- seq(-0.6, 0.6, length.out = 500)
y <- dnorm(x, mean = 0, sd = 0.15)
plot(x, y, type = "l", lwd = 2,
     xlab = "", ylab = "",
     axes = FALSE, frame.plot = TRUE)
mtext("Permuted correlation distribution", side = 1, line = 1, cex = 0.9)
# Observed correlation dashed line - in the tail
obs <- 0.4
segments(obs, 0, obs, max(y) * 0.95, lty = 2, lwd = 1.5)
text(obs + 0.05, max(y) * 1.0, "Observed correlation", pos = 2, cex = 0.7)
dev.off()