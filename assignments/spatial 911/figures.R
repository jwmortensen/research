library(ggplot2)
source("spatial911.R")

ndraws <- 1000

draws <- mh.gibbs(ndraws, 1, 0.01, 0.01)
burnin <- 20000
lambda.post.mn <- apply(draws$lambda[-c(1:burnin),], 2, mean)
plot.grid <- matrix(NA, nrow=125, ncol=125)
plot.grid[kp.gp] <- lambda.post.mn

pdf("posteriorLocationDensity.pdf")
image.plot(x.grid, y.grid, plot.grid, axes=FALSE, frame.plot=TRUE, ann=FALSE,
  xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15)) 
plot(houston, add=TRUE)
dev.off()

delta.density <- data.frame(x=density(draws$delta[-c(1:burnin)])$x, y=density(draws$delta[-c(1:burnin)])$y)
q.delta <- quantile(draws$delta[-c(1:burnin)], c(0.025, 0.975))
middle95 <- subset(delta.density, x > q.delta[1] & x < q.delta[2])

ggplot(delta.density, aes(x=x, y=y)) +
         geom_line() +
         labs(x=bquote(delta), y="Density") +
         geom_ribbon(data=middle95, aes(x=x, ymax=y), ymin=0, fill="blue", alpha=0.5) +
         geom_vline(xintercept=mean(draws$delta[-c(1:burnin)]))
ggsave("delta_density.pdf")


s2.density <- data.frame(x=density(draws$sigma2[-c(1:burnin)])$x, y=density(draws$sigma2[-c(1:burnin)])$y)
q.s2 <- quantile(draws$sigma2[-c(1:burnin)], c(0.025, 0.975))
middle95 <- subset(s2.density, x > q.s2[1] & x < q.s2[2])

ggplot(s2.density, aes(x=x, y=y)) +
         geom_line() +
         labs(x=bquote(sigma^2), y="Density") +
         geom_ribbon(data=middle95, aes(x=x, ymax=y), ymin=0, fill="blue", alpha=0.5) +
         geom_vline(xintercept=mean(draws$sigma2[-c(1:burnin)]))
ggsave("s2_density.pdf")

pdf("delta_trace.pdf")
plot(1:ndraws, draws$delta, type="l", xlab="Iteration", ylab=bquote(delta), main=bquote("Trace for "~delta))
dev.off()

pdf("sigma_trace.pdf")
plot(1:ndraws, draws$sigma2, type="l", xlab="Iteration", ylab=bquote(sigma^2), main=bquote("Trace for"~sigma^2))
dev.off()

pdf("lambda_trace.pdf")
plot(1:ndraws, draws$lambda[, 707], type="l", xlab="Iteration", ylab=bquote(lambda[707]), main=bquote("Trace for"~lambda[707]))
dev.off()
