library(ggplot2)
source("spatial911.R")

ndraws <- 1000

draws <- mh.gibbs(ndraws, 1, 0.01, 0.01)
burnin <- 1000
lambda.post.mn <- apply(draws$lambda[-c(1:burnin),], 2, mean)
lambda.post.low <- apply(draws$lambda[-c(1:burnin),], 2, quantile, probs=c(0.025))
lambda.post.high <- apply(draws$lambda[-c(1:burnin),], 2, quantile, probs=c(0.975))
plot.grid <- matrix(NA, nrow=125, ncol=125)
low.grid <- matrix(NA, nrow=125, ncol=125)
high.grid <- matrix(NA, nrow=125, ncol=125)
plot.grid[kp.gp] <- lambda.post.mn
low.grid[kp.gp] <- lambda.post.low
high.grid[kp.gp] <- lambda.post.high

pdf("posteriorLocationDensity.pdf")
image.plot(x.grid, y.grid, plot.grid, axes=FALSE, frame.plot=TRUE, ann=FALSE,
  xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15))#, zlim=c(0, 0.012)) 
dev.off()
pdf("posteriorLow.pdf")
image.plot(x.grid, y.grid, low.grid, axes=FALSE, frame.plot=TRUE, ann=FALSE,
  xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15), zlim=c(0, 0.012))
dev.off()
pdf("posteriorHigh.pdf")
image.plot(x.grid, y.grid, high.grid, axes=FALSE, frame.plot=TRUE, ann=FALSE,
  xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15), zlim=c(0, 0.014))
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

qplot(1:ndraws, draws$delta, geom="line", xlab="Iteration", ylab=bquote(delta))
ggsave("delta_trace.png")

qplot(1:ndraws, draws$lambda[, 707], geom="line", xlab="Iteration", ylab=bquote(lambda[707]))
ggsave("lambda_trace.png")
