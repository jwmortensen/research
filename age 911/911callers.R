library(ggplot2)
load("AgePP.RData")

N <- length(age911)

# Generate plot of posterior delta density
qplot(c(1000, 2000), 
      stat="function", 
      fun=dgamma, 
      geom="line", 
      args=list(shape=N+5, rate=1.2),
      main=bquote("Posterior Density for "~delta),
      ylab="Density",
      xlab=bquote(delta)) +
      geom_vline(aes(xintercept=(N+5)/(1.2))) # This is the MLE for delta
ggsave("delta.pdf")

# Set up variables in order to plot probabilites for the different ages
N_i <- rbind(data.frame(table(age911)),
             data.frame(age911=c("7","101","104","106","107","108","109","110"), Freq=rep(0,8)))
N_i$age911 <- as.numeric(as.character(N_i$age911))
N_i <- N_i[order(N_i$age911),]
alpha <- 2

# Calculate expected values and 95% CI for each probability
p_i <- data.frame(x=0:110, posterior_mean=(N_i$Freq+alpha)/(110*alpha+N))
p_i$lower.bound <- qbeta(0.025, N_i$Freq+alpha, 110*alpha+N-N_i$Freq+alpha)
p_i$upper.bound <- qbeta(0.975, N_i$Freq+alpha, 110*alpha+N-N_i$Freq+alpha)


ggplot(p_i, aes(x=x, y=posterior_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower.bound, ymax=upper.bound)) +
  ggtitle("Expected Value for Dirichlet p's") +
  xlab("p") +
  ylab("Posterior Mean")
ggsave("dirichlet.pdf")
