oneSurvSampleSize <- function(
    survTime, p1, p2,accrualTime, followTime, 
     alpha, power, side = c("two.sided", "one.sided"), method = c("arcsin", "log-log", "logit", "log", "log-swog", "identity")
){
  h1 = -log(p1) / survTime
  h2 = -log(p2) / survTime
  
  beta <- 1 - power
  
  
  a <- accrualTime
  f <- followTime
  b <- a + f
  s <- survTime
  nk <- 5000
  result <- numeric(2)
  
  if (side == "one.sided") {
    za <- qnorm(1 - alpha)
  } else {
    za <- qnorm(1 - alpha / 2)
  }
  zb <- qnorm(1 - beta)
  
  if (s < f) {
    avar1 <- exp(h1 * s) - 1
    avar2 <- exp(h2 * s) - 1
  } else {
    w <- (s - f) / nk
    t_seq <- seq(f + w, s - w, length.out = nk - 1)
    avar1 <- exp(h1 * f) / (b - f) * 0.5 + sum(exp(h1 * t_seq) / (b - t_seq)) + exp(h1 * s) / (b - s) * 0.5
    avar2 <- exp(h2 * f) / (b - f) * 0.5 + sum(exp(h2 * t_seq) / (b - t_seq)) + exp(h2 * s) / (b - s) * 0.5
    
    avar1 <- w * avar1 * h1 * a + exp(h1 * f) - 1
    avar2 <- w * avar2 * h2 * a + exp(h2 * f) - 1
  }
  
  method <- tolower(method)
  if (method == "arcsin") {
    ncp <- abs(asin(sqrt(p1)) - asin(sqrt(p2)))
    avar2 <- avar2 * exp(-h2 * s) / (1 - exp(-h2 * s)) * 0.25
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "log-log") {
    ncp <- abs(log(-log(p1)) - log(-log(p2)))
    avar2 <- avar2 / (h2^2 * s^2)
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "logit") {
    ncp <- abs(log(p1 / (1 - p1)) - log(p2 / (1 - p2)))
    avar2 <- avar2 / (1 - exp(-h2 * s))^2
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "log") {
    ncp <- abs(log(p1) - log(p2))
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  } else if (method == "log-swog") {
    ncp <- abs(log(p1) - log(p2))
    n <- ceiling(((sqrt(avar2) * za + sqrt(avar1) * zb) / ncp)^2)
    power <- pnorm(-za * sqrt(avar2) / sqrt(avar1) + ncp * sqrt(n) / sqrt(avar1))
  } else {
    # identity
    ncp <- abs(p1 - p2)
    avar2 <- avar2 * exp(-2 * h2 * s)
    n <- ceiling(avar2 * ((za + zb) / ncp)^2)
    power <- pnorm(-za + ncp * sqrt(n) / sqrt(avar2))
  }
  
  result[1] <- n
  result[2] <- power
  names(result) <- c("SampleSize", "Power")
  return(result)
}

oneSurvSampleSize(12, 0.305, 0.435, 24, 24, 0.05, 0.8, "two.sided", "logit")

