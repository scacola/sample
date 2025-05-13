qnorm_custom <- function(p) {
  return(qnorm(p))
}

twoSurvSampleSizeNI <- function(accrualTime, followTime, alloc, h1, h2, alpha, beta, margin) {
  totalTime <- accrualTime + followTime
  hr1 <- h2 / h1
  hr0 <- margin
  
  p2 <- alloc / (1 + alloc)
  p1 <- 1 - p2
  
  za <- qnorm_custom(1 - alpha)
  zb <- qnorm_custom(1 - beta)
  
  nk <- 5000
  w <- totalTime / nk
  
  i0 <- (p1 * h1 + p2 * h2) / (p1 + hr0 * p2)^2 * 0.5 * w
  i1 <- (p1 * h1 + p2 * h2) / (p1 + hr1 * p2)^2 * 0.5 * w
  om <- (p1 * h1 + p2 * h2) / ((p1 + hr0 * p2) * (p1 + hr1 * p2)) * 0.5 * w
  d1 <- 0
  d2 <- 0
  
  for (k in 1:(nk - 1)) {
    t <- k * w
    s1 <- exp(-h1 * t)
    s2 <- s1^hr1
    f1 <- h1 * s1
    f2 <- hr1 * h1 * s2
    
    if (t <= followTime) {
      g <- 1.0
      dg <- 0.0
    } else {
      g <- -t / accrualTime + totalTime / accrualTime
      dg <- -1.0 / accrualTime
    }
    
    i0 <- i0 + g * s1 * s2 * (p1 * f1 + p2 * f2) / (p1 * s1 + hr0 * p2 * s2)^2 * w
    i1 <- i1 + g * s1 * s2 * (p1 * f1 + p2 * f2) / (p1 * s1 + hr1 * p2 * s2)^2 * w
    om <- om + g * s1 * s2 * (p1 * f1 + p2 * f2) / ((p1 * s1 + hr0 * p2 * s2) * (p1 * s1 + hr1 * p2 * s2)) * w
    d1 <- d1 + dg * s1 * w
    d2 <- d2 + dg * s2 * w
  }
  
  i0 <- hr0 * p1 * p2 * i0
  i1 <- hr1 * p1 * p2 * i1
  om <- (hr0 - hr1) * p1 * p2 * om
  
  n <- ((sqrt(i0) * za + sqrt(i1) * zb) / om)^2
  sample_std <- ceiling(n * p1)
  sample_test <- ceiling(sample_std * alloc)
  
  list(
    total_n = sample_std + sample_test,
    standard_group_n = sample_std,
    test_group_n = sample_test,
    expected_events_std = round(n * p1 * (1 + d1), 1),
    expected_events_test = round(n * p2 * (1 + d2), 1)
  )
}
