library(shiny)
library(DT)
# Lakatos sample size function for two-sample survival test


rcode1 <- 'twoSurvSampleSizeNI <- function(syear, yrsurv1, yrsurv2, alloc, accrualTime, followTime, alpha, power, margin) {
  
  h1 <- -log(yrsurv1) / syear
  h2 <- -log(yrsurv2) / syear
  beta <- 1 - power
  
  totalTime <- accrualTime + followTime
  hr1 <- h2 / h1
  hr0 <- margin
  
  p2 <- alloc / (1 + alloc)
  p1 <- 1 - p2
  
  za <- qnorm(1 - alpha)
  zb <- qnorm(1 - beta)
  
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
    Sample_size_of_standard_group = sample_std,
    Sample_size_of_test_group = sample_test,
    Total_sample_size = sample_std + sample_test,
    Expected_event_numbers_of_standard_group = round(n * p1 * (1 + d1), 1),
    Expected_event_numbers_of_test_group = round(n * p2 * (1 + d2), 1),
    Total_expected_event_numbers = round(n * p1 * (1 + d1), 1) + round(n * p2 * (1 + d2), 1)
  )
  
}'

rcode2 <- '
lakatosSampleSize <- function(syear, yrsurv1, yrsurv2, alloc,
    accrualTime, followTime, 
    alpha, power,
    method = c("logrank", "gehan", "tarone-ware"),
    side = c("two.sided", "one.sided"),
    b = 24
) {
  
  h1 <- -log(yrsurv1) / syear
  h2 <- -log(yrsurv2) / syear
  beta <- 1 - power
  
  
  method <- match.arg(method)
  side <- match.arg(side)
  totalTime <- accrualTime + followTime
  m <- floor(totalTime * b)
  
  ti <- seq(0, totalTime, length.out = m)
  n1 <- numeric(m)
  n2 <- numeric(m)
  
  allocRatio <- 1 / (1 + alloc)
  hr <- h2 / h1
  numer <- denom <- phi <- di <- wi <- e <- n <- 0
  eEvt1 <- eEvt2 <- 0
  
  #i in seq_along(ti)
  for (i in 1:m) {
    if (i == 1) {
      n1[i] <- allocRatio
      n2[i] <- 1 - allocRatio
    } else {
      if (ti[i - 1] <= followTime) {
        n1[i] <- n1[i - 1] * (1 - h1 / b)
        n2[i] <- n2[i - 1] * (1 - h2 / b)
      } else {
        n1[i] <- n1[i - 1] * (1 - h1 / b - 1 / (b * (totalTime - ti[i - 1])))
        n2[i] <- n2[i - 1] * (1 - h2 / b - 1 / (b * (totalTime - ti[i - 1])))
      }
    }
    
    phi <- n2[i] / n1[i]
    di <- (n1[i] * h1 + n2[i] * h2) / b
  
    # if(i<=10 | i>=1045) print(c(i, n1[i]+n2[i]))
    
    # if(1){
    #   print(i)
    #   print(n1[i]+n2[i])
    # }
    
    wi <- switch(method,
                 logrank = 1,
                 gehan = n1[i] + n2[i],
                 `tarone-ware` = sqrt(max( n1[i] + n2[i], 0)))
    
    numer <- numer + di * wi * ((hr * phi) / (1 + hr * phi) - phi / (1 + phi))
    denom <- denom + di * wi^2 * (phi / (1 + phi)^2)
    eEvt1 <- eEvt1 + n1[i] * h1 / b
    eEvt2 <- eEvt2 + n2[i] * h2 / b
  }
  
  e <- numer / sqrt(denom)
  result <- list()
  
  if (is.finite(e) && abs(e) > 0) {
    target_beta <- 1 - power
    
    repeat {
      n <- n + 1
      if (side == "two.sided") {
        pow <- pnorm(-sqrt(n) * e - qnorm(1 - alpha / 2)) +
          pnorm(sqrt(n) * e - qnorm(1 - alpha / 2))
      } else {
        pow <- pnorm(-sqrt(n) * e - qnorm(1 - alpha))
      }
      if (pow >= power) break
    }
    
    std_n <- ceiling(n * allocRatio)
    test_n <- ceiling(std_n * alloc)
    total_n <- std_n + test_n
    
    result <- list(
      Sample_size_of_standard_group = std_n,
      Sample_size_of_test_group = test_n,
      Total_sample_size = total_n,
      Expected_event_numbers_of_standard_group = round(n * eEvt1, 1),
      Expected_event_numbers_of_test_group = round(n * eEvt2, 1),
      Total_expected_event_numbers = round(n * (eEvt1 + eEvt2), 1),
      Actual_power = round(pow, 4)
    )
  } else {
    result <- list(
      error = "Non-finite or zero effect size detected. Unable to compute sample size."
    )
  }
  
  return(result)
}
'

rcode3 <- '
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
}'


twoSurvSampleSizeNI <- function(syear, yrsurv1, yrsurv2, alloc, accrualTime, followTime,  alpha, power, margin) {
  
  h1 <- -log(yrsurv1) / syear
  h2 <- -log(yrsurv2) / syear
  beta <- 1 - power
  
  totalTime <- accrualTime + followTime
  hr1 <- h2 / h1
  hr0 <- margin
  
  p2 <- alloc / (1 + alloc)
  p1 <- 1 - p2
  
  za <- qnorm(1 - alpha)
  zb <- qnorm(1 - beta)
  
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
    Sample_size_of_standard_group = sample_std,
    Sample_size_of_test_group = sample_test,
    Total_sample_size = sample_std + sample_test,
    Expected_event_numbers_of_standard_group = round(n * p1 * (1 + d1), 1),
    Expected_event_numbers_of_test_group = round(n * p2 * (1 + d2), 1),
    Total_expected_event_numbers = round(n * p1 * (1 + d1), 1) + round(n * p2 * (1 + d2), 1)
  )
  
}

lakatosSampleSize <- function(syear, yrsurv1, yrsurv2, alloc,
    accrualTime, followTime,
    alpha, power,
    method = c("logrank", "gehan", "tarone-ware"),
    side = c("two.sided", "one.sided"),
    b = 24
) {
  
  h1 <- -log(yrsurv1) / syear
  h2 <- -log(yrsurv2) / syear
  beta <- 1 - power
  
  
  method <- match.arg(method)
  side <- match.arg(side)
  totalTime <- accrualTime + followTime
  m <- floor(totalTime * b)
  
  ti <- seq(0, totalTime, length.out = m)
  n1 <- numeric(m)
  n2 <- numeric(m)
  
  allocRatio <- 1 / (1 + alloc)
  hr <- h2 / h1
  numer <- denom <- phi <- di <- wi <- e <- n <- 0
  eEvt1 <- eEvt2 <- 0
  
  #i in seq_along(ti)
  for (i in 1:m) {
    if (i == 1) {
      n1[i] <- allocRatio
      n2[i] <- 1 - allocRatio
    } else {
      if (ti[i - 1] <= followTime) {
        n1[i] <- n1[i - 1] * (1 - h1 / b)
        n2[i] <- n2[i - 1] * (1 - h2 / b)
      } else {
        n1[i] <- n1[i - 1] * (1 - h1 / b - 1 / (b * (totalTime - ti[i - 1])))
        n2[i] <- n2[i - 1] * (1 - h2 / b - 1 / (b * (totalTime - ti[i - 1])))
      }
    }
    
    phi <- n2[i] / n1[i]
    di <- (n1[i] * h1 + n2[i] * h2) / b
  
    # if(i<=10 | i>=1045) print(c(i, n1[i]+n2[i]))
    
    # if(1){
    #   print(i)
    #   print(n1[i]+n2[i])
    # }
    
    wi <- switch(method,
                 logrank = 1,
                 gehan = n1[i] + n2[i],
                 `tarone-ware` = sqrt(max( n1[i] + n2[i], 0)))
    
    numer <- numer + di * wi * ((hr * phi) / (1 + hr * phi) - phi / (1 + phi))
    denom <- denom + di * wi^2 * (phi / (1 + phi)^2)
    eEvt1 <- eEvt1 + n1[i] * h1 / b
    eEvt2 <- eEvt2 + n2[i] * h2 / b
  }
  
  e <- numer / sqrt(denom)
  result <- list()
  
  if (is.finite(e) && abs(e) > 0) {
    target_beta <- 1 - power
    
    repeat {
      n <- n + 1
      if (side == "two.sided") {
        pow <- pnorm(-sqrt(n) * e - qnorm(1 - alpha / 2)) +
          pnorm(sqrt(n) * e - qnorm(1 - alpha / 2))
      } else {
        pow <- pnorm(-sqrt(n) * e - qnorm(1 - alpha))
      }
      if (pow >= power) break
    }
    
    std_n <- ceiling(n * allocRatio)
    test_n <- ceiling(std_n * alloc)
    total_n <- std_n + test_n
    
    result <- list(
      Sample_size_of_standard_group = std_n,
      Sample_size_of_test_group = test_n,
      Total_sample_size = total_n,
      Expected_event_numbers_of_standard_group = round(n * eEvt1, 1),
      Expected_event_numbers_of_test_group = round(n * eEvt2, 1),
      Total_expected_event_numbers = round(n * (eEvt1 + eEvt2), 1),
      Actual_power = round(pow, 4)
    )
  } else {
    result <- list(
      error = "Non-finite or zero effect size detected. Unable to compute sample size."
    )
  }
  
  return(result)
}


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



ui <- fluidPage(
  includeCSS("www/style.css"),  
  titlePanel("Two-sample Survival Sample Size Calculator"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("test_type", "Test Type:",
                   choices = c("Non-Inferiority" = "ni", "Superiority" = "sup", "One-sample" = "one"),
                   selected = "ni", inline = T),
      numericInput("syear", "Survival Time :", value = 12),
      uiOutput("surv_ui"),
      uiOutput("alloc_ui"),
      fluidRow(
        column(6,
               numericInput("accrual", "Accrual Time :", value = 24)
        ),
        column(6,
               numericInput("follow", "Follow-up Time :", value = 24)
        )
      ),
      uiOutput("alpha_ui"),
      numericInput("power", "Power (1 - Beta):", value = 0.8),
      conditionalPanel(
        condition = "input.test_type == 'ni'",
        numericInput("margin", "Non-inferiority Margin :", value = 1.3)
      ),
      conditionalPanel(
        condition = "input.test_type == 'sup'",
        selectInput("method", "Test Method:", choices = c("logrank", "gehan", "tarone-ware")),
        selectInput("side", "Hypothesis:", choices = c("two.sided", "one.sided"))
      ),
      conditionalPanel(
        condition = "input.test_type=='ni'",
        selectInput("side", "Hypothesis:", choices = "one.sided")
      ),
      conditionalPanel(
        condition = "input.test_type=='one'",
        selectInput("side", "Hypothesis:", choices = c("two.sided", "one.sided")),
        selectInput("trans", "Transformation:", choices = c("arcsin", "log-log", "logit", "log", "log-swog", "identity"))
      ),
      actionButton("calc", "Calculate")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("analysis",
                 DTOutput("result_table"),
              tags$h3("Reference"),
                 tags$div("Jung SH, Chow SC. On sample size calculation for comparing survival curves under general hypothesis testing. Journal of Biopharmaceutical Statistics 2012; 22(3):485–495."),
   
                 tags$div(" Lakatos E. Sample sizes based on the log-rank statistic in complex clinical trials. Biometrics 1988; 44:229–241."),
                
                 tags$div(" Lakatos E, Lan KK. A comparison of sample size methods for the logrank statistic. Statistics in Medicine 1992; 11(2):179–191."),
                
                 tags$div(" Fleming TR, Harrington DP. Counting Processes and Survival Analysis. New York: Wiley, 1991, 236–237, Example 6.3.1."),
               
                 tags$div(" Andersen PK, Borgan Ø, Gill RD, Keiding N. Statistical Models Based on Counting Processes. New York: Springer-Verlag, 1993, 176–287, Section IV.1–3."),
                 
                 tags$div(" Bie O, Borgan Ø, Liestøl K. Confidence intervals and confidence bands for the cumulative hazard rate function and their small sample properties. Scandinavian Journal of Statistics 1987; 14(3): 221–233."),
                 
                 tags$div(" Borgan Ø, Liestøl K. A note on confidence intervals and bands for the survival function based on transformations. Scandinavian Journal of Statistics 1990; 17(1): 35–41."),
                 
                 tags$div(" Nagashima K, Noma H, Sato Y, Gosho M. Sample size calculations for single-arm survival studies using transformations of the Kaplan–Meier estimator. Pharmaceutical Statistics 2020. DOI: 10.1002/pst.2090. [arXiv:2012.03355].")
        ),
        tabPanel("R code",
                 uiOutput("codeBlock")
        )
      )
      
  )
)
)

server <- function(input, output) {
  observeEvent(input$calc, {
    
    
    res_df <- NULL
    
    if (input$test_type == "ni") {
      res <- tryCatch({
        twoSurvSampleSizeNI(
          syear = input$syear,
          yrsurv1 = input$yrsurv1,
          yrsurv2 = input$yrsurv2, 
          accrualTime = input$accrual,
          followTime = input$follow,
          alloc = input$alloc,
          alpha = input$alpha,
          power = input$power,
          margin = input$margin
        )
      }, error = function(e) NULL)
      
      if (!is.null(res)) {
        
        res_df <- data.frame(
          Metric = names(res),
          Value = unname(unlist(res))
        )
      } else {
        res_df <- data.frame(Metric = "Error", Value = "Invalid or missing output from NI function")
      }
      
    } else if(input$test_type == "sup") {
      res <- lakatosSampleSize(
        syear = input$syear,
        yrsurv1 = input$yrsurv1,
        yrsurv2 = input$yrsurv2,
        accrualTime = input$accrual,
        followTime = input$follow,
        alloc = input$alloc,
        alpha = input$alpha,
        power = input$power,
        method = input$method,
        side = input$side
      )
      if (is.null(res$error)) {
        res_df <- data.frame(
          Metric = names(res),
          Value = unname(unlist(res))
        )
      } else {
        res_df <- data.frame(Metric = "Error", Value = res$error)
      }
    }else if(input$test_type =="one"){
      res <- tryCatch({oneSurvSampleSize(
        survTime = input$syear,
        p1 = input$p1,
        p2 = input$p2,
        accrualTime = input$accrual,
        followTime = input$follow,
        alpha = input$alpha,
        power = input$power,
        side = input$side,
        method = input$trans
      )}, error = function(e) NULL)
      
      if (!is.null(res)) {
        res_df <- data.frame(
          Metric = names(res),
          Value = unname(unlist(res))
        )
      } else {
        res_df <- data.frame(Metric = "Error", Value = "Invalid or missing output from NI function")
      }
    }
    
    output$result_table <- renderDT({
      datatable(
        res_df,
        colnames = c("Metric", "N"),
        options = list(
          dom = 't',
          ordering = FALSE,
          columnDefs = list(list(className = 'dt-center', targets = '_all'))
        ),
        rownames = FALSE
      )
    })
  })
  
  output$codeBlock <- renderUI({
    if (input$test_type == "ni") {
      tags$pre(tags$code(rcode1))
    } else if(input$test_type == "sup"){
      tags$pre(tags$code(rcode2))
    }else{
      tags$pre(tags$code(rcode3))
    }
  })
    
  output$alpha_ui <- renderUI({
    if (input$test_type == "ni") {
      numericInput("alpha", "Significance Level (alpha):", value = 0.025)
    } else {
      numericInput("alpha", "Significance Level (alpha):", value = 0.05)
    }
  })
  
  output$alloc_ui <- renderUI({
    if (input$test_type %in% c("ni","sup")){
      numericInput("alloc", "Allocation Ratio :", value = 1)
    }
  })
  
  output$surv_ui <- renderUI({
    if (input$test_type %in% c("ni", "sup")){
      fluidRow(
        column(6,
               numericInput("yrsurv1", "Survival Probability (Standard Group):", value = 0.305)
        ),
        column(6,
               numericInput("yrsurv2", "Survival Probability (Test Group):", value = 0.435)
        )
      )
    }else{
      fluidRow(
        column(6,
               numericInput("p1", "Null survival probability:", value = 0.305)
        ),
        column(6,
               numericInput("p2", "Alternative survival probability", value = 0.435)
        )
      )
    }
  })
  
}

shinyApp(ui, server)
