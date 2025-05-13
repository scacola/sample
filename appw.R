library(shiny)
library(DT)
# Lakatos sample size function for two-sample survival test


twoSurvSampleSizeNI <- function(accrualTime, followTime, alloc, h1, h2, alpha, beta, margin) {
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
    total_n = sample_std + sample_test,
    standard_group_n = sample_std,
    test_group_n = sample_test,
    expected_events_std = round(n * p1 * (1 + d1), 1),
    expected_events_test = round(n * p2 * (1 + d2), 1),
    total_expected_events = round(n * p1 * (1 + d1), 1) + round(n * p2 * (1 + d2), 1)
  )
  
}
lakatosSampleSize <- function(
    accrualTime, followTime, alloc,
    h1, h2, alpha, power,
    method = c("logrank", "gehan", "tarone-ware"),
    side = c("two.sided", "one.sided"),
    b = 24
) {
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
      standard_group_n = std_n,
      test_group_n = test_n,
      total_n = total_n,
      actual_power = round(pow, 4),
      expected_events_std = round(n * eEvt1, 1),
      expected_events_test = round(n * eEvt2, 1),
      total_expected_events = round(n * (eEvt1 + eEvt2), 1)
    )
  } else {
    result <- list(
      error = "Non-finite or zero effect size detected. Unable to compute sample size."
    )
  }
  
  return(result)
}


ui <- fluidPage(
  titlePanel("Two-sample Survival Sample Size Calculator"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("test_type", "Test Type:",
                   choices = c("Non-Inferiority" = "ni", "Inferiority" = "sup"),
                   selected = "ni"),
      numericInput("syear", "Survival Time (years):", value = 12),
      numericInput("yrsurv1", "Survival Probability (Standard Group):", value = 0.5),
      numericInput("yrsurv2", "Survival Probability (Test Group):", value = 0.3),
      numericInput("alloc", "Allocation Ratio (Test / Standard):", value = 1),
      numericInput("accrual", "Accrual Time (months):", value = 24),
      numericInput("follow", "Follow-up Time (months):", value = 24),
      numericInput("alpha", "Significance Level (alpha):", value = 0.025),
      numericInput("power", "Power (1 - Beta):", value = 0.8),
      conditionalPanel(
        condition = "input.test_type == 'ni'",
        numericInput("margin", "Non-inferiority Margin (HR):", value = 1.3)
      ),
      conditionalPanel(
        condition = "input.test_type == 'sup'",
        selectInput("method", "Test Method:", choices = c("logrank", "gehan", "tarone-ware")),
        selectInput("side", "Test Direction:", choices = c("two.sided", "one.sided"))
      ),
      actionButton("calc", "Calculate")
    ),
    mainPanel(
      DTOutput("result_table"),
      tags$hr(),
      tags$div("Reference: Jung SH, Chow SC. On sample size calculation for comparing survival curves under general hypothesis testing. Journal of Biopharmaceutical Statistics 2012; 22(3):485–495."),
      tags$hr(),
      tags$div("Reference: Lakatos E. Sample sizes based on the log-rank statistic in complex clinical trials. Biometrics 1988;44:229–241."),
      tags$hr(),
      tags$div("Reference: Lakatos E, Lan KK. A comparison of sample size methods for the logrank statistic. Statistics in Medicine 1992;11(2):179–191.")
          )
  )
)

server <- function(input, output) {
  observeEvent(input$calc, {
    h1 <- -log(input$yrsurv1) / input$syear
    h2 <- -log(input$yrsurv2) / input$syear
    beta <- 1 - input$power
    
    res_df <- NULL
    
    if (input$test_type == "ni") {
      res <- tryCatch({
        twoSurvSampleSizeNI(
          accrualTime = input$accrual,
          followTime = input$follow,
          alloc = input$alloc,
          h1 = h1,
          h2 = h2,
          alpha = input$alpha,
          beta = beta,
          margin = input$margin
        )
      }, error = function(e) NULL)
      
      if (!is.null(res) && all(c("total_n", "standard_group_n", "test_group_n", "expected_events_std", "expected_events_test", "total_expected_events") %in% names(res))) {
        res_df <- data.frame(
          Metric = c(
            "Total Sample Size",
            "Standard Group Sample Size",
            "Test Group Sample Size",
            "Expected Events (Standard Group)",
            "Expected Events (Test Group)",
            "Total Expected Events"
          ),
          Value = c(
            res$total_n,
            res$standard_group_n,
            res$test_group_n,
            res$expected_events_std,
            res$expected_events_test,
            res$total_expected_events
          )
        )
      } else {
        res_df <- data.frame(Metric = "Error", Value = "Invalid or missing output from NI function")
      }
      
    } else {
      res <- lakatosSampleSize(
        accrualTime = input$accrual,
        followTime = input$follow,
        alloc = input$alloc,
        h1 = h1,
        h2 = h2,
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
    }
    
    output$result_table <- renderDT({
      datatable(
        res_df,
        options = list(
          dom = 't',
          ordering = FALSE,
          columnDefs = list(list(className = 'dt-center', targets = '_all'))
        ),
        rownames = FALSE
      )
    })
  })
}

shinyApp(ui, server)
