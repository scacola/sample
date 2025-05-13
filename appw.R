library(shiny)
# Lakatos sample size function for two-sample survival test
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
  
  for (i in seq_along(ti)) {
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
    wi <- switch(method,
                 logrank = 1,
                 gehan = n1[i] + n2[i],
                 `tarone-ware` = sqrt(n1[i] + n2[i]))
    
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

# Example usage:
# h1 <- -log(0.305) / 12
# h2 <- -log(0.435) / 12
# lakatosSampleSize(24, 24, alloc = 1, h1, h2, alpha = 0.05, power = 0.8, method = "logrank", side = "two.sided")

# UI
ui <- fluidPage(
  titlePanel("Lakatos Sample Size Calculator"),
  sidebarLayout(
    sidebarPanel(
      numericInput("syear", "Survival Time (years):", value = 12),
      numericInput("yrsurv1", "Survival Probability (Standard Group):", value = 0.305),
      numericInput("yrsurv2", "Survival Probability (Test Group):", value = 0.435),
      numericInput("alloc", "Allocation Ratio (Test / Standard):", value = 1),
      numericInput("accrual", "Accrual Time (months):", value = 24),
      numericInput("follow", "Follow-up Time (months):", value = 24),
      numericInput("alpha", "Significance Level (alpha):", value = 0.05),
      numericInput("power", "Power (1 - Beta):", value = 0.8),
      selectInput("method", "Test Method:", choices = c("logrank", "gehan", "tarone-ware")),
      selectInput("side", "Test Type:", choices = c("two.sided", "one.sided")),
      actionButton("calc", "Calculate")
    ),
    mainPanel(
      verbatimTextOutput("result")
    )
  )
)

# Server
server <- function(input, output) {
  
  observeEvent(input$calc, {
    h1 <- -log(input$yrsurv1) / input$syear
    h2 <- -log(input$yrsurv2) / input$syear
    
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
    
    output$result <- renderPrint({ res })
  })
}

shinyApp(ui = ui, server = server)
