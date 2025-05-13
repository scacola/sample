library(shiny)
library(DT)

# --- Sample Size Calculation Function ---
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
  
  data.frame(
    Metric = c(
      "Total Sample Size",
      "Standard Group Sample Size",
      "Test Group Sample Size",
      "Expected Events (Standard Group)",
      "Expected Events (Test Group)",
      "Total Expected Events"
    ),
    Value = c(
      sample_std + sample_test,
      sample_std,
      sample_test,
      round(n * p1 * (1 + d1), 1),
      round(n * p2 * (1 + d2), 1),
      round(n * p1 * (1 + d1) + n * p2 * (1 + d2), 1)
    )
  )
}

# --- Shiny App ---
ui <- fluidPage(
  titlePanel("Sample Size for Log-rank Non-Inferiority Test"),
  sidebarLayout(
    sidebarPanel(
      numericInput("syear", "Survival Time (year):", value = 16),
      numericInput("yrsurv1", "Survival Probability (Standard Group):", value = 0.5, min = 0.01, max = 0.99),
      numericInput("yrsurv2", "Survival Probability (Test Group):", value = 0.3, min = 0.01, max = 0.99),
      numericInput("alloc", "Allocation Ratio (Test / Standard):", value = 2),
      numericInput("accrual", "Accrual Time (months):", value = 24),
      numericInput("follow", "Follow-up Time (months):", value = 24),
      numericInput("alpha", "One-sided Alpha:", value = 0.025),
      numericInput("power", "Power (1 - Beta):", value = 0.6, min = 0.01, max = 1),
      numericInput("margin", "Non-inferiority Margin (HR):", value = 1.3),
      actionButton("calc", "Calculate")
    ),
    mainPanel(
      DTOutput("result_table"),
      tags$hr(),
      tags$div("Reference: Jung SH, Chow SC. On sample size calculation for comparing survival curves under general hypothesis testing. Journal of Biopharmaceutical Statistics 2012; 22(3):485–495.")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$calc, {
    h1 <- -log(input$yrsurv1) / input$syear
    h2 <- -log(input$yrsurv2) / input$syear
    beta <- 1 - input$power
    
    res <- twoSurvSampleSizeNI(input$accrual, input$follow, input$alloc, h1, h2, input$alpha, beta, input$margin)
    
    output$result_table <- renderDT({
      datatable(res,
                colnames = NULL,
                options = list(
                  dom = 't',
                  ordering = FALSE,
                  columnDefs = list(list(className = 'dt-left', targets = '_all'))
                ),
                rownames = FALSE
      )
    })
  })
}

shinyApp(ui, server)
