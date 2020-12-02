library(shiny)
library(dplyr)
library(ggplot2)
library(shinyjs)

basic.model = function(n, beta, gamma, S.rate, I, day, var){
  # Calculate day 0 statistics
  curr.day = 0
  healthy = n - I
  susceptible = round(S.rate * healthy)
  infected = I
  infected.next = round(I * beta * susceptible)
  treat = round(gamma * I)
  recovered = 0
  
  # Table storing S, I, R compartments over time
  df = data_frame(Day = curr.day, Healthy = healthy, Susceptible = susceptible, Infected = infected, Recovered = recovered)
  
  # Predict S, I, R compartments by day
  while(infected*healthy > 0){
    # if number of days is specified, calculate up until that day
    if(day != -1 & curr.day >= day)
      break()
    curr.day = curr.day + 1
    healthy = healthy - infected.next
    susceptible = round(S.rate * healthy)
    infected = infected + infected.next - treat
    if(var){
      infected.next = rbinom(1, infected*susceptible, beta)
      treat = rbinom(1, infected, gamma)
    } else {
      infected.next = round(infected * beta * susceptible)
      treat = round(gamma * infected)
    }
    recovered = recovered + treat
    data = c(curr.day, healthy, susceptible, infected, recovered)
    df = rbind(df, data)
  }
  return(df)
}

model.plot = function(model){
  ggplot(data  = model, aes(x = Day)) +
    geom_line(aes(y = Healthy, color = "Healthy")) +
    geom_line(aes(y = Infected, color = "Infected")) +
    geom_line(aes(y = Recovered, color = "Recovered")) +
    scale_color_discrete(name = "Compartment") +
    labs(y = "Population", title = "Basic SIR Model") + 
    theme_minimal()
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("SIR Modeling - Britney, Linh, and Bowen"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("Pop",
                  "Population of our town",
                  min = 50,
                  max = 2000,
                  value = 500),
      checkboxInput("WatchTilEnd", "Watch Disease till end", TRUE),
      sliderInput("Days",
                   "Number of Days you want to see watch the disease",
                   min = 1,
                   max = 200,
                   step = 1, 
                   value = 10),
      sliderInput("Initialinf",
                  "Initial People infected",
                  min = 1,
                  max = 100,
                  step = 1, 
                  value = 5),
      sliderInput("Spread_Rate",
                  "Spread Rate", 
                  min = 0.01, 
                  max = 1.00, 
                  step = 0.01,
                  value = 0.01),
      sliderInput("Treat_Rate",
                  "Recovery Rate", 
                  min = 0.01, 
                  max = 1.00, 
                  step = 0.01,
                  value = 0.7),
      sliderInput("Susceptible_Rate",
                  "Susceptible Rate",
                  min = 0.01,
                  max = 1,
                  step = 0.01,
                  value = 0.2),
      checkboxInput("Var", "Add randomness into the model", FALSE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", 
                 plotOutput("Plot")),

        tabPanel("Data",
                 dataTableOutput("df")))
    )
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #Dynamic Model Input 
  observe({
    #If regression model is not selected
    if(input$WatchTilEnd){
      enable("Days")
    } else {
      disable("Days")
    }
  })

  output$Plot <- renderPlot({
    validate(
      need(input$Initialinf < input$Pop, "Initial infected number cannot be more than the population")
    )
    if(input$WatchTilEnd){
      day = -1
    } else {
      day = input$Days
    }
    df = basic.model(n = input$Pop, 
                     beta = input$Spread_Rate, 
                     gamma = input$Treat_Rate, 
                     S.rate = input$Susceptible_Rate, 
                     I = input$Initialinf, 
                     day = day,
                     var = input$Var)
    plot = model.plot(df)
    output$df = renderTable(df)
    return(plot)
 })
}

# Run the application 
shinyApp(ui = ui, server = server)

