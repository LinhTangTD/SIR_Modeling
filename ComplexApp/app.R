library(shiny)
library(dplyr)
library(ggplot2)
library(shinyjs)
library(plotly) 

basic.model = function(n, beta, gamma, S.rate, I, day, var){
    # Calculate day 0 statistics
    curr.day = 1
    healthy = n - I
    susceptible = round(S.rate * healthy)
    infected = I
    infected.next = round(I * beta * susceptible)
    treat = round(gamma * I)
    recovered = 0
    
    # Table storing S, I, R compartments over time
    df = data_frame(Day = curr.day, 
                    Healthy = healthy, 
                    Susceptible = susceptible, 
                    Infected = infected, 
                    Infected.Next = infected.next,
                    Recovering = treat,
                    Recovered = recovered)
    
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
        # data = c(curr.day, healthy, susceptible, infected, recovered)
        data = c(curr.day, healthy, susceptible, infected, infected.next, treat, recovered)
        df = rbind(df, data)
    }
    return(df)
}

complex.model = function(N, r, g, day, d1, c1, u1, PS, tp, tn, Test.cost, PTA, A.eff, A.cost, PTB, B.eff, B.cost){
    
    # Initial conditions
    curr.day = 1
    d = d1
    c = c1
    u = u1
    # Calculate initial statistics
    h = N - d - c
    hs = round(g * h)
    s = hs + d
    a = s - u
    ns = round(PS * a)
    ad = round(PS * (d - u))
    as = ns - ad
    tps = rbinom(1, as, 1-tn)
    tpd = rbinom(1, ad, tp)
    tpt = tps + tpd
    tpt1 = tpt
    NDAT = tpd + u
    NDTA = 0
    NDTB = 0
    DCA = 0
    DCB = 0
    NAT = 0
    NTA = 0
    NTB = 0
    SDC = 0
    New = min(N - c - d, round(r * hs * d))
    
    # Table storing S, I, R compartments over time
    df = data_frame(Day = curr.day, 
                    Healthy = h,
                    Symptopms = hs,
                    Infected = d,
                    Tot.Symp = s,
                    Avail.screen = a,
                    Untreated = u,
                    Screen.Tot = ns,
                    Screen.Sym = as,
                    Screen.Dis = ad,
                    Test.Pos.Sym = tps,
                    Test.Pos.Dis = tpd,
                    Test.Pos.Total = tpt,
                    Avail.TRT.Dis = NDAT,
                    TRT.A.Dis = NDTA, 
                    TRT.B.Dis = NDTB, 
                    TRT.A.Cure = DCA, 
                    TRT.B.Cure = DCB,
                    Avail.TRT.Tot = NAT,
                    TRT.A.Tot = NTA,
                    TRT.B.Tot = NTB,
                    Sum.Dis.Cured = SDC,
                    Recovered = c,
                    Catch.Dis = New)
    
    # Predict S, I, R compartments by day
    while(d*h > 0 | (day != -1 & curr.day < day)){
        # if number of days is specified, calculate up until that day
        # else calculate until there is no more diseased individual
        
        # Update current day and c, u, d
        curr.day = curr.day + 1
        c = c + SDC
        u = d - SDC
        d = d + New - SDC
        
        # Calculate other variables
        h = N - d - c
        hs = round(g * h)
        s = hs + d
        a = s - u
        ns = round(PS * a)
        ad = round(PS * (d - u))
        as = ns - ad
        tps = rbinom(1, as, 1-tn)
        tpd = rbinom(1, ad, tp)
        tpt = tps + tpd
        NDAT = tpd + u
        NDTA = round(PTA * NDAT)
        NDTB = min(NDAT - NDTA, round(PTB * NDAT))
        DCA = min(NDTA, round(NDTA * A.eff))
        DCB = min(NDTB, round(NDTB * B.eff))
        if (curr.day == 2){
            NAT = tpt + tpt1
        }
        else{
            NAT = tpt + u
        }
        NTA = round(PTA * NAT)
        NTB = NAT - NTA
        SDC = DCA + DCB
        New = min(N - c - d, round(r * hs * d))
        
        data = c(curr.day, h, hs, d, s, a, u,
                 ns, as, ad, tps, tpd, tpt,
                 NDAT, NDTA, NDTA, DCA, DCB, NAT, NTA, NTB,
                 SDC, c, New)
        df = rbind(df, data)
    }
    df$Cost = df$TRT.A.Tot * A.cost + df$TRT.B.Tot * B.cost + df$Screen.Tot * Test.cost
    return(df)
}

model.plotly = function(model){
  accumulate_by <- function(dat, var) {
    var <- lazyeval::f_eval(var, dat)
    lvls <- plotly:::getLevels(var)
    dats <- lapply(seq_along(lvls), function(x) {
      cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
    })
    dplyr::bind_rows(dats)
  }
  visual_data = select(model, c(Day, Healthy, Infected, Recovered)) %>%
    tidyr::gather(key = Compartment, value = Population, Healthy, Infected, Recovered) %>% 
    accumulate_by(~Day)
  fig = plot_ly(data = visual_data, type = "scatter", mode = "lines",
                x = ~Day, y = ~Population, color = ~Compartment, frame = ~frame) %>% 
          layout(hovermode = 'compare') %>%
          animation_opts(frame = 100, transition = 0, redraw = FALSE) %>%
          animation_slider(hide = T) %>%
          animation_button(x = 1, xanchor = "right", y = -0.1, yanchor = "bottom", label = "View Plot")
  return(fig)
}

model.plot = function(model){
  ggplot(data  = model, aes(x = Day)) +
    geom_line(aes(y = Healthy, color = "Healthy")) +
    geom_line(aes(y = Infected, color = "Infected")) +
    geom_line(aes(y = Recovered, color = "Recovered")) +
    scale_color_discrete(name = "Compartment") +
    labs(y = "Population", title = "SIR Model") + 
    theme_minimal()
}

model.rep.plot = function(reps){
  line.plot = function(model){
    list(geom_point(data = model, aes(x = Day, y = Healthy, color = "Healthy"), size = 0.8),
         geom_point(data = model, aes(x = Day, y = Diseased, color = "Diseased"), size = 0.8),
         geom_point(data = model, aes(x = Day, y = Cured, color = "Cured"), size = 0.8),
         geom_line(data = model, aes(x = Day, y = Healthy, color = "Healthy"), size = 0.1),
         geom_line(data = model, aes(x = Day, y = Diseased, color = "Diseased"), size = 0.1),
         geom_line(data = model, aes(x = Day, y = Cured, color = "Cured"), size = 0.1))
  }
  plot = ggplot() + labs(y = "Population", title = paste("Complex SIR Model with", reps, "reps")) + theme_minimal()
  for(i in 1:reps){
      df = complex.model.rand()
      plot = plot + line.plot(df)
  }
  plot = plot + scale_color_discrete(name = "Compartment")
  plot
  return(plot)
}

ui <- fluidPage(
    
    # Application title
    titlePanel("SIR Modeling"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            id = "tPanel", 
            style = "overflow-y:scroll; position:relative; max-height:900px",
            radioButtons("Model", "Type of Model",
                         c("Basic SIR" = "basic", "Complex SIR with treatment & testing" = "complex")),
                         # c("Basic SIR" = "basic")),
            radioButtons("Watch", "Length of Observation:",
                         c("Watch disease till end" = "end", "Watch for specific number of days" = "day")),
            conditionalPanel(
                condition = "input.Watch == 'day'",
                sliderInput("Days",
                            NULL,
                            min = 1,
                            max = 200,
                            step = 1,
                            value = 38)
            ),
            sliderInput("Pop",
                        "Population",
                        min = 50,
                        max = 1000,
                        value = 300),
            sliderInput("InitialInfected",
                        "Initial Infected Individuals",
                        min = 1,
                        max = 100,
                        step = 1, 
                        value = 5),
            sliderInput("InfectionRate",
                        "Tranmission/Infection Rate", 
                        min = 0.01, 
                        max = 0.4, 
                        step = 0.01,
                        value = 0.02),
            sliderInput("SusceptibleRate",
                        "Percentage of healthy individuals with symptoms",
                        min = 0.01,
                        max = 1,
                        step = 0.01,
                        value = 0.2),
            conditionalPanel(
                condition = "input.Model == 'basic'",
                sliderInput("RecoveryRate",
                            "Recovery Rate",
                            min = 0.01,
                            max = 1.00,
                            step = 0.01,
                            value = 0.7),
                checkboxInput("Var", "Add randomness into the model", FALSE)
            ),
            conditionalPanel(
                condition = "input.Model == 'complex'",
                sliderInput("InitialCured",
                            "Initial Cured Individuals",
                            min = 0,
                            max = 100,
                            step = 1,
                            value = 0),
                sliderInput("InitialKnown",
                            "Initial Known Cases",
                            min = 0,
                            max = 100,
                            step = 1,
                            value = 0),
                sliderInput("PS",
                            "Percentage of individuals being tested daily",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 1),
                sliderInput("tp",
                            "Test sensitivity (%)",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 0.99),
                sliderInput("tn",
                            "Test specificity (%)",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 0.95),
                numericInput("TestCost",
                             "Testing Cost ($)",
                             min = 1,
                             max = 1000,
                             step = 1,
                             value = 1),
                sliderInput("PTA",
                            "Percentage of individuals receiving Treatment A",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 1),
                sliderInput("PTB",
                            "Percentage of individuals receiving Treatment B",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 0),
                sliderInput("Aeff",
                            "Treatment A effectiveness (%)",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 0.75),
                sliderInput("Beff",
                            "Treatment B effectiveness (%)",
                            min = 0.1,
                            max = 1,
                            step = 0.01,
                            value = 0.25),
                numericInput("Acost",
                             "Treatment A cost ($)",
                             min = 1,
                             max = 100,
                             step = 1,
                             value = 20),
                numericInput("Bcost",
                             "Treatment B cost ($)",
                             min = 1,
                             max = 100,
                             step = 1,
                             value = 15)
            )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          id = "tPanel", 
          style = "overflow-y:scroll; position:relative; max-height:900px",
            tabsetPanel(
                tabPanel("Plot", 
                         plotOutput("Plot"),
                         # plotlyOutput("Plot1", width = "100%", height= "100%"),
                         htmlOutput("model_info")),
                tabPanel("Data",
                         checkboxInput("FullTable", "See Full Table", TRUE),
                         dataTableOutput("df")), 
                tabPanel("Basic Model Info",
                         includeHTML("basic.model.html")),
                tabPanel("Complex Model Info",
                         includeHTML("complex.model.html")))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #Dynamic Model Input 
    observeEvent(input$PTA, 
                 {updateSliderInput(session, "PTB", value = 1 - input$PTA)
                 })
    
    observeEvent(input$PTB,
                 {updateSliderInput(session, "PTA", value = 1 - input$PTB)
                 })
    
    observeEvent(input$Pop,
                 {updateSliderInput(session, "InitialInfected", max = input$Pop)
                     observeEvent(input$InitialInfected,
                                  {updateSliderInput(session, "InitialCured", max = input$InitialInfected)
                                      updateSliderInput(session, "InitialKnown", max = input$InitialInfected)
                                  }
                     )
                 })
    
    output$Plot <- renderPlot({
        if(input$Watch == "end"){
            day = -1
        } else {
            day = input$Days
        }
        
        if(input$Model == "basic"){
            df = basic.model(n = input$Pop, 
                             beta = input$InfectionRate, 
                             gamma = input$RecoveryRate, 
                             S.rate = input$SusceptibleRate, 
                             I = input$InitialInfected, 
                             day = day,
                             var = input$Var)
        } else {
            df = complex.model(N = input$Pop, r = input$InfectionRate, g = input$SusceptibleRate, day = day,
                               d1 = input$InitialInfected, c1 = input$InitialCured, u1 = input$InitialKnown,
                               PS = input$PS, tp = input$tp, tn = input$tn, Test.cost = input$TestCost,
                               PTA = input$PTA, A.eff = input$Aeff, A.cost = input$Acost,
                               PTB = input$PTB, B.eff = input$Beff, B.cost = input$Bcost)
        }
        
        observeEvent(input$FullTable,
                     {if(input$FullTable){
                         output$df = renderDataTable({df}, options = list(scrollX = TRUE))
                     } else {
                         small_df = select(df, c(Day, Healthy, Infected, Recovered))
                         output$df = renderDataTable({small_df}, options = list(scrollX = TRUE))
                     }}
                     )
        
        output$model_info = renderUI({
            peak_infected = max(df$Infected)
            peak_day = df[df$Infected == peak_infected,]$Day[1]
            total_infected = sum(df$Infected)
            line1 = paste("The peak number of infected individuals is ", peak_infected, " at day ", peak_day, ".", sep = "")
            line2 = paste("The total number of infected individuals is ", total_infected, ".", sep = "")
            content = paste(line1, line2, sep = "<br/>")
            if(input$Model == "complex"){
                total_cost = sum(df$Cost)
                line3 = paste("The total cost of treatment and testing is $", total_cost, sep = "")
                content = paste(content, line3, sep = "<br/>")
            }
            HTML(content)
        })
        
        plot = model.plot(df)
        # output$Plot1 = renderPlotly({model.plotly(df)})
        return(plot)
    })
    
    # output$markdown <- renderUI({
    #     HTML(markdown::markdownToHTML(knitr::knit('../basic.model.rmd', quiet = TRUE)))
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)