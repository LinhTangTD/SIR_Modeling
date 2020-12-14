library(shiny)
library(dplyr)
library(ggplot2)
library(shinyjs)
library(plotly) 

SIR.model = function(n, beta, gamma, S.rate, I, day, var){
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

SIR.plot = function(model){
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

#--------NEW MODEL------------#
SIHRD.data = function(N, beta, I.0, p, epsilon, theta, gamma, TimeH, TimeNH, watch.day){
  I.rate = beta/N
  psi = 1/TimeH
  lambda = 1/TimeNH
  curr.day = 0
  S = p * N - I.0
  NS.0 = N - S - I.0
  NS = NS.0
  NewInfected = I.0
  Total.Hospital = 0
  Total.NonHospital = 0
  Total.Infected = NewInfected + Total.Hospital + Total.NonHospital
  Total.Infective = NewInfected + Total.NonHospital
  Hospital.Recovered = 0
  NonHospital.Recovered = 0
  Hospital.Dead = 0
  NonHospital.Dead = 0
  Total.Recovered = Hospital.Recovered + NonHospital.Recovered
  Total.Dead = Hospital.Dead + NonHospital.Dead
  Pop = S + NS + Total.Infected + Total.Dead
  
  # If you want a full table like the excel, uncomment all lines in the df below
  df = data_frame(Day = curr.day,
                  Susceptible = S,
                  Non.Susceptible = NS,
                  Daily.New.Case = NewInfected,
                  Total.In.Hospital = Total.Hospital,
                  # Total.NonHospital = Total.NonHospital,
                  Infected = Total.Infected,
                  # Total.Infective = Total.Infective,
                  # Recovered.In.Hospital = Hospital.Recovered,
                  # Dead.In.Hospital = Hospital.Dead,
                  Total.Recovered = Total.Recovered,
                  Total.Dead = Total.Dead)
  
  # Watch the disease until it ends (daily new case = 0) or at max 200 days
  while(NewInfected * S > 0 & curr.day < watch.day){
    curr.day = curr.day + 1
    
    # Number of individuals recovered from the infected in hospital (previous day)
    Hospital.Recovered = Hospital.Recovered + round(theta * round(Total.Hospital * psi))
    Hospital.Dead = Hospital.Dead + round(Total.Hospital * psi) - round(theta * round(Total.Hospital * psi))
    
    # Number of individuals dead from the infected in hospital (previous day)
    NonHospital.Recovered = NonHospital.Recovered + round(gamma * round(Total.NonHospital * lambda))
    NonHospital.Dead = NonHospital.Dead + round(Total.NonHospital * lambda) - round(gamma * round(Total.NonHospital * lambda))
    
    # Total number of recovered or dead so far (cumulative)
    Total.Recovered = Hospital.Recovered + NonHospital.Recovered
    Total.Dead = Hospital.Dead + NonHospital.Dead
    
    # Number of individuals hospitalized from the infected in the previous day
    Total.Hospital = round(NewInfected * epsilon) + Total.Hospital - round(Total.Hospital * psi)
    Total.NonHospital = NewInfected - round(NewInfected * epsilon) + Total.NonHospital - round(Total.NonHospital * lambda)
    
    # Newly infected case today and total actively infected so far
    NewInfected = round(I.rate * S * Total.Infective)
    Total.Infected = NewInfected + Total.Hospital + Total.NonHospital
    
    # Susceptible population left
    S = S - NewInfected
    
    # Current non-susceptible group
    NS = NS.0 + Total.Recovered
    
    # Infective group that can spread the disease the next day
    Total.Infective = NewInfected + Total.NonHospital
    
    # Confirm population
    Pop = S + NS + Total.Infected + Total.Dead
    
    # If you want a full table like the excel, uncomment this
    # data = c(curr.day, S, NS, NewInfected, Total.Hospital, Total.NonHospital, Total.Infected, Total.Infective, Hospital.Recovered, Hospital.Dead, Total.Recovered, Total.Dead, Pop)
    data = c(curr.day, S, NS, NewInfected, Total.Hospital, Total.Infected, Total.Recovered, Total.Dead)
    df = rbind(df, data)
  }
  return(df)
}

SIHRD.plot = function(data){
  p = ggplot(data  = data, aes(x = Day)) +
    geom_line(aes(y = Susceptible, color = "Susceptible")) +
    geom_line(aes(y = Daily.New.Case, color = "Daily New Case")) +
    geom_line(aes(y = Infected, color = "Infected")) +
    geom_line(aes(y = Total.In.Hospital, color = "Currently In Hospital")) +
    geom_line(aes(y = Total.Recovered, color = "Total Recoverd")) +
    geom_line(aes(y = Total.Dead, color = "Total Dead")) +
    scale_color_discrete(name = "Compartment") +
    labs(y = "Population", title = "SIHRD Model") + 
    theme_minimal()
  return(p)
}

ui <- fluidPage(
    
    # Application title
    titlePanel("Infectious Disease Modeling"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            id = "tPanel", 
            style = "overflow-y:scroll; position:relative; max-height:900px",
            radioButtons("Model", 
                         "Type of Model",
                         c("SIHRD" = "SIHRD", 
                           # "SIHRD (+TT)" = "TT",
                           "SIR" = "SIR")),
            # radioButtons("Watch", 
            #              "Length of Observation:",
            #              c("Watch disease till end" = "end", "Watch for specific number of days" = "day")),
            # conditionalPanel(
            #     condition = "input.Watch == 'day'",
            #     sliderInput("Days",
            #                 NULL,
            #                 min = 1,
            #                 max = 200,
            #                 step = 1,
            #                 value = 38)
            # ),
            sliderInput("Days",
                        "Watch the disease for number of days",
                        NULL,
                        min = 1,
                        max = 200,
                        step = 1,
                        value = 100),
            sliderInput("Pop",
                        "Population",
                        min = 0,
                        max = 10000,
                        step = 1000,
                        value = 3000),
            sliderInput("InitialInfected",
                        "Initial Infected Individuals",
                        min = 1,
                        max = 100,
                        step = 1, 
                        value = 5),
            sliderInput("SusceptibleRate",
                        "Percentage of population susceptible to getting disease",
                        min = 0,
                        max = 1,
                        step = 0.05,
                        value = 0.5),
            conditionalPanel(
                condition = "input.Model == 'SIR'",
                sliderInput("RecoveryRate",
                            "Recovery Rate",
                            min = 0,
                            max = 1.00,
                            step = 0.01,
                            value = 0.7),
                sliderInput("InfectionRate",
                            "Infection Rate", 
                            min = 0.01, 
                            max = 0.3, 
                            step = 0.001,
                            value = 0.02),
                checkboxInput("Var", "Add randomness into the model", FALSE)
            ),
            conditionalPanel(
              condition = "input.Model == 'SIHRD'",
              # beta, I.0, p, epsilon, theta, gamma, TimeH, TimeNH
              sliderInput("beta",
                          "Number of people an infected person can transmit the disease to", 
                          min = 1, 
                          max = 10, 
                          step = 1,
                          value = 3),
              sliderInput("epsilon",
                          "Hospitalization Rate",
                          min = 0.1,
                          max = 1,
                          step = 0.05,
                          value = 0.6),
              sliderInput("theta",
                          "Recovery Rate in Hospital",
                          min = 0.1,
                          max = 1,
                          step = 0.05,
                          value = 0.9),
              sliderInput("gamma",
                          "Recovery Rate without Hospitalization",
                          min = 0.1,
                          max = 1,
                          step = 0.05,
                          value = 0.8),
              sliderInput("TimeH",
                          "Infection Time in Hospitalization",
                          min = 1,
                          max = 20,
                          step = 1,
                          value = 10),
              sliderInput("TimeNH",
                          "Infection Time without Hospitalization",
                          min = 0.1,
                          max = 20,
                          step = 1,
                          value = 5)
            )
            # conditionalPanel(
            #     condition = "input.Model == 'TT'",
            #     sliderInput("InitialCured",
            #                 "Initial Cured Individuals",
            #                 min = 0,
            #                 max = 100,
            #                 step = 1,
            #                 value = 0),
            #     sliderInput("InitialKnown",
            #                 "Initial Known Cases",
            #                 min = 0,
            #                 max = 100,
            #                 step = 1,
            #                 value = 0),
            #     sliderInput("PS",
            #                 "Percentage of individuals being tested daily",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 1),
            #     sliderInput("tp",
            #                 "Test sensitivity (%)",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 0.99),
            #     sliderInput("tn",
            #                 "Test specificity (%)",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 0.95),
            #     numericInput("TestCost",
            #                  "Testing Cost ($)",
            #                  min = 1,
            #                  max = 1000,
            #                  step = 1,
            #                  value = 1),
            #     sliderInput("PTA",
            #                 "Percentage of individuals receiving Treatment A",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 1),
            #     sliderInput("PTB",
            #                 "Percentage of individuals receiving Treatment B",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 0),
            #     sliderInput("Aeff",
            #                 "Treatment A effectiveness (%)",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 0.75),
            #     sliderInput("Beff",
            #                 "Treatment B effectiveness (%)",
            #                 min = 0.1,
            #                 max = 1,
            #                 step = 0.01,
            #                 value = 0.25),
            #     numericInput("Acost",
            #                  "Treatment A cost ($)",
            #                  min = 1,
            #                  max = 100,
            #                  step = 1,
            #                  value = 20),
            #     numericInput("Bcost",
            #                  "Treatment B cost ($)",
            #                  min = 1,
            #                  max = 100,
            #                  step = 1,
            #                  value = 15)
            # )
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
                         conditionalPanel(
                           condition = "input.Model == 'SIR'",
                           checkboxInput("FullTable", "See Full Table", TRUE)),
                         dataTableOutput("df")),
                tabPanel("SIR Model Info",
                         includeHTML("basic.model.html")))
                # tabPanel("SIHRD Model Info",
                #          includeHTML("SIRHD.model.html")),
                # tabPanel("SIHRD (+TT) Model Info",
                #          includeHTML("complex.model.html")))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #Dynamic Complex Model Input 
    # observeEvent(input$PTA, 
    #              {updateSliderInput(session, "PTB", value = 1 - input$PTA)
    #              })
    # 
    # observeEvent(input$PTB,
    #              {updateSliderInput(session, "PTA", value = 1 - input$PTB)
    #              })
    
    observeEvent(input$Pop,
                 {updateSliderInput(session, "InitialInfected", max = input$Pop)
                     observeEvent(input$InitialInfected,
                                  {updateSliderInput(session, "InitialCured", max = input$InitialInfected)
                                      updateSliderInput(session, "InitialKnown", max = input$InitialInfected)
                                  }
                     )
                     updateSliderInput(session, "beta", max = input$Pop)
                 })
    
    output$Plot <- renderPlot({
        # if(input$Watch == "end"){
        #     day = -1
        # } else {
        #     day = input$Days
        # }
        
        if(input$Model == "SIR"){
            df = SIR.model(n = input$Pop, 
                           beta = input$InfectionRate, 
                           gamma = input$RecoveryRate, 
                           S.rate = input$SusceptibleRate, 
                           I = input$InitialInfected, 
                           day = input$Days,
                           var = input$Var)
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
            plot = SIR.plot(df)
        } else if (input$Model == "SIHRD"){
          df = SIHRD.data(N = input$Pop, 
                          beta = input$beta, 
                          I.0 = input$InitialInfected, 
                          p = input$SusceptibleRate, 
                          epsilon = input$epsilon, 
                          theta = input$theta, 
                          gamma = input$gamma, 
                          TimeH = input$TimeH, 
                          TimeNH = input$TimeNH,
                          watch.day = input$Days)
          output$df = renderDataTable({df}, options = list(scrollX = TRUE))
          plot = SIHRD.plot(df)
          output$model_info = renderUI({
            Last.Day = nrow(df)
            Status = ifelse(df[Last.Day,]$Daily.New.Case == 0, "Ended", "Going On")
            Total.Death = df[Last.Day,]$Total.Dead
            Total.Infected = sum(df$Infected)
            Max.New.Case = max(df$Daily.New.Case)
            Peak.Day.New.Case = df[df$Daily.New.Case == Max.New.Case,]$Day
            Max.In.Hospital = max(df$Total.In.Hospital)
            Peak.Day.Hospital= df[df$Total.In.Hospital == Max.In.Hospital,]$Day
            
            peak_infected = max(df$Infected)
            peak_day = df[df$Infected == peak_infected,]$Day[1]
            total_infected = sum(df$Infected)
            line1 = paste("The peak number of infected individuals is ", Max.New.Case, " at day ", Peak.Day.New.Case, ".", sep = "")
            line2 = paste("The peak number of individuals in the hospital is ", Max.In.Hospital, " at day ", Peak.Day.Hospital, ".", sep = "")
            line3 = paste("The total number of infected individuals is ", Total.Infected, ".", sep = "")
            line4 = paste("The total number of death is ", Total.Death, ".", sep = "")
            line5 = ifelse(df[Last.Day,]$Daily.New.Case == 0, 
                           paste("The disease ended on day ", Last.Day, ".", sep = ""),
                           paste("On the last day being observed, day ", Last.Day, ", the disease is still going on.", sep = ""))
            content = paste(line1, line2, line3, line4, line5, sep = "<br/>")
            HTML(content)
          })
        } 
      # else {
      #       df = complex.model(N = input$Pop, r = input$InfectionRate, g = input$SusceptibleRate, day = day,
      #                          d1 = input$InitialInfected, c1 = input$InitialCured, u1 = input$InitialKnown,
      #                          PS = input$PS, tp = input$tp, tn = input$tn, Test.cost = input$TestCost,
      #                          PTA = input$PTA, A.eff = input$Aeff, A.cost = input$Acost,
      #                          PTB = input$PTB, B.eff = input$Beff, B.cost = input$Bcost)
      #   }
        return(plot)
    })
    
    # output$markdown <- renderUI({
    #     HTML(markdown::markdownToHTML(knitr::knit('../basic.model.rmd', quiet = TRUE)))
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)