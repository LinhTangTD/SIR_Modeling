library(shiny)
library(dplyr)
library(ggplot2)
library(shinyjs)

ui <- fluidPage(
    
    # Application title
    titlePanel("Complex Model"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(id = "tPanel",style = "overflow-y:scroll; position:relative; max-height:600px",
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
            sliderInput("PS",
                        "Percentage of individuals being tested daily",
                        min = 0.1,
                        max = 1,
                        step = 0.1, 
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
                        value = 15),
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
            )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
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
    
    observeEvent(input$PTB,
                 {updateSliderInput(session, "PTA", value = 100 - input$PTB)
                 })
    
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
                        Diseased = d,
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
                        Cured = c,
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
    
    model.plot = function(model){
        cost = sum(model$Cost)
        sick = sum(model$Diseased)
        peak = max(model$Diseased)
        ggplot(data  = model, aes(x = Day)) +
            geom_line(aes(y = Healthy, color = "Healthy")) +
            geom_line(aes(y = Diseased, color = "Diseased")) +
            geom_line(aes(y = Cured, color = "Cured")) +
            scale_color_discrete(name = "Compartment") +
            labs(y = "Population", title = paste("Complex SIR Model: infected = " , sick, ", peak = ", peak, ", cost = $", cost, sep = "")) +
            theme_minimal()
    }
    
    output$distPlot <- renderPlot({
        if(input$Watch == "end"){
            day = -1
        } else {
            day = input$Days
        }
        df = complex.model(N = input$Pop, r = input$InfectionRate, g = input$SusceptibleRate, day = day,
                           d1 = input$InitialInfected, c1 = input$InitialCured, u1 = input$InitialKnown,
                           PS = input$PS, tp = input$tp, tn = input$tn, Test.cost = input$TestCost,
                           PTA = input$PTA, A.eff = input$Aeff, A.cost = input$Acost,
                           PTB = input$PTB, B.eff = input$Beff, B.cost = input$Bcost)
        model.plot(df)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)