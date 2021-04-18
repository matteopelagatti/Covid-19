###############
# Covid-19 CH #
###############


library(curl)
library(shiny)
library(dplyr)
library(dygraphs)
library(KFAS)
library(shinythemes)


dt <- read.csv("https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv",
               stringsAsFactors = FALSE)
regioni <- names(dt)[2:28]


# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),
                
                # Application title
                titlePanel("Report sui dati Covid-19 in Svizzera"),
                
                # Sidebar with a slider input for number of bins 
                sidebarLayout(
                    sidebarPanel(
                        #actionButton("carica_dati", label = "Ricarica i dati"),
                        textOutput("ora"),
                        selectInput("regione",
                                    label    = "Seleziona il cantone",
                                    choices  = regioni,
                                    selected = "CH"),
                        selectInput("serie",
                                    label    = "Seleziona serie storica",
                                    choices  = c("Nuovi positivi",
                                                 "Decessi",
                                                 "Positivi in terapia intensiva",
                                                 "Positivi ospedalizzati"),
                                    selected = "Nuovi positivi"),
                        sliderInput("conf",
                                    label = "Livello di confidenza",
                                    min   = 0,
                                    max   = 100,
                                    value = 50,
                                    step  = 5,
                                    round = TRUE),
                        h3("Previsioni"),
                        helpText("Queste previsioni proiettano nel futuro
                                  l'ultima stima del tasso di crescita. Esistono modelli
                                  epidemici (come il SIR) per prevedere l'andamento delle
                                  epidemie. Qui ci limitiamo ad estrarre regolarità
                                  dai dati senza un vero e proprio modello di sviluppo
                                  dell'epidemia. Queste proiezioni sono corredate di
                                  intervalli di confidenza, utilissimi a evidenziare
                                  l'enorme incertezza che la previsione porta con sé."),
                        sliderInput("orizzonte",
                                    label = "Orizzonte previsivo in giorni",
                                    min = 0,
                                    max = 30,
                                    value = 0,
                                    round = TRUE),
                        h4("Data di azzeramento"),
                        textOutput("data_zero"),
                        helpText("Elaborazioni di Matteo Pelagatti sui dati collezionati da
                   Daniel Probst e pubblicati qui: https://github.com/daenuprobst/covid19-cases-switzerland/")
                ),
                    
                    # Show a plot of the generated distribution
                    mainPanel(
                        dygraphOutput("line_plot"),
                        dygraphOutput("slope_plot")
                    )
                )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    dataInput <- reactive({
        file <- switch(input$serie,
                       "Nuovi positivi" = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv",
                       "Decessi" = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_fatalities_switzerland_openzh.csv",
                       "Positivi in terapia intensiva" = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_icu_switzerland_openzh.csv",
                       "Positivi ospedalizzati" = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_hospitalized_switzerland_openzh.csv"
        )
        dt <- read.csv(file,
                       stringsAsFactors = FALSE)
        dt$Date <- as.Date(dt$Date)
        dt
    })
    
    dataSelect <- reactive({
        dt <- dataInput()[, c("Date", input$regione)]
        if (input$serie %in% c("Nuovi positivi", "Decessi")) {
            dt[, 2] <- c(NA, diff(dt[, 2]))
        }
        dt
    })
    
    ssModel <- reactive({
        dt <- dataSelect()
        ly <- log(dt[, 2] + 0.1)
        if (input$orizzonte) ly <- c(ly, rep(NA, input$orizzonte))
        if (input$regione == "CH") {
            mod <- SSModel(ly ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(7, NA),
                           H = NA)
            fit <- fitSSM(mod, c(0.1, 2*log(sd(ly, na.rm = TRUE)/2), 0))
        } else {
            mod <- SSModel(ly ~ SSMtrend(2, Q = list(0, NA)),
                           H = NA)
            fit <- fitSSM(mod, c(0.1, 2*log(sd(ly, na.rm = TRUE)/2)))
        }
        smo <- KFS(fit$model, smoothing = c("state", "signal"))
        smo$date <- dt$Date
        if (input$orizzonte) {
            smo$date <- c(smo$date, dt$Date[nrow(dt)] + 1:input$orizzonte)
        }
        smo
    })
    
    # observeEvent(input$carica_dati, {
    #     dataInput()
    #     output$ora <- renderText(as.character(Sys.time()))
    # })
    
    output$line_plot <- renderDygraph({
        z <- -qnorm((100 - input$conf)/200)
        smo <- ssModel()
        oss <- exp(as.numeric(smo$model$y)) - 0.1
        med <- exp(as.numeric(smo$alphahat[, "level"]))
        lwr <- med * exp(-z * sqrt(as.numeric(smo$V[1, 1, ])))
        upr <- med * exp( z * sqrt(as.numeric(smo$V[1, 1, ])))
        xts::xts(
            x = cbind(inf    = lwr,
                      attesi = med,
                      sup    = upr,
                      oss    = oss),
            order.by = smo$date
        ) %>%
            dygraph(main = input$serie) %>%
            dySeries(c("inf", "attesi", "sup"), label = "livello") %>%
            dySeries("oss", label = "osservati") %>%
            dyShading(from = smo$date[length(smo$date)] - input$orizzonte,
                      to   = smo$date[length(smo$date)])
    })
    
    output$slope_plot <- renderDygraph({
        z <- -qnorm((100 - input$conf)/200)
        smo <- ssModel()
        mediana <- exp(as.numeric(smo$alphahat[, "slope"]))
        lwr <- mediana * exp(-z * sqrt(as.numeric(smo$V[2, 2, ])))
        upr <- mediana * exp( z * sqrt(as.numeric(smo$V[2, 2, ])))
        xts::xts(
            x = cbind(inf = lwr, tasso = mediana, sup = upr),
            order.by = smo$date
        ) %>%
            dygraph(main = "Tasso lordo giornaliero") %>%
            dySeries(c("inf", "tasso", "sup"), label = "tasso") %>%
            dyShading(from = smo$date[length(smo$date)] - input$orizzonte,
                      to   = smo$date[length(smo$date)])
    })
    
    output$data_zero <- renderText({
        smo <- ssModel()
        ultima_data  <- smo$date[length(smo$date)]
        ultimo_livel <- exp(smo$alphahat[nrow(smo$alphahat), "level"])
        ultimo_tasso <- exp(smo$alphahat[nrow(smo$alphahat), "slope"])
        if (ultimo_tasso >= 1) {
            out <- paste0("Se il tasso rimane l'ultimo stimato (", round(ultimo_tasso, 2),
                          ") la data di azzeramento è indeterminata")
        } else {
            h <- round(-log(ultimo_livel) / log(ultimo_tasso), 0)
            out <- paste0("Se il tasso rimane l'ultimo stimato (", round(ultimo_tasso, 2),
                          ") la data di azzeramento è ", as.character(ultima_data + h))
        }
        out
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)