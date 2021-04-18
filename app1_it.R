###############
# Covid-19 IT #
###############

library(curl)
library(shiny)
library(dplyr)
library(dygraphs)
library(KFAS)
library(shinythemes)

dt <- jsonlite::fromJSON(txt = "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-json/dpc-covid19-ita-regioni.json")
regioni <- unique(dt$denominazione_regione)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),
   
   # Application title
   titlePanel("Report sui dati Covid-19 in Italia"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         selectInput("regione",
                     label    = "Seleziona regione",
                     choices  = c("Italia", regioni),
                     selected = "Italia"),
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
                  l'ultima stima del tasso di crescita. Esistono,",
                  a("modelli epidemici",
                    href = "https://it.wikipedia.org/wiki/Modelli_matematici_in_epidemiologia"),
                  "per prevedere l'andamento delle
                  epidemie. Tuttavia, dato che i dati sono parziali
                  e incerti, e dato che i comportamenti degli Italiani e
                  delle istituzioni sanitarie sono cambiati diverse
                  volte per via di innovazioni nei decreti, regolamenti e
                  politiche sanitarie, qui ci limitiamo ad estrarre regolaritÃ 
                  dai dati senza un vero e proprio modello di sviluppo
                  dell'epidemia. Queste proiezioni sono corredate di
                  intervalli di confidenza, che evidenziano
                  l'enorme incertezza delle previsioni a lungo termine."),
         sliderInput("orizzonte",
                     label = "Orizzonte previsivo in giorni",
                     min = 0,
                     max = 30,
                     value = 0,
                     round = TRUE),
         h4("Data di azzeramento"),
         textOutput("data_zero"),
         helpText("Elaborazioni di Matteo Pelagatti sui dati della
                   da Protezione Civile pubblicati su",
                   a("https://github.com/pcm-dpc/COVID-19",
                     href = "https://github.com/pcm-dpc/COVID-19"),
                  ". Si consulti",
                  a("https://matteopelagatti.github.io/Covid-19/rapporto_italia_short",
                    href = "https://matteopelagatti.github.io/Covid-19/rapporto_italia_short"),
                  "per maggiori dettagli sul modello usato.")
      ),
      
      # Show the plots of series, level and rate
      mainPanel(
         dygraphOutput("line_plot"),
         dygraphOutput("slope_plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  dataInput <- reactive({
    dt <- jsonlite::fromJSON(txt = "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-json/dpc-covid19-ita-regioni.json")
    dt$data <- as.Date(strtrim(dt$data, 10))
    dt
  })
  
  dataSelect <- reactive({
    if (input$regione == "Italia") {
      dt <- dataInput() %>%
        group_by(data) %>%
        select_if(is.numeric) %>%
        summarise_all(sum) %>% as.data.frame()
    } else {
      dt <- dataInput() %>%
        dplyr::filter(denominazione_regione == input$regione)
    }
    if (input$serie == "Nuovi positivi") {
      dt <- dt[, c("data", "nuovi_positivi")]
    } else if (input$serie == "Decessi") {
      dt <- data.frame(data = dt$data, decessi = c(0, diff(dt$deceduti)))
    } else if (input$serie == "Positivi in terapia intensiva") {
      dt <- dt[, c("data", "terapia_intensiva")]
    } else if (input$serie == "Positivi ospedalizzati") {
      dt <- dt[, c("data", "totale_ospedalizzati")]
    }
    dt
  })
  
  ssModel <- reactive({
    dt <- dataSelect()
    ly <- log(dt[, 2] + 0.1)
    if (input$orizzonte) ly <- c(ly, rep(NA, input$orizzonte))
    mod <- SSModel(ly ~ SSMtrend(2, Q = list(0, NA)) +
                     SSMseasonal(7, NA),
                   H = NA)
    vly <- var(diff(ly, 7), na.rm = TRUE)
    init <- log(c(vly / 10, vly / 10, vly / 2))
    fit <- fitSSM(mod, init)
    smo <- KFS(fit$model, smoothing = c("state", "signal"))
    smo$date <- dt$data
    if (input$orizzonte) {
      smo$date <- c(smo$date, dt$data[nrow(dt)] + 1:input$orizzonte)
    }
    smo
  })
  
  observeEvent(input$carica_dati, {
    dataInput()
    output$ora <- renderText(as.character(Sys.time()))
  })
  
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
                   ") la data di azzeramento Ã¨ indeterminata")
    } else if (ultimo_livel < 1) {
      out <- paste("L'azzeramento Ã¨ avvenuto in data",
                   smo$date[-(1:35)][
                     which(smo$alphahat[-(1:35), "level"] < 0)[1]
                     ])
    } else {
      h <- round(-log(ultimo_livel) / log(ultimo_tasso), 0)
      out <- paste0("Se il tasso rimane l'ultimo stimato (", round(ultimo_tasso, 2),
                   ") la data di azzeramento Ã¨ ", as.character(ultima_data + h))
    }
    out
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
