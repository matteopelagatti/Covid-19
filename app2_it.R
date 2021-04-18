library(shiny)
library(shinythemes)
library(tidyverse)
library(KFAS)
library(xts)
library(dygraphs)
library(DT)

pop <- tribble(
    ~Regione, ~Popolazione,
    "Abruzzo",                1293941,
    "Basilicata",              553254,
    "Calabria",               1894110,
    "Campania",               5712143,
    "Emilia-Romagna",         4464119,
    "Friuli Venezia Giulia",  1206216,
    "Lazio",                  5755700,
    "Liguria",                1524826,
    "Lombardia",             10027602,
    "Marche",                 1512672,
    "Molise",                  300516,
    "P.A. Bolzano",            532644,
    "P.A. Trento",             545425,
    "Piemonte",               4311217,
    "Puglia",                 3953305,
    "Sardegna",               1611621,
    "Sicilia",                4875290,
    "Toscana",                3692555,
    "Umbria",                  870165,
    "Valle d'Aosta",           125034,
    "Veneto",                 4879133,
    "Italia",                59641488
)

# --- aux functions ---
# x is a data.frame with dates in first position
get_cmp <- function(x) {
    x[x < 0] <- NA
    k <- dim(x)[2] - 1
    n <- dim(x)[1]
    levels <- matrix(0, n, k, dimnames = list(NULL, names(x)[-1]))
    rates  <- matrix(0, n, k, dimnames = list(NULL, names(x)[-1]))
    for (i in 1:k) {
        y <- as.numeric(log(1 + x[[1+i]]))
        mod <- SSModel(y ~ SSMtrend(2, list(0, NA)) +
                           SSMseasonal(7, NA),
                       H = NA)
        vy <- var(diff(y), na.rm = TRUE)
        fit <- fitSSM(mod, log(c(vy / 100, vy / 100, vy / 100)))
        kfs <- KFS(fit$model, smoothing = "state")
        levels[, i] <- as.numeric(exp(kfs$alphahat[, 1] + 0.5 * kfs$V[1, 1, ])) - 1
        levels[levels < 0] <- 0
        rates[, i]  <- as.numeric(kfs$alphahat[, 2])
    }
    list(livello = as_tibble(data = x[[1]], levels),
         crescita = as_tibble(data = x[[1]], rates*100))
}
# ---------------------


dt <- jsonlite::fromJSON(txt = "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-json/dpc-covid19-ita-regioni.json")
regioni <- unique(dt$denominazione_regione)
start_date <- as.Date("2020-09-15")

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),

    # Application title
    titlePanel("COVID-19 Monitor Italia"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput("regioni",
                        "Selezionare regioni:",
                        choices = pop$Regione,
                        selected = "Italia"),
            helpText("Elaborazioni di Matteo Pelagatti sui dati della
                   da Protezione Civile pubblicati su",
                     a("https://github.com/pcm-dpc/COVID-19",
                       href = "https://github.com/pcm-dpc/COVID-19"),
                     ". Si consulti",
                     a("https://matteopelagatti.github.io/Covid-19/rapporto_italia_short",
                       href = "https://matteopelagatti.github.io/Covid-19/rapporto_italia_short"),
                     "per maggiori dettagli sul modello usato.")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Nuovi casi",
                    dygraphOutput("livello_casi"),
                    dygraphOutput("crescita_casi")
                ),
                tabPanel("Decessi",
                         dygraphOutput("livello_decessi"),
                         dygraphOutput("crescita_decessi")
                ),
                tabPanel("Tabelle",
                         textOutput("tab_testo"),
                         DTOutput("tabella")
                )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    # ---- casi    
    dati_casi <- reactive({
        casi <- dt %>%
            transmute(data = as.Date(strtrim(dt$data, 10)),
                      regione = denominazione_regione,
                      positivi = nuovi_positivi) %>%
            pivot_wider(names_from = regione, values_from = positivi) %>%
            arrange(data) %>%
            mutate(Italia = rowSums(.[, -1])) %>%
            filter(data >= start_date)
        
        casi_pop <- casi
        for (i in names(casi)[-1]) {
            casi_pop[, i] <- casi[, i] / pop$Popolazione[pop$Regione == i] * 100000
        }
        
        cmp_casi <- get_cmp(casi_pop)
        list(livello = xts(cmp_casi$livello, casi$data),
             crescita = xts(cmp_casi$crescita, casi$data))
    })
    
    output$livello_casi <- renderDygraph({
        dygraph(dati_casi()$livello[, input$regioni],
                main = "Nuovi casi ogni 100.000 abitanti")
    })
    output$crescita_casi <- renderDygraph({
        dygraph(dati_casi()$crescita[, input$regioni],
                main = "Tasso di crescita giornaliero dei nuovi casi (%)") %>%
            dyLimit(limit = 0, strokePattern = "dashed")
    })
    
    # ---- decessi    
    dati_decessi <- reactive({
        decessi <- dt %>%
            transmute(data = as.Date(strtrim(dt$data, 10)),
                      regione = denominazione_regione,
                      decessi = deceduti) %>%
            pivot_wider(names_from = regione, values_from = decessi) %>%
            arrange(data) %>%
            mutate(Italia = rowSums(.[, -1])) %>%
            mutate_at(vars(Abruzzo:Italia), function(x) x - lag(x)) %>%
            filter(data >= start_date)
        
        decessi_pop <- decessi
        for (i in names(decessi)[-1]) {
            decessi_pop[, i] <- decessi[, i] / pop$Popolazione[pop$Regione == i] * 1000000
        }
        
        cmp_decessi <- get_cmp(decessi_pop)
        list(livello  = xts(cmp_decessi$livello, decessi$data),
             crescita = xts(cmp_decessi$crescita, decessi$data))
    })
    
    output$livello_decessi <- renderDygraph({
        dygraph(dati_decessi()$livello[, input$regioni],
                main = "Decessi ogni 1.000.000 abitanti")
    })
    output$crescita_decessi <- renderDygraph({
        dygraph(dati_decessi()$crescita[, input$regioni],
                main = "Tasso di crescita giornaliero dei decessi (%)")
    })
    
    # ---- tabella
    output$tab_testo <- renderText({
        paste("Stime per l'ultimo dato disponbile:", max(as.Date(dt$data)))
    })

    output$tabella <- renderDT({
        tibble(Regione = pop$Regione,
               `Casi per 100k abitanti` = dati_casi()$livello %>% last() %>% as.numeric(),
               `Crescita casi (%)` = dati_casi()$crescita %>% last() %>% as.numeric(),
               `Decessi per 1M abitanti` = dati_decessi()$livello %>% last() %>% as.numeric(),
               `Crescita decessi (%)` = dati_decessi()$crescita %>% last() %>% as.numeric()) %>%
            datatable(rownames = FALSE, options = list(pageLength = 22, dom = "t")) %>%
            formatStyle(columns = c(3, 5),
                        color = styleInterval(0, c("black", "red"))) %>%
            formatRound(columns = 2:5, digits = 1)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
