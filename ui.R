### Salmonella Shiny App
# 
# This is the user-interface definition of a Shiny web application.
#
###


library(shiny)
library(plotly)

library(ISOweek)
# library(lubridate)

# Define UI for application
shinyUI(
  fluidPage(
    # Application title
    # titlePanel("Salmonella Shiny App"),
    
    # Sidebar with a input controls 
    sidebarLayout(
      sidebarPanel(width = 2, 
                   fileInput('file', 'Choose CSV file:', accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                   # fileInput('file1', 'Choose data file:', accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv', ".xls")),
                   
                   radioButtons("datesType", "Type of dates:", c("Date of reception", "Date of isolation")),
                   
                   selectInput("year", "Year:", c()),
                   tags$div(
                     style = "margin-left:30px;",
                     selectInput("week", "ISO week:", c())
                   ),
                   
                   selectInput("region", "Region:", c()),
                   tags$div(
                     style = "margin-left:30px;",
                     selectInput("department", "Department:", c())
                   ),
                   
                   selectInput("serotype", "Serotype:", c()),
                   tags$div(
                     style = "margin-left:30px;",
                     selectInput("subtype", "Subtype:", c())
                   ),
                   
                   selectInput("ageGroup", "Age group:", c()),
                   
                   # actionButton("button", "Detect outbreaks", width = 130),
                   
                   checkboxInput("detectReg", "Regional alerts", FALSE),
                   checkboxInput("detectSer", "Serotype alerts", FALSE)
      ),
      
      mainPanel(
        tableOutput("contents"),
        
        tabsetPanel(
          tabPanel("Distributions",
                   
                   column(6, 
                          tags$div(
                            style = "margin-top:30px;margin-bottom:20px;",
                            plotlyOutput("temporalPlotly", height = 300)
                          ),
                          tags$div(
                            style = "margin-top:30px;margin-bottom:20px;",
                            plotlyOutput("serotypePlotly", height = 300)
                          )
                   ),
                   
                   column(6,
                          tags$div(
                            style = "margin-top:30px;margin-bottom:20px;",
                            plotlyOutput("regionalPlotly", height = 300)
                          ),
                          tags$div(
                            style = "margin-top:30px;margin-bottom:20px;",
                            plotlyOutput("dynamicPlotly", height = 300)
                          )
                   ),
                   
                   verbatimTextOutput("info"),
                   
                   column(6, div(dataTableOutput("table"), style = "font-size:80%"))
          ),
          
          tabPanel("Projections",
                   tags$div(
                     style = "margin-top:20px;margin-bottom:20px;",
                     plotlyOutput("projectionPlot")
                   ), 
                   sliderInput("highSeasonRange", "High season (weeks):",
                               min = 0, max = 52, value = c(29, 40)),
                   sliderInput("coeffForHighSeason",
                               "Percent of samples sequenced in high season:",
                               min = 0,
                               max = 1,
                               value = 0.5),
                   
                   sliderInput("coeffDefault",
                               "Percent of samples sequenced by default:",
                               min = 0,
                               max = 1,
                               value = 0.7,
                               step = 0.01)
          ),
          
          tabPanel("Case clusters", 
                   fluidRow(
                     splitLayout(cellWidths = c("30%", "70%"),
                                 tags$div(
                                   style = "margin-top:20px;margin-bottom:20px;",
                                   plotlyOutput("clusterHist")
                                 ),
                                 tags$div(
                                   style = "margin-top:20px;margin-bottom:20px;",
                                   plotlyOutput("clusterPlotly")
                                 )
                     )),
                   # plotOutput("clusterPlot"),
                   sliderInput("window", "Window (days):",
                               min = 0, max = 50, value = 30),
                   sliderInput("threshold", "Threshold (cases):",
                               min = 0, max = 30, value = 3)
          ),
          
          tabPanel("Monthly cases", 
                   tags$div(
                     style = "margin-top:20px;margin-bottom:20px;",
                     plotlyOutput("monthlyHist")
                   )
          )
        )
      )
      
    )
  )
)
