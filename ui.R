#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


library(shiny)
library(plotly)
library(shinyjs)

library(ISOweek)
# library(lubridate)

# Define UI for application
shinyUI(
  fluidPage(
    # Application title
    # titlePanel("Estimation of WGS capacity"),
    
    useShinyjs(),
    # code to reset plotlys event_data("plotly_click", source="A") to NULL -> executed upon action button click
    # note that "A" needs to be replaced with plotly source string if used
    extendShinyjs(text = "shinyjs.resetSerClick = function() { Shiny.onInputChange('.clientValue-plotly_click-serPlot', 'null'); }"),
    extendShinyjs(text = "shinyjs.resetRegClick = function() { Shiny.onInputChange('.clientValue-plotly_click-regPlot', 'null'); }"),
    
    
    extendShinyjs(text = "shinyjs.resetAgeClick = function() { Shiny.onInputChange('.clientValue-plotly_click-agePlot', 'null'); }"),
    extendShinyjs(text = "shinyjs.resetDepClick = function() { Shiny.onInputChange('.clientValue-plotly_click-depPlot', 'null'); }"),
    
    extendShinyjs(text = "shinyjs.resetWeekClick = function() { Shiny.onInputChange('.clientValue-plotly_click-mainPlot', 'null'); }"),
    
    # Sidebar with a slider input for number of bins 
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
                   
                   
                   # selectInput("curWeek", "Current week:",
                   #             c(NA, 0:52), selected = NA),
                   
                   # sliderInput("highSeasonRange", "High season (weeks):",
                   #             min = 0, max = 52, value = c(29, 40)),
                   
                   
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        # plotOutput("distPlot"),,
        
        tableOutput("contents"),
        
        tabsetPanel(
          # tabPanel("Total", plotOutput("distPlot")),
          # tabPanel("Summary", verbatimTextOutput("summary")),
          # tabPanel("Table", tableOutput("table"))
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
          # tabPanel("Distributions (old)",
          #          
          #          column(6, 
          #                 tags$div(
          #                   style = "margin-top:20px;margin-bottom:20px;",
          #                   plotlyOutput("mainPlotly", height = 300)
          #                 ),
          #                 tags$div(
          #                   style = "margin-top:30px;margin-bottom:20px;",
          #                   plotlyOutput("regionalPlotly", height = 150)
          #                 ),
          #                 tags$div(
          #                   style = "margin-top:30px;margin-bottom:20px;",
          #                   plotlyOutput("mapPlotly", height = 200)
          #                 )
          #          ),
          #          
          #          column(6, 
          #                 tags$div(
          #                   style = "margin-top:20px;margin-bottom:20px;",
          #                   plotlyOutput("mainRegPlotly", height = 300)
          #                 ),
          #                 tags$div(
          #                   style = "margin-top:30px;margin-bottom:20px;",
          #                   plotlyOutput("serotypePlotly", height = 150)
          #                 ),
          #                 tags$div(
          #                   style = "margin-top:30px;margin-bottom:20px;",
          #                   plotlyOutput("stPlotly", height = 200)
          #                 )
          #          )
          # ),
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
          # ,
          # tabPanel("Test", 
          #          # plotOutput("testPlot"),
          #          verbatimTextOutput("info")
          # )
        )
      )
      
    )
  )
)
