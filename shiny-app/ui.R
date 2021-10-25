#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source('loadShinyData.R')
library(shiny)
library(plotly)
library(tidyverse)

fileName <- "html_about.html"
html_str <- readChar(fileName, file.info(fileName)$size)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    tabsetPanel(
        tabPanel("Trait plots", fluid = TRUE,
                 sidebarPanel(
                     h3("Overreporting indices by traits"),
                     p("Choose a trait to visualize how overreporting varies 
                       with it. Mouse over species indices to learn more, or 
                       select species or orders to focus on. Note that 
                       confidence intervals are symmetric but truncated at (-10, 10)."),
                     radioButtons(inputId = "traitChoice", 
                                  label = "Trait", 
                                  choiceNames = c(
                                      "Log mass", "Color", 
                                      "Log prevalence", "Log range size"),
                                  choiceValues = c(
                                      "mass", "color", "prevalence", "range"
                                  ),
                                  inline = FALSE),
                     
                     selectInput(inputId = "highlightSpecies_Traits", 
                                 label = "Select species", 
                                 choices = sort(allData$common_name),
                                 multiple = TRUE),
                     selectInput(inputId = "highlightOrder_Traits", 
                                 label = "Select order", 
                                 choices = sort(unique(allData$order)),
                                 multiple = TRUE),
                     radioButtons(inputId = "highlightArg_Traits", 
                                  label = "Highlight selection",
                                  choices = c("None", "Species", "Order"))
                     
                 ),
            mainPanel(plotlyOutput("traitPlot"))
        ),
        tabPanel("Browse species", fluid = TRUE,
                 
                 sidebarPanel(
                     h3("Overreporting indices by species"),
                     p("Mouse over species indices to learn more, or select 
                       species or orders to focus on. Note that confidence 
                       intervals are symmetric but truncated at (-10, 10)."),
                     style = "position:fixed;width:inherit;",
                     selectInput(inputId = "highlightSpecies_Browse", 
                                 label = "Select species", 
                                 choices = sort(allData$common_name),
                                 multiple = TRUE),
                     selectInput(inputId = "highlightOrder_Browse", 
                                 label = "Select order", 
                                 choices = sort(unique(allData$order)),
                                 multiple = TRUE),
                     radioButtons(inputId = "highlightArg_Browse", 
                                  label = "Highlight selection",
                                  choices = c("None", "Species", "Significant", "Order")),
                     radioButtons(inputId = "browseSortBy", 
                                  label = "Sort by",
                                  choices = c("Index", "Lower bound", 
                                              "Upper bound", "Alphabetical"))
                     
                 ),
                 mainPanel(uiOutput("browseplot.ui"))
        ),
        tabPanel("Summaries", fluid = TRUE,
                 sidebarPanel(
                     h3("Over- and underreporting summaries"),
                     p("View counts of over or underreported species across all orders
                       or grouped by taxonomic order."),
                     selectInput(inputId = "highlightOrder_Counts", 
                                 label = "Select order", 
                                 choices = sort(unique(allData$order)),
                                 multiple = TRUE),
                 ),
                 mainPanel(plotlyOutput("summaryPlot"))
        ),
        # tabPanel("Data maps", fluid = TRUE,
        #          sidebarPanel(
        #              h3("Visualize reporting patterns in space"),
        #              p("View distributions of reporting rates of each species in space.
        #                In our analysis, we fit spatial splines to these surfaces 
        #                and calculated the typical difference between them."),
        #              selectInput(inputId = "map_species", 
        #                          label = "Select species", 
        #                          choices = sort(unique(allData$common_name)),
        #                          multiple = FALSE),
        #          ),
        #          mainPanel(uiOutput("mapplot.ui"))
        # ),
        tabPanel("About", fluid = TRUE,
                 mainPanel(HTML(html_str))
        )
    )
))