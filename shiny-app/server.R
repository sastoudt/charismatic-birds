#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

if (!"scroller" %in% rownames(installed.packages())) {
    devtools::install_github("lgnbhl/scroller")
}
library(scroller)
source('loadShinyData.R')
library(shiny)
library(plotly)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    output$text <- renderPrint("
This tool is designed to visualize the results published in 'Identifying charismatic bird species and traits with community science observations' by Stoudt, Goldstein, and de Valpine.
    ")
    
    getPlotData <- function(input, output, plotCat, selectSigs = TRUE) {
        
        if (plotCat == "Traits") {
            highlights <- input$highlightSpecies_Traits
            orders <- input$highlightOrder_Traits
            hArg <- input$highlightArg_Traits
        } else if (plotCat == "Browse") {
            highlights <- input$highlightSpecies_Browse
            orders <- input$highlightOrder_Browse
            hArg <- input$highlightArg_Browse
        } else if (plotCat == "Counts") {
            highlights <- c()
            orders <- input$highlightOrder_Counts
            if (length(orders) > 0) {
                hArg <- "Order"
            } else {
                hArg <- "None"
            }
        }
        
        allData <- allData %>% 
            filter(abs(difference) <= 10) %>%
            mutate(ub95CI = pmin(ub95CI, 10), lb95CI = pmax(lb95CI, -10)) %>% 
            mutate(selected = 
                       switch(hArg,
                              None = {
                                  if (selectSigs) {
                                      sigdir != "Non-significant"
                                  } else {
                                      TRUE
                                  }
                              },
                              Significant = {
                                  sigdir != "Non-significant"
                              },
                              Species = {
                                  common_name %in% highlights
                              },
                              Order = {
                                  order %in% orders
                              }))
        return(allData)
    }
    
    output$traitPlot <- renderPlotly({
        allData <- getPlotData(input, output, plotCat = "Traits")
        
        switch(input$traitChoice,
               mass = {
                   massplot <- allData %>%                       
                       mutate(
                           scaled_mass = scale(log(mass)),
                       ) %>% 
                       # filter(diff_SE <= 5) %>%
                       ggplot(aes(scaled_mass, difference, col = sigdir, alpha = selected,
                                  text=sprintf("Species: %s<br>Difference: %s<br>SE: %s<br>%s", 
                                               common_name, 
                                               round(difference, 3), 
                                               round(diff_SE, 3), 
                                               sigdir))) + 
                       geom_point() +
                       geom_hline(yintercept = 0) +
                       geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI)) +
                       scale_color_manual("", values = c("#888888", poscol, negcol)) +
                       scale_alpha_manual("", values = c(0.22, 1), guide = F) +
                       theme_minimal() +
                       theme(panel.grid = element_blank()) +
                       guides(alpha=FALSE) +
                       # ylim(c(-10, 10)) +
                       xlab("Scaled log mass") +
                       ylab("Overreporting index") %>% 
                       add_trace()
                   
                   ggplotly(massplot, tooltip = "text")
               },
               color = {
                   colorplot <- allData %>%                         
                       mutate(
                           scaled_color = scale(max.color.contrast),
                       ) %>% 
                       # filter(diff_SE <= 5) %>%
                       ggplot(aes(scaled_color, difference, col = sigdir, alpha = selected,
                                  text=sprintf("Species: %s<br>Difference: %s<br>SE: %s<br>%s", 
                                               common_name, 
                                               round(difference, 3), 
                                               round(diff_SE, 3), 
                                               sigdir))) + 
                       geom_point() +
                       geom_hline(yintercept = 0) +
                       geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = selected)) +
                       scale_color_manual("", values = c("#888888", poscol, negcol)) +
                       scale_alpha_manual("", values = c(0.22, 1)) +
                       theme_minimal() +
                       theme(panel.grid = element_blank()) +
                       # ylim(c(-10, 10)) +
                       xlab("Scaled color contrast") +
                       ylab("Overreporting index")
                   
                   ggplotly(colorplot, tooltip = "text")
               },
               prevalence =  {
                   prevplot <- allData %>% 
                       mutate(
                           scaled_rate = scale(log(rate))
                       ) %>% 
                       # filter(diff_SE <= 5) %>%
                       ggplot(aes(scaled_rate, difference, col = sigdir, alpha = selected,
                                  text=sprintf("Species: %s<br>Difference: %s<br>SE: %s<br>%s", 
                                               common_name, 
                                               round(difference, 3), 
                                               round(diff_SE, 3), 
                                               sigdir))) + 
                       geom_point() +
                       geom_hline(yintercept = 0) +
                       geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = selected)) +
                       scale_color_manual("", values = c("#888888", poscol, negcol)) +
                       scale_alpha_manual("", values = c(0.22, 1)) +
                       theme_minimal() +
                       theme(panel.grid = element_blank()) +
                       # ylim(c(-10, 10)) +
                       xlab("Scaled log reporting rate") +
                       ylab("Overreporting index")
                   
                   ggplotly(prevplot, tooltip = "text")
               },
               range = {
                   rangeplot <- allData %>% 
                       # left_join(species_rarity) %>% 
                       mutate(
                           scaled_nhex = scale(log(nhex))
                       ) %>% 
                       ggplot(aes(scaled_nhex, difference, col = sigdir, alpha = selected,
                                  text=sprintf("Species: %s<br>Difference: %s<br>SE: %s<br>%s", 
                                               common_name, 
                                               round(difference, 3), 
                                               round(diff_SE, 3), 
                                               sigdir))) + 
                       geom_point() +
                       geom_hline(yintercept = 0) +
                       geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = selected)) +
                       scale_color_manual("", values = c("Non-significant" = "#888888", 
                                                         "Overreported" = poscol, 
                                                         "Underreported" = negcol)) +
                       scale_alpha_manual("", values = c(0.22, 1)) +
                       theme_minimal() +
                       theme(panel.grid = element_blank()) +
                       xlab("Scaled log range size") +
                       ylab("Overreporting index")
                   
                   ggplotly(rangeplot, tooltip = "text")
               }
        )

        # generate bins based on input$bins from ui.R
        # x    <- allData$difference
        # bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins

    })

    
    output$browsePlot <- renderPlotly({
        allData <- getPlotData(input, output, plotCat = "Browse", selectSigs = F)
        
        if (all(!allData$selected)) {
            allData$selected <- T
        }
        
        
        allData <- allData %>%
            filter(selected) %>% 
            mutate(common_name = gsub("_", " ", common_name))
        
        if (input$browseSortBy == "Index") {
            allData <- allData %>% arrange(difference)
        } else if (input$browseSortBy == "Lower bound") {
            allData <- allData %>% arrange(lb95CI)
        } else if (input$browseSortBy == "Upper bound") {
            allData <- allData %>% arrange(-ub95CI)
        } else if (input$browseSortBy == "Alphabetical") {
            allData <- allData %>% arrange(common_name)
        }
        
        browsePlot <- allData %>% 
            mutate(common_name = factor(common_name, levels = common_name)) %>%
            ggplot(aes(common_name, difference, col = sigdir,
                       text=sprintf("Species: %s<br>Difference: %s<br>SE: %s<br>%s", 
                                    common_name, 
                                    round(difference, 3), 
                                    round(diff_SE, 3), 
                                    sigdir))) +
            geom_linerange(aes(ymin = lb95CI, ymax = ub95CI)) +
            geom_hline(yintercept = 0) +
            geom_point(cex = 3) +
            coord_flip() +
            xlab("") +
            ylab("") +
            # ylab("Overreporting index (95% CI)") +
            theme_minimal() +
            scale_color_manual("", values = c("Non-significant" = "#888888", 
                                              "Overreported" = poscol, 
                                              "Underreported" = negcol)) +
            scale_x_discrete(position = "bottom") +
            scale_y_continuous(position = "right") +
            theme(panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank())
        ggplotly(browsePlot, tooltip = "text")
    })
    
    
    
    
    output$summaryPlot <- renderPlotly({
        allData <- getPlotData(input, output, plotCat = "Counts", selectSigs = F)
        
        if (all(!allData$selected)) {
            allData$selected <- T
        }

        if (length(input$highlightOrder_Counts) > 1) {
            title <- NULL
        } else if (length(input$highlightOrder_Counts) == 1) {
            title <- paste0("Order ", input$highlightOrder_Counts)
        } else  {
            title <- "All species"
        }
        
        allData <- allData %>%
            filter(selected) %>% 
            mutate(
                sigdir = factor(sigdir, levels = c(
                    "Underreported", "Non-significant", "Overreported"
                ))
            )
        
        if (length(input$highlightOrder_Counts) > 0) {
            plotData <- allData %>% 
                count(sigdir, order) %>%
                mutate(hovertext = sprintf("%s<br>n = %s<br>", sigdir, n))
            
            for (i in 1:nrow(plotData)) {
                orderSpecs <- allData %>% 
                    filter(order == plotData$order[i]) %>% 
                    filter(sigdir == plotData$sigdir[i])
                
                if (nrow(orderSpecs) > 20) {
                    plotData$hovertext[i] <- paste0(
                        plotData$hovertext[i],
                        paste(c(sort(orderSpecs$common_name)[1:20], "..."), 
                            collapse = "<br>")
                    )
                } else {
                    plotData$hovertext[i] <- paste0(
                        plotData$hovertext[i],
                        paste(sort(orderSpecs$common_name), 
                              collapse = "<br>")
                    )
                }
            }
        } else {
            plotData <- allData %>% 
                count(sigdir) %>% 
                mutate(hovertext = sprintf("%s<br>n = %s", sigdir, n))
        }
            
        
            
        summaryPlot <- plotData %>% 
            ggplot(aes(sigdir, n,
                       fill = sigdir, 
                       text = hovertext)) +
            xlab("") + ylab("") +
            geom_col() +
            scale_fill_manual(values = c("Non-significant" = "#888888", 
                                         "Overreported" = poscol, 
                                         "Underreported" = negcol)) +
            # geom_text(data = sigcounts, aes(sigdir, n, label = lab), vjust = "bottom", nudge_y = 5) +
            # theme_minimal() +
            ggtitle(title)
        
        if (length(input$highlightOrder_Counts) > 1) {
            summaryPlot <- summaryPlot + facet_wrap(~order)
        } else {
            summaryPlot <- summaryPlot + theme_minimal()
        }
        summaryPlot <- summaryPlot + 
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  # axis.text.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank())
        
        ggplotly(summaryPlot, tooltip = "text")
    })
    
    # output$mapPlot <- renderPlotly({
    #     
    #     if (length(input$map_species) == 0) {
    #         mapPlot <- ggplot()
    #     } else {
    #         this_map_data <- map_data %>% 
    #             filter(name_clean == input$map_species, success > 0) %>% 
    #             mutate(dataset = ifelse(dataset == "ebird_n", "eBird", "iNaturalist"))
    #         
    #         mapPlot <- this_map_data %>% 
    #             ggplot(aes(lon, lat, col = log(success/total),
    #                        text = sprintf("Rate = %s", round(success/total, digits = 3)))) +
    #             geom_point() +
    #             theme_minimal() +
    #             theme(axis.text = element_blank()) +
    #             xlab("Longitude") + ylab("Latitude") +
    #             labs(col = "Log reporting rate") +
    #             scale_color_gradient(low = "#2c7fb8", high = "#edf8b1") +
    #             facet_wrap(~dataset, ncol = 1)
    #     }
    #     
    #     ggplotly(mapPlot, tooltip = "text")
    # })
    
    
    
    output$mapplot.ui <- renderUI({
        plotlyOutput("mapPlot", height = 800)
    })
    
    plotHeight <- function() {
        dat <- getPlotData(input, output, plot = "Browse", selectSigs = F)
        numSpec <- sum(dat$selected)
        if (numSpec == 0) numSpec <- nrow(dat)
        return(100 + 20 * numSpec)
    }
    
    output$browseplot.ui <- renderUI({
        plotlyOutput("browsePlot", height = plotHeight())
    })
})
