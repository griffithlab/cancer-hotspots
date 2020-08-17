library(shiny)
library(ggplot2)


ui <- fluidPage(
    #selectInput(inputId = "box", label = "Choose gene", choices = <genes>
    #), 
    plotOutput(outputId = "scatter")
)

server <- function(input, output) {
    #output$scatter <- renderPlot({ <function(data)> })
    #output$ui_output_name <- render*({ plot_function(input$ui_input_name) })
}

shinyApp(ui = ui, server = server)