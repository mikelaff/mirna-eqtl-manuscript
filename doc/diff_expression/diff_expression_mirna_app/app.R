#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(mikelaff)

#source(here("src/utils/lafferty_utils.R"))

# basic plot theme
plotTheme <- theme(axis.title = element_text(size = 16, face = "bold"),
                   axis.text = element_text(size = 14),
                   title = element_text(size = 18, face = "bold"),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 14),
                   plot.caption = element_text(size = 14, face = "plain", hjust = 0.5, margin = unit(c(5,0,0,0), "mm")),
                   panel.background = element_rect(fill = "white", linetype = "solid", color = "black", size = 1),
                   panel.grid = element_blank())

# INPUT FILES #########################################################################################################
# diff expression data frame, mirna
df.de.mirna.rds <- here("doc/diff_expression/rdata/20190519_diff_expression_df.rds")


# load data
df <- readRDS(df.de.mirna.rds)

# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Fetal Cortical Tissue miRNA Differential Expression"),

    # Copy the line below to make a text input box
    textInput(inputId = "mir",
              label = h3("miRNA"),
              value = "hsa-miR-124-3p"
              ),

    # Show a plot
    mainPanel(

        plotOutput(outputId = "MAplot_gzcp", height = 400, width = 700),
        plotOutput(outputId = "MAplot_gw", height = 400, width = 700)

    )

)

# Define server logic required to draw a histogram
server <- function(input, output) {

    datadata <- df


    output$MAplot_gzcp <- renderPlot({



        datadata %>%
            dplyr::filter(baseMean.gzcp > 0) %>%
            arrange(sig.gzcp) %>%
            ggplot(aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp)) +
            geom_point() +
            labs(x="Mean of Normalized Counts",
                 y="Shrunken Log2 Fold Change",
                 title="GZ/CP Differential Expression",
                 color=paste("p.adj < 0.1", sep="")) +
            scale_x_log10() +
            theme(legend.position = "bottom") +
            geom_hline(yintercept = 0, size=1, color="blue") +
            plotTheme +
            scale_color_manual(values=c("grey70", "blue")) +
            geom_point(data = dplyr::filter(datadata, grepl(input$mir, Name)),
                       mapping = aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp),
                       shape=21, size=4, stroke=2) +
            geom_label(data = dplyr::filter(datadata, grepl(input$mir, Name)),
                       mapping = aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp, label=Name),
                       hjust=1, vjust=1.5, size=5, show.legend = FALSE)

    })

    output$MAplot_gw <- renderPlot({



        datadata %>%
            dplyr::filter(baseMean.gw > 0) %>%
            arrange(sig.gw) %>%
            ggplot(aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw)) +
            geom_point() +
            labs(x="Mean of Normalized Counts",
                 y="Shrunken Log2 Fold Change",
                 title="Gest. Week Differential Expression",
                 color=paste("p.adj < 0.1", sep="")) +
            scale_x_log10() +
            theme(legend.position = "bottom") +
            geom_hline(yintercept = 0, size=1, color="blue") +
            plotTheme +
            scale_color_manual(values=c("grey70", "blue")) +
            geom_point(data = dplyr::filter(datadata, grepl(input$mir, Name)),
                       mapping = aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw),
                       shape=21, size=4, stroke=2) +
            geom_label(data = dplyr::filter(datadata, grepl(input$mir, Name)),
                       mapping = aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw, label=Name),
                       hjust=1, vjust=1.5, size=5, show.legend = FALSE)

    })
}

# Run the application
shinyApp(ui = ui, server = server)
