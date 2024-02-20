#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Comparative Mus Musculus RNA/Protein Expression Over Time"),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "proteinoi", 
                label = "Protein(s) of Interest: (comma separated)", 
                value = "", width = NULL, placeholder = NULL
                ),
      checkboxGroupInput(inputId = "timepoint",
        label = "Time Points", 
        choices = c("E09", "E10", "E11", "E12", "E13", "E14", "E15", "E16"),
        selected = c("E09", "E10", "E11", "E12", "E13", "E14", "E15", "E16"),
        inline = FALSE,
        width = NULL,
        choiceNames = NULL,
        choiceValues = NULL
      ),
      selectInput("normalize", "Select Time Point to Normalize to", choices = c("E09", "E10", "E11", "E12", "E13", "E14", "E15", "E16")),
      actionButton(inputId = "update", label = "update")
    ),
    mainPanel(
      plotOutput(outputId = "scatterplot"),
      downloadButton("DownloadPlot", "Download Plot")
    )
  )
)

server <- function(input, output, session) {
  #Functions
  ## Get Metadata
  get_metadata<- reactive({
    proteinstouse<- input$proteinoi
    proteinstouse<-unlist(strsplit(proteinstouse, ","))
    #print(proteinstouse)
    return(proteinstouse)
  })
  
  get_timepoints<- reactive({
    timepointstouse <- input$timepoint
    #print(timepointstouse)
    return(timepointstouse)
  })
  
  get_normalized_timepoint<- reactive({
    normtimepoint<- input$normalize
    #print(normtimepoint)
    return(normtimepoint)
  })
  
  df2=NULL
  
  #When "update" is pressed:
  shiny::observeEvent(input$update,{
    #Update the proteins list
    proteinsupdate<- get_metadata()
    output$proteinstouse<- shiny::renderText({proteinsupdate})
    #Update the time points to include
    timepointsupdate<- get_timepoints()
    output$timepointstouse<- shiny::renderText({timepointsupdate})
    #Update the time point to normalize to
    normalizeupdate<- get_normalized_timepoint()
    output$normalizedtimepoint<- shiny::renderText({normalizeupdate})
    #Get the data for the genes we want and the time points we want
    #A "for loop" iterates through all x in y. In this example, it will perform the code inside the {} for every gene in the variable 'genes' that I defined above.
    for (p in {proteinsupdate}){
      #Check to see if protein is in the database
      if (p %in% rownames(rawvalues)){
        a<- getdata(p,{normalizeupdate})
        df2=rbind(df2,a)
      }

    }
    #Subset only the timepoints you want
    df2<- df2[df2$day %in% {timepointsupdate},]
    #Plot the scatterplot
    output$scatterplot <- renderPlot({
      p<- ggplot(df2, aes(x=day, y=Normalized, group=protein, color=protein)) + 
        geom_line() +
        geom_point()+
        geom_errorbar(aes(ymin=Normalized-sd, ymax=Normalized+sd), width=.2,
                      position=position_dodge(0.05))
      p=p+labs(title="Protein Expression", x="Day", y = "Avg. Abundance (Normalized)")+
        theme_classic()
      plot(p)
      observeEvent(input$update, print(as.numeric(input$update)))
    })
    p<- ggplot(df2, aes(x=day, y=Normalized, group=protein, color=protein)) + 
      geom_line() +
      geom_point()+
      geom_errorbar(aes(ymin=Normalized-sd, ymax=Normalized+sd), width=.2,
                    position=position_dodge(0.05))
    p=p+labs(title="Protein Expression", x="Day", y = "Avg. Abundance (Normalized)")+
      theme_classic()
    
    output$DownloadPlot = downloadHandler(
      file= paste({proteinsupdate}, {timepointsupdate}, ".png", sep=""),
      content = function(file){
        ggsave(p, filename = file)
      }
    )
  })
  
  #Load in the data frame
  masterdf <- read_excel("~/Documents/UNC Consulting/Protein_Expression_Project/1-s2.0-S1534580723001818-mmc2.xlsx", sheet = "Master Protein Table")
  genes<- masterdf$`Gene Symbol`
  rawvalues<- as.data.frame(masterdf[,29:52])
  rownames(rawvalues)<- make.names(genes, unique = TRUE)
  
  #Get a data frame for the proteins of interest
  getdata<- function(protein,timepoint){
    tp<- c("E09", "E10", "E11", "E12", "E13", "E14", "E15", "E16")
    normtp<- tp[tp %in% timepoint]
    protein<- protein[protein %in% rownames(rawvalues)]
    df<- as.data.frame(t(rawvalues[rownames(rawvalues)==protein,]))
    df$day<- c(rep("E09",3), rep("E10",3), rep("E11",3), rep("E12",3), rep("E13",3), rep("E14",3), rep("E15",3), rep("E16",3))
    df$protein<- protein
    colnames(df)[1]<-"Expression"
    #Get all of the average expressions for each time point
    avg<- c(mean(df$Expression[1:3]), mean(df$Expression[4:6]), mean(df$Expression[7:9]), mean(df$Expression[10:12]), mean(df$Expression[13:15]), mean(df$Expression[16:18]), mean(df$Expression[19:21]), mean(df$Expression[22:24]))
    names(avg)<- tp
    #Normalize to the time point we want
    df$Normalized<- df$Expression/(avg[names(avg) %in% normtp])
    
    data_summary <- function(data, varname, groupnames){
      require(plyr)
      summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum<-ddply(data, groupnames, .fun=summary_func,
                      varname)
      data_sum <- rename(data_sum, c("mean" = varname))
      return(data_sum)
    }
    df2 <- data_summary(df, varname="Normalized", 
                        groupnames=c("day", "protein"))
    head(df2)
    return(df2)
  }
  
}

# Run the application 
shinyApp(ui = ui, server = server)
