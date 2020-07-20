library(shiny)
library(feather)

# Read in data - stored in feather format
mas <- read_feather("mas.feather")

library(data.table)

# Convert dataframe to data.table
# Using set to avoid copying of dataframe
# Setting the reference key (for joins/subsets) to Bradi id column
setDT(mas, key="CycID")


library(ggplot2)
library(viridisLite)

# Create elements of user interface
ui <- fluidPage(
  
  # title of page
  title="Bd21 time course data - Transcript Mapped",
  # header at top of page
  h1("Bd21 Time Course Data - Transcript Mapped"),	
  
  br(),
  
  #sidebar element
  sidebarLayout(sidebarPanel (width = 3,
  fluidRow(
          # selectize input allows the selection of ids as the user types
          # Initally set choices to null - updated in the server section to speed up page load
           selectizeInput(inputId = "bradi", label="Transcript identifier:", choices = NULL, selected=NULL, multiple=TRUE, options = list( create=TRUE,selectOnTab = TRUE, closeAfterSelect = TRUE)) ,

           br(),
           
           # file upload button - specifies allowed file types
           fileInput("newFile", "Upload list of ids - 1 per line.", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv", ".txt")),
           br(),
           br(),
           hr(),
           
           # slider input to adjust plot width
           sliderInput("plotWidth", "Width of plot (in pixels)", min=200, max=4000, value=1200, step=100),
           # slider input to adjust plot height
           sliderInput('plotHeight', 'Height of plot (in pixels)', min = 200, max = 4000, value = 600, step=100),
           br(),
           br(),
           
           # Section to fill in experimental details in the sidebar
           # eventually planning to make this section update as each tab is clicked
           h6("We explored two environmental cues: photocycle (L:light, and D:dark, 12h:12h) and thermocycle (H:hot 28C, and C:cold 12C, 12h:12h). Accession Bd21 seeds were grown on 100 mL 1x MS + 1% Agar in Magenta boxes. All seeds were germinated in dark for 3 days then transferred to LDHC for 10 days. Seedlings were then released into one of four conditions: LDHC, LDHH, LLHC or LLHH. 12 hours after placing the samples in each condition the currently emerging leaf was sampled and pooled from at least 3 separate plants. 13 additional time points were collected following this first sample (total 14 time points), one every 3.5 hours. RNA-seq was performed on all samples to quantify transcript abundace.")
    )

    ),
  # end of side bar element, now working on the main center panel
  mainPanel(
  hr(),
  fluidRow(
    
           # using a verbatim text output element to have an auto updating 
           # text section that displays the defline from phytozome annotation file
           h6(verbatimTextOutput("defline"))
           
  ),
  hr(),
  
  # Creating a tabbed section underneath the verbatim text output
  tabsetPanel(type="tabs",
               # First tab panel displays the standard plot
                                tabPanel("Plot", value=1, br(), 
                                         conditionalPanel( condition = "input.bradi != ''",
                                         h4("Values shown are DESeq2 normalized read counts."),                      
                                         br(),
                      # second tabbed section within the original tab - made it easier to deal with multiple plots
                                         tabsetPanel(type = "pills" ,
                      # first tab panel within nested tabbed section
                                                     tabPanel("All", align="center", br(), br(), plotOutput("graph", height = "auto")),
                       # second tab panel within nested tabbed section
                       tabPanel("By Condition",  br(), br(),plotOutput("faceted", height = "auto")),
                       # third tab panel within nested tabbed section
                       tabPanel("By ID",  br(), br(),plotOutput("face2", height = "auto"))
                      )
                                )
                      ),
              # end of first nested tabbed section - new tab panel next to "Plot"
              tabPanel("Scaled", value=2,
                       h4("Normalized read counts scaled by z-score transform."),
                       # Nested tabbed section under the scaled panel
                       tabsetPanel(type = "pills" ,
                       # first nested tab panel
                                   tabPanel("All",plotOutput("sgraph", height = "auto" )),
           
                        # second nested tab panel
                                   tabPanel("By Condition",plotOutput("sfaceted", height = "auto" )),
           
                        #third nested tab panel
                                   tabPanel("By ID",plotOutput("sface2" , height = "auto"))
                       )
              ),
              
              
              
              
              # end of second nested tab section 
              # Third main tab panel
              tabPanel("MetaCycle Output",
                       h3("MetaCycle Output:"),
                       br(),
                       tableOutput("mytable")
              ),
              # Fourth main tab panel
              tabPanel("Variance",
                       h3("Coefficient of Variance:"),
                       br(),
                       h6("(Standard deviation / mean) reported as percentage."),
                       tableOutput("variance")
              ),
              
             
              # Fifth main tab panel
              tabPanel("Expression Values", 
                       h3("Expression Values:"), 
                       br(),
                       h4("Values shown are DESeq2 normalized read counts."),
                       br(),
                       tableOutput("expressionTable")
              )
  , id = "tabselected")
  
  
  
  )
  )
)
	
# End of user interface elements - functional aspects go in server section

server <- function(input, output, session) {
  
  # Observe function to be running basically constantly
	observe({
	  # parseQueryString will parse the url to accept passing a bradi id in the url to the app
	  # https://hazenlab.shinyapps.io/shinyapp/?id=Bradi3g33170.2
	  query <- parseQueryString(session$clientData$url_search)
	  # If query isn't null then update the selectize input - changing the choices, and the selected
	  # to be the query
	  if (!is.null(query[['id']])) {
	    updateSelectizeInput(session, "bradi", choices = query[['id']], selected= query[['id']], server = TRUE, options = list(create=TRUE))
	  }
	  # Check to see if the file upload section is empty or not, if a file is loaded then try to read it in as a dataframe, then update the selectize input to reflect the selections
	  if ( !is.null(input$newFile))
	  {
	    try(ids <- fread(input$newFile$datapath, header=FALSE))
	    updateSelectizeInput(session, "bradi",  choices=mas$CycID, selected= ids$V1, server = TRUE, options = list(create=TRUE))
	    
	  }
	 
	})
	
# After graphical elements have loaded - update the selectize input to use the data.table as all the choices for input
	  updateSelectizeInput(session, "bradi", choices = mas$CycID, server = TRUE ,options = list(create=TRUE))
	  

	  
	  
# Once a search term has been input - the initial tidying of the data
	searchTerm <- reactive({
	  #prevents the function from running if search is empty or too short (will not update as person is typing)
	  if (is.null(input$bradi)) return()
	  if (nchar(input$bradi) < 11) return()
	  # trim the white space from the id
	  trimmed <- trimws(input$bradi, which="both")
	  trimmed <- gsub(' ', '', trimmed)
	  # subset the data.table into a temporary data.table that only contains selected ids
	  temp <- mas[paste0(trimmed)]
	  
	  # If someone enters a gene id instead of a transcript id - do a search for matching ids 
	  if ( is.na(temp$Condition) == TRUE ) {
	    temp <- mas[ CycID %ilike% paste0(trimmed)]
	  }
	  
	  # convert from wide data.table to narrow for plotting
	  temp <- melt.data.table(temp, measure.vars = patterns("T_") )
	  # get rid of characters in hour data - they come from column titles
	  temp[, variable := gsub("T_", "", variable)]
	  # convert hour data to numeric
	  temp[, variable := as.numeric(variable)]
	  # create new column for scaled data
	  temp[, scaled := scale(value), by=.(CycID, Condition)]
	  # return new data.table - it will be stored as "searchTerm()"
	  temp
	  })


# Section to render the output for the defline verbatim text output
	output$defline <- renderText({
	  if ( is.null(searchTerm()$Condition) ) return() 
	  # Return unique deline/arabidopsis id matches from searchTerm data.table
		unique(paste(searchTerm()$CycID, searchTerm()$`Best-hit-arabi-name`, searchTerm()$`arabi-defline`, "\n", sep="\t"))
	  
	 })
	
# Section to render metacycle table	
	output$mytable <- renderTable({
	  if ( is.null(searchTerm()$Condition)) return() 
	  # Returning unique CycID Condition and all outputs from metacycle
	  temp <- unique(searchTerm()[ , names(searchTerm()) %like% "CycID|Condition|meta", with=FALSE])
	  # Formatting meta2d_pvalue column to fewer digits for readability
	  temp[, meta2d_pvalue := format(as.numeric(meta2d_pvalue), digits=3, scientific = TRUE)]
	  temp
	  
	})	

# Section to render variance table		
	output$variance <- renderTable({
	  if (is.null(searchTerm()$Condition)) return()
	  # going from narrow data.table to wide while simultaneously calculating coefficient of variance
	  	  dcast.data.table(searchTerm()[, .(cv = (sd(value) / mean(value) * 100)), by=.(CycID, Condition)], CycID ~ Condition, value.var = "cv")
	  
	  
	})
	
# Section to render expression table	
	output$expressionTable <- renderTable({
	  if (is.null(searchTerm()$Condition)) return() 
	  # going from narrow data to wide and selecting which columns I want to return
	  dcast.data.table(searchTerm(), CycID + Condition ~ variable, value.var = "value") 

	  })
	
	
# the following sections are all to render various plots - each one is generated individually based on the user interface element that called it
	# this could be improved to be more concise
	
	output$graph <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return()
	    # fairly standard ggplot call to generate the plot
	    # notable difference is the syntax for searchTerm() - your data.table is stored as a function, it looks different but functions in the same way
		  ggplot(searchTerm(), aes(x=variable, y=value, group=interaction(CycID, Condition), color=Condition, shape=CycID)) + 
	    theme_bw() + 
	    ylab("Normalized counts") + 
	    xlab("Time (h)") + 
	    theme(panel.grid=element_blank()) + 
	    scale_color_viridis_d() + 
	    geom_line() + 
	    geom_point() + 
	    annotate("rect", xmin=12, xmax=24, ymin=0, ymax=Inf , fill="gray", color=NA, alpha=0.30) + 
	    annotate("rect", xmin=36, xmax=48, ymin=0, ymax=Inf, fill="gray", color=NA, alpha=0.30) + 
	    scale_x_continuous(breaks= c(0,12,24,36)) + 
	    theme() +
	    labs(color="", shape="")
	  # another significant point - in order to refer to the slider generated height/width values
		# you need to encase them in a exprToFunction call
		}, res=300, height = exprToFunction(input$plotHeight), width = exprToFunction(input$plotWidth))
	
	# The remainder are just specific plots - subtlety adjusted to meet the criteria for the tab
	output$faceted <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return() 
	  
	  ggplot(searchTerm(), aes(x=variable, y=value, group=interaction(CycID,Condition) )) + 
	    theme_bw() + 
	    ylab("Normalized counts") + 
	    xlab("Time (h)") + 
	    theme(panel.grid=element_blank()) + 
	    scale_color_viridis_d() + 
	    geom_line(aes(color=CycID)) + 
	    geom_point(aes(shape=CycID, color=CycID)) + 
	    annotate("rect", xmin=12, xmax=24, ymin=0, ymax=Inf , fill="gray", color=NA, alpha=0.30) + 
	    annotate("rect", xmin=36, xmax=48, ymin=0, ymax=Inf, fill="gray", color=NA, alpha=0.30) + 
	    scale_x_continuous(breaks= c(0,12,24,36)) + 
	    facet_wrap(~Condition, scales="free_y", ncol=1) +
	    theme(legend.position = "bottom") +
	    labs(color="", shape="")
	  
	}, res=300 , height = exprToFunction(4 * input$plotHeight), width = exprToFunction(input$plotWidth))
		
	output$sgraph <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return() 
	  
	  ggplot(searchTerm(), aes(x=variable, y=scaled, group=interaction(CycID, Condition), color=Condition, shape=CycID)) + 
	    theme_bw() + 
	    ylab("Scaled expression") + 
	    xlab("Time (h)") + 
	    theme(panel.grid=element_blank()) + 
	    scale_color_viridis_d() + 
	    geom_line() + 
	    geom_point() + 
	    annotate("rect", xmin=12, xmax=24, ymin=-Inf, ymax=Inf , fill="gray", color=NA, alpha=0.30) + 
	    annotate("rect", xmin=36, xmax=48, ymin=-Inf, ymax=Inf, fill="gray", color=NA, alpha=0.30) + 
	    scale_x_continuous(breaks= c(0,12,24,36)) + 
	    theme() +
	    labs(color="", shape="")
	  
	}, res=300, height = exprToFunction(input$plotHeight), width = exprToFunction(input$plotWidth))
	
	output$sfaceted <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return() 
	  ggplot(searchTerm(), aes(x=variable, y=scaled, group=interaction(CycID,Condition) )) + 
	    theme_bw() + 
	    ylab("Scaled expression") + 
	    xlab("Time (h)") + 
	    theme(panel.grid=element_blank()) + 
	    scale_color_viridis_d() + 
	    geom_line(aes(color=CycID)) + 
	    geom_point(aes(shape=CycID, color=CycID)) + 
	    annotate("rect", xmin=12, xmax=24, ymin=-Inf, ymax=Inf , fill="gray", color=NA, alpha=0.30) + 
	    annotate("rect", xmin=36, xmax=48, ymin=-Inf, ymax=Inf, fill="gray", color=NA, alpha=0.30) + 
	    scale_x_continuous(breaks= c(0,12,24,36)) + 
	    facet_wrap(~Condition, scales="free_y", ncol=1) +
	    theme(legend.position = "bottom") +
	    labs(color="", shape="")
	  
	}, res=300 , height = exprToFunction(4 * input$plotHeight), width = exprToFunction(input$plotWidth))

	output$face2 <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return() 
	  
	  ggplot(searchTerm(), aes(x=variable, y=value, group=interaction(CycID,Condition) )) + 
	    theme_bw() + 
	    ylab("Normalized counts") + 
	    xlab("Time (h)") + 
	    theme(panel.grid=element_blank()) + 
	    scale_color_viridis_d() + 
	    geom_line(aes(color=Condition)) + 
	    geom_point(aes(shape=Condition, color=Condition)) + 
	    annotate("rect", xmin=12, xmax=24, ymin=-Inf, ymax=Inf , fill="gray", color=NA, alpha=0.30) + 
	    annotate("rect", xmin=36, xmax=48, ymin=-Inf, ymax=Inf, fill="gray", color=NA, alpha=0.30) + 
	    scale_x_continuous(breaks= c(0,12,24,36)) + 
	    facet_wrap(~CycID, scales="free_y", ncol=1) +
	    theme(legend.position = "bottom") +
	    labs(color="", shape="")
	  
	}, res=300, height = exprToFunction(length(unique(searchTerm()$CycID)) * input$plotHeight), width = exprToFunction(input$plotWidth) )
	
	output$sface2 <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return() 

	  ggplot(searchTerm(), aes(x=variable, y=scaled, group=interaction(CycID,Condition) )) + 
	    theme_bw() + 
	    ylab("Scaled expression") + 
	    xlab("Time (h)") + 
	    theme(panel.grid=element_blank()) + 
	    scale_color_viridis_d() + 
	    geom_line(aes(color=Condition)) + 
	    geom_point(aes(shape=Condition, color=Condition)) +
	    annotate("rect", xmin=12, xmax=24, ymin=-Inf, ymax=Inf , fill="gray", color=NA, alpha=0.30) + 
	    annotate("rect", xmin=36, xmax=48, ymin=-Inf, ymax=Inf, fill="gray", color=NA, alpha=0.30) + 
	    scale_x_continuous(breaks= c(0,12,24,36)) + 
	    facet_wrap(~CycID, scales="free_y", ncol=1) +
	    theme(legend.position = "bottom") +
	    labs(color="", shape="")
	  
	}, res=300 , height = exprToFunction(length(unique(searchTerm()$CycID)) * input$plotHeight), width = exprToFunction(input$plotWidth))



	
	
}

# the final call to run the shinyApp
shinyApp(ui,server)
