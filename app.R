library(shiny)

library(feather)




mas <- read_feather("mas.feather")

library(data.table)
setDT(mas, key="CycID")
#setkey(mas, CycID)

library(ggplot2)
library(viridisLite)

#firstlink <- "https://phytozome.jgi.doe.gov/pz/portal.html#!gene?search=1&detail=1&method=5010&searchText=transcriptid:"

ui <- fluidPage(
  
  title="Bd21 time course data - Transcript Mapped",
  h1("Bd21 Time Course Data - Transcript Mapped"),	
  
  br(),
  
  sidebarLayout(sidebarPanel (width = 3,
  fluidRow(
    
           selectizeInput(inputId = "bradi", label="Transcript identifier:", choices = NULL, selected=NULL, multiple=TRUE, options = list( create=TRUE,selectOnTab = TRUE, closeAfterSelect = TRUE)) ,

           br(),
           fileInput("newFile", "Upload list of ids - 1 per line.", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv", ".txt")),
           br(),
           br(),
           hr(),
           sliderInput("plotWidth", "Width of plot (in pixels)", min=200, max=4000, value=1200, step=100),
           sliderInput('plotHeight', 'Height of plot (in pixels)', min = 200, max = 4000, value = 600, step=100),
           br(),
           br(),
           h6("We explored two environmental cues: photocycle (L:light, and D:dark, 12h:12h) and thermocycle (H:hot 28C, and C:cold 12C, 12h:12h). Accession Bd21 seeds were grown on 100 mL 1x MS + 1% Agar in Magenta boxes. All seeds were germinated in dark for 3 days then transferred to LDHC for 10 days. Seedlings were then released into one of four conditions: LDHC, LDHH, LLHC or LLHH. 12 hours after placing the samples in each condition the currently emerging leaf was sampled and pooled from at least 3 separate plants. 13 additional time points were collected following this first sample (total 14 time points), one every 3.5 hours. RNA-seq was performed on all samples to quantify transcript abundace.")
    )

    ),
  
  mainPanel(
  hr(),
  fluidRow(
    
           
           h6(verbatimTextOutput("defline"))
           
  ),
  hr(),
  
  tabsetPanel(type="tabs",
               
                                tabPanel("Plot", value=1, br(), 
                                         conditionalPanel( condition = "input.bradi != ''",
                                         h4("Values shown are DESeq2 normalized read counts."),                      
                                         br(),
                      
                                         tabsetPanel(type = "pills" ,
                      tabPanel("All", align="center", br(), br(), plotOutput("graph", height = "auto")),
                       
                       tabPanel("By Condition",  br(), br(),plotOutput("faceted", height = "auto")),
                       
                       tabPanel("By ID",  br(), br(),plotOutput("face2", height = "auto"))
                      )
                                )
                      ),
              
              tabPanel("Scaled", value=2,
                       h4("Normalized read counts scaled by z-score transform."),
                       tabsetPanel(type = "pills" ,
                       
                                   tabPanel("All",plotOutput("sgraph", height = "auto" )),
           
                        
                                   tabPanel("By Condition",plotOutput("sfaceted", height = "auto" )),
           
                        
                                   tabPanel("By ID",plotOutput("sface2" , height = "auto"))
                       )
              ),
              
              
              
              
              
              tabPanel("MetaCycle Output",
                       h3("MetaCycle Output:"),
                       br(),
                       tableOutput("mytable")
              ),
              tabPanel("Variance",
                       h3("Coefficient of Variance:"),
                       br(),
                       h6("(Standard deviation / mean) reported as percentage."),
                       tableOutput("variance")
              ),
              
             
              
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
	


server <- function(input, output, session) {
  
  
	observe({
	  query <- parseQueryString(session$clientData$url_search)
	  if (!is.null(query[['id']])) {
	    updateSelectizeInput(session, "bradi", choices = query[['id']], selected= query[['id']], server = TRUE, options = list(create=TRUE))
	  }
	  if ( !is.null(input$newFile))
	  {
	    try(ids <- fread(input$newFile$datapath, header=FALSE))
	    updateSelectizeInput(session, "bradi",  choices=mas$CycID, selected= ids$V1, server = TRUE, options = list(create=TRUE))
	    
	  }
	 
	})
	

	  updateSelectizeInput(session, "bradi", choices = mas$CycID, server = TRUE ,options = list(create=TRUE))
	  

	
	searchTerm <- reactive({
	  
	  if (is.null(input$bradi)) return()
	  if (nchar(input$bradi) < 11) return()
	  
	  trimmed <- trimws(input$bradi, which="both")
	  trimmed <- gsub(' ', '', trimmed)
	  temp <- mas[paste0(trimmed)]
	  if ( is.na(temp$Condition) == TRUE ) {
	    temp <- mas[ CycID %ilike% paste0(trimmed)]
	  }
	  
	  temp <- melt.data.table(temp, measure.vars = patterns("T_") )
	  temp[, variable := gsub("T_", "", variable)]
	  temp[, variable := as.numeric(variable)]
	  temp[, scaled := scale(value), by=.(CycID, Condition)]
	  temp
	  })

	
	# output$caption <- renderText({
	#   if ( is.null(searchTerm()$Condition) ) return() 
	#   paste("<h3> <a href=\"", firstlink, searchTerm()$`#pacId`[1], "\" target=\"_blank\">", input$bradi, "</a> </h3>", sep='')
	#   
	#   })

	output$defline <- renderText({
	  if ( is.null(searchTerm()$Condition) ) return() 
		unique(paste(searchTerm()$CycID, searchTerm()$`Best-hit-arabi-name`, searchTerm()$`arabi-defline`, "\n", sep="\t"))
	  
	 })
	
	output$mytable <- renderTable({
	  if ( is.null(searchTerm()$Condition)) return() 
	  temp <- unique(searchTerm()[ , names(searchTerm()) %like% "CycID|Condition|meta", with=FALSE])
	  temp[, meta2d_pvalue := format(as.numeric(meta2d_pvalue), digits=3, scientific = TRUE)]
	  temp
	  
	})	
	
	output$variance <- renderTable({
	  if (is.null(searchTerm()$Condition)) return()
	  dcast.data.table(searchTerm()[, .(cv = (sd(value) / mean(value) * 100)), by=.(CycID, Condition)], CycID ~ Condition, value.var = "cv")
	  
	  
	})
	
	output$expressionTable <- renderTable({
	  if (is.null(searchTerm()$Condition)) return() 
	  dcast.data.table(searchTerm(), CycID + Condition ~ variable, value.var = "value") 

	  })
	
	
	
	output$graph <- renderPlot({ 
	  if ( is.null(searchTerm()$Condition)) return()
	  
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
	  
		}, res=300, height = exprToFunction(input$plotHeight), width = exprToFunction(input$plotWidth))
	
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

shinyApp(ui,server)
