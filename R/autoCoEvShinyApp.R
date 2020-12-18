# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

#install and load necessary packages
#if packages already installed and loaded, skip.
#But if packages installed and not loaded then load
#If packages not installed, install and load.
packages = c("igraph","shiny","visNetwork","DT","dplyr")


## load or install packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)




###in stall and load bioconductor package
InstalledPackages <- as.data.frame(installed.packages())
ToInstall<-as.data.frame(c("clusterProfiler","STRINGdb","GeneNetworkBuilder","org.Mm.eg.db"))
BiocManager::install(setdiff(ToInstall[,1],InstalledPackages[,1]),update=FALSE)
#BiocManager::install("clusterProfiler","STRINGdb","GeneNetworkBuilder","org.Mm.eg.db")
library(clusterProfiler)
library(STRINGdb)
library(GeneNetworkBuilder)
library(org.Mm.eg.db)



###############Global
filePath<-filePath
data<-read.table(file = filePath, header = TRUE)%>%
  dplyr::select(Name1,msa1,seq1,Name2,msa2,seq2,bonferroni)%>%
  dplyr::rename("from"="Name1","to"="Name2","UniprotId1"="msa1", "UniprotId2"="msa2")
string_db <- STRINGdb$new( version="11", species=10090,
                           score_threshold=400)

##Nodes
vert<- dplyr::select(data,from,bonferroni)
vert2<-dplyr::select(data,to,bonferroni)%>%dplyr::rename("from"="to")
nodes<-rbind(vert,vert2)%>%dplyr::rename("id"="from")

##Edges
edge<-dplyr::select(data, from, to, bonferroni,UniprotId1,UniprotId2)%>%dplyr::rename("width"="bonferroni")
edge$title=edge$width

###############ui

  shinyApp(
ui <- fluidPage(
  
  
  #Application title
  titlePanel("Predicted co-evolving proteins"),
  

  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}","h3{
      color:blue;
    }" ))
    ),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
      #side panel
      sidebarPanel(
        tags$h3("Network Settings"),
        sliderInput("bins",
          "Filter with bonferroni Adjusted P-values:",
          min = min(data$bonferroni),
          max = max(data$bonferroni),
          value = 0.05,
          #step =  0.000001,
          #width = "700px"
          ),
        #Dropdown selection input for node shapes
        selectInput(inputId="shape",
          label = "Select node shape:",
          choices=c('ellipse', 'circle', 'box','text'),selected='ellipse',width = "400px"),
        
        #Dropdown selection input for node color
        selectInput(inputId="nodeColor",
          label = "Select node color:",
          choices=c('lightblue','blue','lightgreen','green','red','grey', 'pink',
                    'magenta','yellow','orange','darkred', 'darkblue','purple',
                    'brown'),
          selected='lightgreen',
          width = "400px"),
        
        #Dropdown selection input for edge color
        selectInput(inputId="edgeColor",
          label = "Select edge color:",
          choices=c('lightblue','blue','lightgreen','green','red','grey', 'pink',
                    'magenta','yellow','orange','darkred', 'darkblue','purple',
                    'brown'),
          selected='grey',
          width = "400px"),
        
        #Dropdown selection input for text color
        selectInput(inputId="textColor",
          label = "Select node text color:",
          choices=c('black','lightblue','blue','lightgreen','green','red','grey',
                    'pink','magenta','yellow','orange','darkred', 'darkblue','purple', 'brown'),
          selected='black',
          width = "400px"),
        
        #Dropdown selection input for node's border color
        selectInput(inputId="nodeBorderColor",
          label = "Select node border color:",
          choices=c('black','lightblue','blue','lightgreen','green','red','grey',
                    'pink','magenta','yellow','orange','darkred', 'darkblue','purple',
                    'brown'),
          selected='green',
          width = "400px"),
        
        #numerical input text size
        numericInput("labelFont", "Label Font:", 20, min = 8, max = 100),
        tags$hr(),
        tags$h3("Functional Analysis Settings"),
        radioButtons(inputId="Specie",
                     label = "Select Specie:",
                     choices=c('Human','Mouse'),
                     selected='Mouse',
                     width = "100px"),
        radioButtons(inputId="analysis",
                    label = "Select analysis type:",
                    choices=c('Biological process','Cellular components', 'KEGG'),
                    selected='Biological process',
                    width = "400px"),
        
        ),
      
      
      
      
      #main panel
      mainPanel(
        
        tabsetPanel(type = "pills",
          tabPanel("AutoCoEv Network",visNetworkOutput("network", width = "100%", height = "600px"),
                   verbatimTextOutput(outputId = "stats"),
                   #download botton 
                   downloadButton(outputId = "downloadTable",label = "Download Table"),
                   
                   #render reactive node and edge dataframe  
                   tabsetPanel(id = 'dataset',
                               tabPanel("Edge Table", DT::dataTableOutput("mytable1")),
                               tabPanel("Node Table", DT::dataTableOutput("mytable2"))
                   )
                   ),
          tabPanel("Merged AutoCoEv/Stringdb Network",visNetworkOutput("networkstring", width = "100%", height = "600px"),
                   # tabsetPanel(id = 'dataset2',
                   #             tabPanel("Edges", DT::dataTableOutput("mytable3")),
                   #             tabPanel("Nodes", DT::dataTableOutput("mytable4"))
                   # )
                   ),
          # tabPanel("Merge Network",visNetworkOutput("networkmerge", width = "100%", height = "600px"),
          #          tabsetPanel(id = 'dataset2')
          # ),
          tabPanel("Functional Information",DT::dataTableOutput("mytable5")
          )
        ),
        #tags$hr(),
          
        # #download botton 
        # downloadButton(outputId = "downloadTable",label = "Download Table"),
        #   
        # #render reactive node and edge dataframe  
        # tabsetPanel(id = 'dataset',
        #   tabPanel("Edges", DT::dataTableOutput("mytable1")),
        #   tabPanel("Nodes", DT::dataTableOutput("mytable2"))
        #   ),
        tags$hr(),
          
        #article link
        uiOutput(outputId = "mainArticle")
  
      )
)),


###############server side
server <- function(input, output) {
  
  #output number of nodes and edges
  output$stats<-renderPrint({
    nodes<-dplyr::filter(nodes,nodes$bonferroni<=input$bins)
    nodes<-nodes[!duplicated(nodes[,'id']),]%>%dplyr::select(id)
    edge<-dplyr::filter(edge,edge$width<=input$bins)
    
    paste("Nodes: ",length(nodes$id)," ;" ," Interactions: ", length(edge$from), sep="" )
  })
  
  #download edge file
  output$downloadTable<-downloadHandler(
    'edges.csv',
    content = function(file){
      write.csv((dplyr::filter(data, data$bonferroni<=input$bins)), quote = FALSE,file)
    }
  )
  
  #output edge table
  output$mytable1 <- DT::renderDataTable({
    filterEdgeTable<-dplyr::filter(data, data$bonferroni<=input$bins)
    #link to uniprot
    filterEdgeTable$UniprotId1 <- paste0("<a href='https://www.uniprot.org/uniprot/",filterEdgeTable$UniprotId1,"'target='_blank'>",filterEdgeTable$UniprotId1,"</a>")
    filterEdgeTable$UniprotId2<- paste0("<a href='https://www.uniprot.org/uniprot/",filterEdgeTable$UniprotId2,"'target='_blank'>",filterEdgeTable$UniprotId2,"</a>")
    
    filterEdgeTable$from <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",filterEdgeTable$from,"'target='_blank'>",filterEdgeTable$from,"</a>")
    filterEdgeTable$to<- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",filterEdgeTable$to,"'target='_blank'>",filterEdgeTable$to,"</a>")
    
    DT::datatable(filterEdgeTable, options = list(orderClasses = TRUE), escape = FALSE
                                    )
  })
  
  #output node table
  output$mytable2 <- DT::renderDataTable({
    nodes<-dplyr::filter(nodes,nodes$bonferroni<=input$bins)
    nodes<-nodes[!duplicated(nodes[,'id']),]%>%dplyr::select(id)%>%dplyr::rename("Gene name"="id")
    #filterNodeTable<-filter(nodes, data$bonferroni<=input$bins)
    DT::datatable(nodes, options = list(orderClasses = TRUE), width='5px',rownames=FALSE)
  })
  
  
  #article hyperlink 
  output$mainArticle<-renderPrint({
    tags$a(href="https://www.biorxiv.org/content/10.1101/2020.09.29.315374v1", 
           "Petrov, P. B., Sustar, V*., Awoniyi, L. O.*, Tolvanen, M., & Mattila, P. K. (2020). AutoCoEv-a high-throughput in silico pipeline for revealing novel protein-protein interactions. BioRxiv, 2020.09.29.315374.https://doi.org/10.1101/2020.09.29.315374"
    , target="_blank")
  } )

  
  #plot string network
  output$networkstring <- renderVisNetwork({
    #nodes<-dplyr::filter(nodes,nodes$bonferroni<=0.05)


    filterEdge<-dplyr::filter(edge,edge$width<=0.05)
    edgeStng<- dplyr::select(filterEdge,from)
    edgeStng2<-dplyr::select(filterEdge,to)%>%dplyr::rename("from"="to")
    edgeStngs<-rbind(edgeStng,edgeStng2)%>%dplyr::rename("id"="from")
    edgeStngs<-data.frame(edgeStngs[!duplicated(edgeStngs[,'id']),])%>%
      dplyr::rename("id"="edgeStngs..duplicated.edgeStngs....id......")

    example1_mapped <- string_db$map( edgeStngs, "id", removeUnmappedRows = TRUE)
    hits<-example1_mapped$STRING_id
    edgesString <- string_db$get_interactions(example1_mapped$STRING_id)
    IDsMap <- example1_mapped$id
    names(IDsMap) <- example1_mapped$STRING_id
    cifNetwork2 <- convertID(edgesString, IDsMap)
    #enrichment <- string_db$get_enrichment( hits, category = "Component", iea = TRUE)


    ##Nodes
    vertStng<- dplyr::select(cifNetwork2,from)
    vertStng2<-dplyr::select(cifNetwork2,to)%>%dplyr::rename("from"="to")
    nodesStng<-rbind(vertStng,vertStng2)%>%dplyr::rename("id"="from")
    nodesStng$label = nodesStng$id
    nodesStng$title = nodesStng$id
    nodesStng$value = nodesStng$id
    nodesStng<-nodesStng[!duplicated(nodesStng[,'id']),]%>%dplyr::select(id)
    
    ###from autocoev edge
    edge<-dplyr::filter(edge,edge$width<=0.05)
    ###
    
    
    
    net1<-dplyr::select(cifNetwork2,from, to)
    net2<-dplyr::select(edge,from, to)
    net2<-net2[!duplicated(net2[,c('from','to')]),]
    
    
    net3<-anti_join(net1,net2)
    net4<-dplyr::select(net3,to,from)%>%dplyr::rename(from=to, to=from)
    g5<-anti_join(net4,net2)%>%dplyr::select(to,from)%>%
      dplyr::rename(from=to, to=from)
    ###################################################incomplete
    # net6<-anti_join(net2,net1)
    # net7<-dplyr::select(net6,to,from)%>%dplyr::rename(from=to, to=from)
    # 
    # g6<-anti_join(net7,net1)%>%dplyr::select(to,from)%>%
    #   dplyr::rename(from=to, to=from)
    # 
    # 
    # g5$color="red"
    # net2$color="blue"
    # g6$color="green"
    #############################################
    
    g5$color="red"
    net2$color="blue"
    
    edgetest<-rbind(g5,net2)


    visNetwork(nodesStng, edgetest, width="100%")%>%
      visIgraphLayout()%>%
      visNodes(
        shape = input$shape,
        color = list(
          background = input$nodeColor,
          border = input$nodeBorderColor,
          highlight = "#FF8000"
          #size=data$Best.P.value
        ),

        shadow = list(enabled = TRUE, size = 10),
        font=list(color=input$textColor,
                  size=input$labelFont)
      ) %>%
      visEdges(
        shadow = FALSE,
        #color = list(color = input$edgeColor, highlight = "#C62F4B")
        #width=edge$width
      )%>%
      visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
      visLayout(randomSeed = 11)%>%
      visInteraction(navigationButtons = TRUE)%>%
      visExport(type = "png", name = "network", float = "left")

  })
  
  # output$mytable3 <- DT::renderDataTable({
  #   DT::datatable(enrichment, options = list(orderClasses = TRUE), width='5px',rownames=FALSE)
  # })
  ##funtional analysis
  output$mytable5 <- DT::renderDataTable({
    filterEdge<-dplyr::filter(edge,edge$width<=input$bins)
    edgeStng<- dplyr::select(filterEdge,UniprotId1)
    edgeStng2<-dplyr::select(filterEdge,UniprotId2)%>%dplyr::rename("UniprotId1"="UniprotId2")
    edgeStngs<-rbind(edgeStng,edgeStng2)%>%dplyr::rename("id"="UniprotId1")
    edgeStngs<-data.frame(edgeStngs[!duplicated(edgeStngs[,'id']),])%>%
      dplyr::rename("id"="edgeStngs..duplicated.edgeStngs....id......")
    protein_list<-edgeStngs$id%>%unlist()%>%as.character()
    ego_cc<-enrichGO(gene = protein_list,
                     OrgDb = org.Mm.eg.db,
                     keyType = "UNIPROT",
                     ont = "BP",
                     minGSSize = 2,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE
                     
                       )
    ego_df<-data.frame(ego_cc)
    DT::datatable( ego_df,options = list(pageLength = 5, scrollX = T))
  })
  
  
  #plot Autocoev network  
  output$network <- renderVisNetwork({
    #define node
    nodes<-dplyr::filter(nodes,nodes$bonferroni<=input$bins)
    nodes<-nodes[!duplicated(nodes[,'id']),]%>%dplyr::select(id)
    nodes$label = nodes$id
    nodes$title = nodes$id
    nodes$value = nodes$id
    #define edge
    edge<-dplyr::filter(edge,edge$width<=input$bins)
   
    
    
    visNetwork(nodes, edge, width="100%")%>%
      visIgraphLayout()%>%
      visNodes(
        shape = input$shape,
        color = list(
          background = input$nodeColor,
          border = input$nodeBorderColor,
          highlight = "#FF8000"
          #size=data$Best.P.value
          ),
        
        shadow = list(enabled = TRUE, size = 10),
        font=list(color=input$textColor,
                  size=input$labelFont)
      ) %>%
      visEdges(
        shadow = FALSE,
        color = list(color = input$edgeColor, highlight = "#C62F4B")
        #width=edge$width
      )%>%
      visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
      visLayout(randomSeed = 11)%>%
    visInteraction(navigationButtons = TRUE)%>%
      visExport(type = "png", name = "network", float = "left")
    
  })
  
  
  # 
  # #plot Autocoev network  
  # output$networkmerge <- renderVisNetwork({
  #   #define node
  #   nodes<-dplyr::filter(nodes,nodes$bonferroni<=input$bins)
  #   nodes<-nodes[!duplicated(nodes[,'id']),]%>%dplyr::select(id)
  #   nodes$label = nodes$id
  #   nodes$title = nodes$id
  #   nodes$value = nodes$id
  #   #define edge
  #   edge<-dplyr::filter(edge,edge$width<=input$bins)
  #   
  #   
  #   
  #   visNetwork(nodes, edge, width="100%")%>%
  #     visIgraphLayout()%>%
  #     visNodes(
  #       shape = input$shape,
  #       color = list(
  #         background = input$nodeColor,
  #         border = input$nodeBorderColor,
  #         highlight = "#FF8000"
  #         #size=data$Best.P.value
  #       ),
  #       
  #       shadow = list(enabled = TRUE, size = 10),
  #       font=list(color=input$textColor,
  #                 size=input$labelFont)
  #     ) %>%
  #     visEdges(
  #       shadow = FALSE,
  #       color = list(color = input$edgeColor, highlight = "#C62F4B")
  #       #width=edge$width
  #     )%>%
  #     visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  #     visLayout(randomSeed = 11)%>%
  #     visInteraction(navigationButtons = TRUE)%>%
  #     visExport(type = "png", name = "network", float = "left")
  #   
  # })
  
    })

# Run the application 
#shinyApp(ui = ui, server = server)
##########################################################
# Pass input and app location as arguments
# args <- commandArgs(trailingOnly = TRUE)
# filePath <- args[1]
# appLocation<-args[2]
# 
# shiny::runApp(appLocation)


