library(shiny)
library(shinydashboard)
#library(rcytoscapejs)
library(WGCNA)
library(igraph)
#library(InteractiveIGraph)
library(networkD3)
library(sendmailR)
#library(mailR)
library(htmlwidgets)
#library(limma)
options(shiny.maxRequestSize=30*1024^2)
uiactual <- dashboardPage(
  dashboardHeader(title="WGCNA shiny"),
  dashboardSidebar(menuItem(
    checkboxInput('header', 'Header', TRUE)),
    menuItem(radioButtons('sep', 'Separator',
                          c(Space=' ',
                            Comma=',',
                            Semicolon=';',
                            Tab='\t'), selected='\t')),
    menuItem(p('If you want a sample .csv or .tsv file to upload,',
               'you can first download the sample',
               downloadLink("table", "table.txt")
    )),
    menuItem(checkboxInput(inputId = "gene_data", label = strong("Gene-expression data"), value = T)
    ),
    menuItem(selectInput(inputId = "sign", label = strong("Signed network"), choices = c("unsigned", "signed"), selected = "unsigned")
    ),
    menuItem(p("Please upload a dataframe file that contains values of type 'numeric'. You can specify the file's format ('.txt' or '.csv') delimiter via the options displyaed in the sidebar menu.")),
    #radioButtons("reading", "File type:",
    #             c("Text file" = "table_file",
    #               "CSV file" = "csv"
    #             ), selected = "Text file"),
    menuItem(fileInput('file1', 'Choose data-file to upload',
                       accept = c(
                         'text/csv',
                         'text/comma-separated-values',
                         'text/tab-separated-values',
                         'text/plain',
                         '.csv',
                         '.tsv'
                       )
    )),
    menuItem(actionButton("choices", "run analyses"))
  ),
  dashboardBody(
    h2("Network Inference and Subsequent Display"),
    p("Before initiating the clustering process, please choose a minimum-allowed cluster size. You can come backand change this later and the graph results will change interactively.", style = "font-family: 'times'; font-size: 12pt"),
    sliderInput(inputId="cluster_size", label="Minimum cluster size", value=10, min=10, max=30),
    numericInput(inputId = "cutoff_filter", label="Choose cutoff for tree clustering distance", value=5000),
    plotOutput("plot_filter"),
    actionButton(inputId="go", label="Run post-filter analyses"),
    plotOutput("plot_power"),
    #verbatimTextOutput("recommended_power"),
    p("Considering the above output, select a soft threshold power for the network to reach scale-free topology.", style = "font-family: 'times'; font-size: 12pt"),
    sliderInput(inputId="num", label="Choose softthreshold for scale-free topology", value=9, min=1, max=20),
    plotOutput("plot_modules"),
    sliderInput(inputId="dendocut", label="CutHeight for 'cutreeDynamic'", value=0.99, min=0, max=0.99),
    #rcytoscapejsOutput("g3plot"),
    selectInput(inputId = "basis_choice", label = "Choose whether the network should be plotted according to the adjacency values, the Topological Overlap Measures (TOM) or the dissimilarity between the TOM values", choices = c("adjacency values"="ad", "TOM values"="tom", "TOM dissimilarity values"="tomdiss"), selected="adjacency values"),
    plotOutput("graph"),
    tableOutput("eigen"),
    p("average degree per node below:",  style = "font-family: 'times'; font-size: 12pt"),
    verbatimTextOutput("average_deg"),
    #actionButton("labeleromama", "show labels"),
    #radioButtons("l", "display labels:",
    #             c("No" = "none",
    #               "Yes" = "lab"), selected = "No"),
    actionButton("send", "Save interactive network"),
    #textInput("to", "Send to:"),
    radioButtons("type", "Distribution type:",
                 c("Auto" = "norm",
                   "Fruchtmann-Reingold" = "r",
                   "DRL" = "drl",
                   "Kamada Kawai" = "layout"), selected = "r"),
    #selectInput(inputId = "adj", "choose graph type", choices=c("layout.fruchterman.reingold", "layout.drl")),
    
    sliderInput(inputId="adjust", label="Choose adjacency threshold", value=0.01, min=0.001, max=1),
    #sliderInput(inputId="adjust2", label="Choose adjacency threshold", value=0.3, min=0, max=1),
    sliderInput(inputId="penalty", label="Choose node-size connectedness enrichment", value=1, min=1, max=50),
    sliderInput(inputId="enrich", label="Choose edge-enrichment", value=20, min=1, max=50),
    plotOutput("MErelations"),
    plotOutput("hist"),
    #plotOutput("scale_free_plot"),
    h2("Correlating traits to Module-eigengenes"),
    p("If available you can upload a second file that gives information on each of the sample rows and this can then be correlated to the module eigengenes to give an estimate which modules are most implicated in whatever particular phenotype etc. Please ensure that the file-format is consistent with the ", style = "font-family: 'times'; font-size: 12pt"),
    fileInput('file2', 'Choose info-file to upload',
              accept = c(
                'text/csv',
                'text/comma-separated-values',
                'text/tab-separated-values',
                'text/plain',
                '.csv',
                '.tsv'
              )
    ),
    actionButton("choices2", "incorporate external information"),
    p("Please select a columns and the graphs will change and analyse a particular column", style = "font-family: 'times'; font-size: 12pt"),
    selectInput("columnsfirst","Select Columns", choices = NULL),
    plotOutput("GS"),
    selectInput(inputId = "map_choice", label = "choose whether to plot as a heatmap or a series of regression analyses", choices = c("heatmap_and_dendogram", "matrix_of_regression_analyses"), selected = "heat"),
    plotOutput("module_heatmap"),
    p("Multipe columns can be selected here for analyis and this gives rise to a heatmap documenting the selected traits correlation with the module eigengenes.", style = "font-family: 'times'; font-size: 12pt"),
    selectInput("columns","Select Columns", choices = NULL, multiple = TRUE),
    actionButton(inputId = "go2", "Render Multiple-trait Heatmap"),
    plotOutput("heatmap_multiple")
    #tags$style(type="text/css",
    #           ".shiny-output-error { visibility: hidden; }",
    #           ".shiny-output-error:before { visibility: hidden; }"
    #)
  )
)

serveractual <- function(input, output, session) {
  #read <- reactive({
  #  read_style <- switch(input$reading,
  #                       table_file = read.table,
  #                       csv = read.csv)
  #})
  dataInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    isolate(t1<-read.table(inFile$datapath, header = input$header,
                           sep = input$sep))
    isolate(t<-apply(t1,2, as.numeric)) #constituents read as character strings need to be converted
    rownames(t)<-rownames(t1)
    if(input$gene_data){
      gsg<-goodSamplesGenes(t, verbose=3)
      if (!gsg$allOK)
      {
        if (sum(!gsg$goodGenes)>0) 
          printFlush(paste("Removing genes:", paste(names(t)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0) 
          printFlush(paste("Removing samples:", paste(rownames(t)[!gsg$goodSamples], collapse = ", ")));
        t = t[gsg$goodSamples, gsg$goodGenes]
      }
    }
    t
  })
  sampletree<-eventReactive(input$choices, {
    hclust(dist(dataInput()), method = "average")
  })
  output$plot_filter<-renderPlot({
    input$choices
    if (is.null(dataInput()))
      return(NULL)
    withProgress(message = 'Please Wait',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   plot(sampletree(), main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
                   abline(h = input$cutoff_filter, col = "red")
                   #legend("bottomleft", scc$csize>1, pt.bg=unique(node_colors), pch=21)
                 })
    
  })
  dataInput1<-eventReactive(input$go, {
    input$file1
    if (is.null(dataInput()))
      return(NULL)
    dataInput<-dataInput()
    withProgress(message = 'Removing Outlying Samples',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   clust = cutreeStatic(sampletree(), cutHeight = input$cutoff_filter, minSize = 10)
                   keepSamples = (clust==1)
                   dataInput[keepSamples, ]
                   #legend("bottomleft", scc$csize>1, pt.bg=unique(node_colors), pch=21)
                 })
    
  })
  adjacentmat<-reactive({
    adj_mat<-adjacency(as.matrix(dataInput1()),power=input$num, type=input$sign)
    if(input$sign=="unsigned"){
      adj_mat[adj_mat < input$adjust] <- 0
      #adj_mat[adj_mat> input$adjust2] <-0
    }
    adj_mat[adj_mat > 0.99] <- 0
    adj_mat
  })
  lay <- reactive({
    layout_style <- switch(input$type,
                           norm = layout.auto,
                           r = layout.fruchterman.reingold,
                           drl = layout.drl,
                           layout=layout.kamada.kawai)
  })
  sft<-reactive({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    withProgress(message = 'Calculating Soft thresholds',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   pickSoftThreshold(dataInput1(), powerVector = c(c(1:10), seq(from = 12, to=20, by=2)), verbose = 5)
                 })
  })
  #output$recommended_power<-renderText({
  #  sft()$powerEstimate
  #})
  output$plot_power<-renderPlot({
    sft<-sft()
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=1.65,col="red")
  })
  TOM<-reactive({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    #adj_mat<-adjacency(as.matrix(dataInput1()),power=input$num, type=input$sign)
    #adj_mat[adj_mat < input$adjust] <- 0
    #adj_mat[adj_mat> input$adjust2] <-0
    #adj_mat[adj_mat > 0.99] <- 0
    withProgress(message = 'Calculating TOMsimilarity and dissimilarity',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   TOM = TOMsimilarity(adjacentmat(), TOMType = input$sign)
                   #legend("bottomleft", scc$csize>1, pt.bg=unique(node_colors), pch=21)
                 })
    
  })
  diss<-reactive({
    1-TOM()
  })
  output$table <- downloadHandler(
    filename = function() {
      "table.txt"
    },
    content = function(file) {
      file.copy("table.txt", file)
    }
  )
  #dataInput2<-reactive({
  #  ADJ1=abs(adjacentmat())^6
  #  dissADJ1<-1-ADJ1
  #  hclust(as.dist(dissADJ1), method = "average")
  #})
  
  #dataInput2<-reactive({
  #  withProgress(message = 'Hierarchical clustering',
  #               detail = 'This may take a while...', value = 0, {
  #                 for (i in 1:15) {
  #                   incProgress(1/15)
  #                   Sys.sleep(0.25)
  #                 }
  #                 hclust(as.dist(diss()), method = "average")
  #               })
  
  #  })
  clusters<-reactive({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    withProgress(message = 'Hierachical Clustering',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   hclust(as.dist(diss()), method = "average")
                 })
  })
  dynamic_tree<-reactive({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    withProgress(message = 'Identifying Modules',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   cutreeDynamic(dendro = clusters(), distM = diss(), cutHeight = input$dendocut,
                                 deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = input$cluster_size)
                 })
    
  })
  
  
  output$plot_modules<-renderPlot({
    plotDendroAndColors(clusters(), colors = labels2colors(dynamic_tree()), dendroLabels = FALSE, hang = 0.03, main = "Gene hierarchical clustering dendrogram and simulated module colors")
  })
  merges<-reactive({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    #minModuleSize = input$cluster_size
    # Module identification using dynamic tree cut: 
    dynamicColors = labels2colors(dynamic_tree())
    MEList = moduleEigengenes(dataInput1(), colors = dynamicColors, softPower = input$num) 
    MEs = MEList$eigengenes 
    # Calculate dissimilarity of module eigengenes 
    MEDiss = 1-cor(MEs);  ### Is this correct?? Should use absolute correlation ????  YES ###
    # Cluster module eigengenes 
    withProgress(message = 'Module-Eigengene summarisation',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   METree = hclust(as.dist(MEDiss), method = "average")
                 })
    
  })
  output$MErelations<-renderPlot({
    plot(merges(), main = "Clustering of module eigengenes", xlab = "", sub = "")
  })
  node_colors<-reactive({
    labels2colors(dynamic_tree())
  })
  node_colors_no_grey<-reactive({node_colors()[node_colors()!="grey"]})
  g_input<-reactive({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    adj_mat<-adjacency(as.matrix(dataInput1()),power=input$num, type=input$sign)
    if(input$sign=="unsigned"){
      adj_mat[adj_mat < input$adjust] <- 0
      #adj_mat[adj_mat> input$adjust2] <-0
    }
    adj_mat[adj_mat > 0.99] <- 0
    TOM<-TOM()
    diss<-diss()
    #g<-g_input()
    if(input$basis_choice=="ad"){
      g <- simplify(graph.adjacency(adj_mat, mode='undirected', weighted=TRUE,  add.colnames=NA))
    }
    else if(input$basis_choice=="tom"){
      g <- simplify(graph.adjacency(TOM, mode='undirected', weighted=TRUE,  add.colnames=NA))
    }
    else if(input$basis_choice=="tomdiss"){
      g <- simplify(graph.adjacency(diss, mode='undirected', weighted=TRUE,  add.colnames=NA))
    }
    #node_colors <- labels2colors(dynamic_tree())
    #g <- simplify(graph.adjacency(adj_mat, mode='undirected', weighted=TRUE,  add.colnames=NA))
    v<-V(g)
    node_colors<-node_colors()
    v$color <- node_colors
    g<-delete_vertices(g, v[node_colors=="grey"])
    withProgress(message = 'Deleting Nodes of grey module from displayed network',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   g=delete_vertices(g,which(degree(g)<input$adjust))
                   #g=delete_vertices(g,which(degree(g)>input$adjust2))
                   g
                 })
  })
  output$eigen<-renderTable({
    as.data.frame(cbind(node_name=rownames(adjacentmat())[order(betweenness(g_input()), decreasing = T)[1:10]], betweenness=betweenness(g_input())[order(betweenness(g_input()), decreasing = T)[1:10]], colours=node_colors_no_grey()[order(betweenness(g_input()), decreasing = T)[1:10]]))
  })
  output$average_deg<-renderText({mean(degree(g_input()))})
  size<-reactive({
    g<-g_input()
    deg <- degree(g, mode="all")
    V(g)$size <- input$penalty*deg
    as.vector(V(g)$size)
  })
  labs<- reactive({
    if(input$l=='none'){ans=NA}
    else{ans=rownames(adjacentmat())[which(node_colors()!="grey")]}                         
    ans
    })
  members<-reactive({
    g<-g_input()
    members<-cluster_walktrap(g)
    members<-membership(members)
  })
  output$graph<-renderPlot({
    input$go
    if (is.null(dataInput1()))
      return(NULL)
    g<-g_input()
    coords_fr = lay()(g)
    igraph.options(vertex.size=9, edge.width=(log(E(g)$weight)*input$enrich))
    E(g)$arrow.mode <- 0
    deg <- degree(g, mode="all")
    V(g)$size <- input$penalty*log(deg)
    #g = InteractiveIGraph.Constructor(g)
    node_colors<-node_colors()
    node_colors<-subset(node_colors, node_colors!="grey")
    withProgress(message = 'Compiling Network for display',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   plot.igraph(g, layout=coords_fr, vertex.color=node_colors, mark.groups = members(), vertex.label = NA)
                   #legend("bottomleft", scc$csize>1, pt.bg=unique(node_colors), pch=21)
                 })
  })
  #output$scale_free_plot<-renderPlot({
  #  ADJ1=abs(cor(dataInput1(),use="p"))^6
  # When you have relatively few genes (<5000) use the following code
  #  k=as.vector(apply(ADJ1,2,sum, na.rm=T))
  # When you have a lot of genes use the following code
  #  k=softConnectivity(datE=dataInput1(),power=6) 
  # Plot a histogram of k and a scale free topology plot
  #  scaleFreePlot(k, main="Check scale free topology\n")
  #})
  data<-reactive({
    adj_mat<-adjacentmat()
    m<-adj_mat
    if(input$sign=="unsigned"){
      idx = m > input$adjust
      data.frame(
        source = row(m)[idx],
        target = col(m)[idx],
        corr = m[idx])
    }
    else{
      #idx = m > input$adjust | m < input$adjust
      data.frame(
        source = row(m),
        target = col(m),
        corr = as.numeric(m))
    }
  })
  output$hist<-renderPlot({hist(as.numeric(data()$corr), main="Histogram documenting the number of edges of a strength above the chosen adjacency threshold.", xlab = "Edge-strength")})
  data2<-reactive({
    adj_mat<-adjacentmat()
    TOM<-TOM()
    diss<-diss()
    g<-g_input()
    if(input$basis_choice=="ad"){
      g <- simplify(graph.adjacency(adj_mat, mode='undirected', weighted=TRUE,  add.colnames=NA))
    }
    else if(input$basis_choice=="tom"){
      g <- simplify(graph.adjacency(TOM, mode='undirected', weighted=TRUE,  add.colnames=NA))
    }
    else if(input$basis_choice=="tomdiss"){
      g <- simplify(graph.adjacency(diss, mode='undirected', weighted=TRUE,  add.colnames=NA))
    }
    v<-V(g)
    node_colors<-labels2colors(dynamic_tree())
    v$color <- node_colors
    g<-delete_vertices(g, v[node_colors=="grey"])
    members<-cluster_walktrap(g)
    members<-membership(members)
    igraph_to_networkD3(g, group = members)
  })
  
  observe({
    if(is.null(input$send) || input$send==0) return(NULL)
    nl<-cbind(data2()$nodes, val=node_colors()[which(node_colors()!="grey")], value=rownames(adjacentmat())[which(node_colors()!="grey")], size=size())
    i<-forceNetwork(Links = data2()$links, Nodes = nl, Source = 'source', Target = 'target', NodeID = 'value', Nodesize = 'size', Group = "val", zoom = T, legend = T, opacityNoHover = TRUE)
    saveWidget(i, file="eractive_net.html")
    #write.csv(group=as.data.frame(unique(node_colors()), colour=unique(dynamic_tree())), "colour_index.csv")
    exportNetworkToCytoscape(TOM(), nodeNames = rownames(adjacentmat()), edgeFile = "edges.txt", nodeFile = "nodes.txt")
    #attachmentObject <- mime_part(x="./sctterplot2.html",name="scatterplot2.html")
    #bodyWithAttachment <- list(body,attachmentObject)
  })
  info <- eventReactive(input$choices2, {
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    isolate(t1<-read.table(inFile$datapath, header = input$header, sep = input$sep))
    vars <- names(t1)
    updateSelectInput(session, "columnsfirst","Select Columns", choices = vars)
    updateSelectInput(session, "columns","Select Columns", choices = vars)
    t1
  })
  t2<-reactive({
    t<-dataInput1()
    withProgress(message = 'Please Wait',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   t2<-t[match(intersect(rownames(t), rownames(info())), rownames(t)), ]
                   t2
                 })
  })
  y<-reactive({
    funcs2<-info()
    req(input$columnsfirst)
    isolate(y<-funcs2[match(intersect(rownames(t2()), rownames(funcs2)), rownames(t2())), ])
    as.vector(as.matrix(y[,input$columnsfirst]))
  })
  
  adjacentmat_rearranged<-reactive({
    adj_mat<-adjacency(as.matrix(t2(), power=input$num))
    adj_mat[adj_mat < 0.01] <- 0
    adj_mat[adj_mat > 0.99] <- 0
    adj_mat
  })
  diss_rearranged<-reactive({
    TOM = TOMsimilarity(adjacentmat_rearranged())
    diss<-1-TOM
  })
  output$GS<-renderPlot({
    clusters<-hclust(as.dist(diss_rearranged()), method="average")
    dynamic_tree<-cutreeDynamic(dendro = clusters, distM = diss_rearranged(), cutHeight = input$dendocut, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = input$cluster_size)
    dynamicColors<-labels2colors(dynamic_tree)
    correlate_external_info= cor(y(), t2(), use="p")
    correlate_external_info_abs<-abs(correlate_external_info)
    plotModuleSignificance(correlate_external_info_abs,dynamicColors)
  })
  
  MEs_rearranged<-reactive({
    clusters<-hclust(as.dist(diss_rearranged()), method="average")
    dynamic_tree<-cutreeDynamic(dendro = clusters, distM = diss_rearranged(), cutHeight = input$dendocut, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = input$cluster_size)
    dynamicColors<-labels2colors(dynamic_tree)
    MEs<-moduleEigengenes(t2(), dynamicColors, softPower = input$num)$eigengenes
    MEs
  })
  
  output$module_heatmap<-renderPlot({ 
    
    clusters<-hclust(as.dist(diss_rearranged()), method="average")
    dynamic_tree<-cutreeDynamic(dendro = clusters, cutHeight = input$dendocut, distM = diss_rearranged(), deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = input$cluster_size)
    dynamicColors<-labels2colors(dynamic_tree)
    correlate_external_info= cor(y(), t2(), use="p")
    correlate_external_info_abs<-abs(correlate_external_info)
    correlate<-orderMEs(cbind(MEs_rearranged(), y()))
    if(input$map_choice=="heatmap_and_dendogram"){plotEigengeneNetworks(correlate, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)}
    if(input$map_choice=="matrix_of_regression_analyses"){
      datME=moduleEigengenes(t2(),dynamicColors)$eigengenes
      signif(cor(datME, use="p"), 2)
      plotMEpairs(datME,y=y())
    }
  })
  
  y2<-eventReactive(input$go2, {
    funcs2<-info()
    req(input$columns)
    isolate(y<-funcs2[match(intersect(rownames(t2()), rownames(funcs2)), rownames(t2())), ])
    y[,input$columns]
  })
  
  output$heatmap_multiple<-renderPlot({
    t2<-t2()
    y<-y2()
    adj_mat<-adjacency(as.matrix(t2), power=input$num)
    adj_mat[adj_mat < 0.001] <- 0
    adj_mat[adj_mat > 0.99] <- 0
    TOM = TOMsimilarity(adj_mat);
    diss<-1-TOM
    clusters<-hclust(as.dist(diss), method="average")
    dynamic_tree<-cutreeDynamic(dendro = clusters, distM = diss, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = input$cluster_size)
    moduleColors = labels2colors(dynamic_tree)
    MEs0 = moduleEigengenes(t2, moduleColors, softPower = input$num)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, y, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(t2))
    #sizeGrWindow(10,6)
    textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    #par(mar = c(6, 8.5, 3, 3));
    withProgress(message = 'Rendering Multi-trait heatmap',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   labeledHeatmap(Matrix = moduleTraitCor,
                                  xLabels = names(y),
                                  yLabels = names(MEs),
                                  ySymbols = names(MEs),
                                  colorLabels = FALSE,
                                  colors = blueWhiteRed(50),
                                  textMatrix = textMatrix,
                                  setStdMargins = T,
                                  cex.text = 0.5,
                                  zlim = c(-1,1),
                                  main = paste("Module-trait relationships"))
                 })
    
  })
}
shinyApp(uiactual, serveractual)