options(shiny.maxRequestSize = 30 * 1024 ^ 2)

library(shiny)
library(shinydashboard)
library(readr)
library(DT)
library(shinyjs)

source("utils.R")

############### Global variables ##############
sess_id = ""
beds_dir <- "beds"
fasta_dir <- "fasta"
meme_out_fold <- "meme_out"



skin <- Sys.getenv("DASHBOARD_SKIN")
skin <- tolower(skin)
if (skin == "")
  skin <- "blue"

sidebar <- dashboardSidebar(sidebarMenu(
  menuItem(
    "Load Files",
    tabName = "loadFiles",
    icon = icon("dashboard")
  ),
  menuItem("Motif Analysis",
           tabName = "mot_sum",
           icon = icon("th")),
  menuItem(
    "Tfbfs motif Analysis",
    tabName = "tfbs_mot",
    icon = icon("bar-chart-o")
  )
))

body <- dashboardBody(tabItems(
  ########## load file tab ###############
  tabItem(
    "loadFiles",
    fluidRow(box(
      selectInput(label = "Select CpG Annotation",inputId = "cpgType",
                  choices = c("","Infinium Human Methylation 450K","Infinium MethylationEPIC","Custom annotation"),selected = ""), 
      shinyjs::disabled(fileInput(inputId = "CustomAnn",
        label = "Upload Custom annotation file"
      ))
    ),
    box(textOutput(outputId = "Ann_sum_info"))),
    fluidRow(box(
      fileInput(
        inputId = "backFile",
        label = "Upload CpG background file"
      )
    ),
    box(textOutput(outputId = "back_sum_info"))),
    fluidRow(box(
      fileInput(
        inputId = "targetFiles",
        label = "Upload targetFiles files",
        multiple = T
      )
    ),
    box(htmlOutput(outputId = "targ_sum_info"))),
    fluidRow(box(
      sliderInput(
        inputId = "flankSize",
        label = "flanking region size",
        min = 1,
        max = 100,
        step = 1,
        value = 20
      )
    ))
  ),
  ########## motif summary tab ###############
  tabItem(
    "mot_sum",
    fluidRow(box(
      width = 12,
      column(width = 4,
      actionButton(inputId = "startMEME", label = "Start motif analysis")),
      column(width = 4,
        numericInput(
          inputId = "unbalance_tresh",
          value = 0.7,
          min = 0.01,
          max = 0.99,
         step = 0.1,
          label = "Choose methylation unnbalannce threshold"
        )
      ),
    column(width = 4,
      downloadButton("exportSum", label = "Export summary")
    ))),
    fluidRow(box(
      width = 12,
    DT::dataTableOutput(outputId = "motifs_summary")
    )),
    fluidRow(box(
      width = 12,
      title = "Distance Heatmap",
      column(width = 12, plotOutput(outputId = "distHeat"))
    )),
    fluidRow(box(
      width = 12,
      title = "Distance dendrogram",
      column(width = 12, plotOutput(outputId = "distDend"))
    ))
  ),
  
  ########## tfbs summary tab ###############
  
  tabItem("tfbs_mot",
          fluidRow(box(
            actionButton(inputId = "startTomTom", label = "Start tfbs analysis"),
            downloadButton("exportSumTfbs", label = "Export summary tfbs")
            
          )),
          fluidRow(
            box(
              width = 12,
              DT::dataTableOutput(outputId = "motifs_summary_tfbs")
            )
          ),
          fluidRow(
            box(
              width = 12,
              title = "Tfbs Presence",
              column(width = 12, plotOutput(outputId = "tfbsHeat"))
            )
          ))
))

header <- dashboardHeader(title = "CpG Motifs")

ui <- dashboardPage(header, sidebar, body, skin = skin, shinyjs::useShinyjs())

server <- function(input, output, session) {
  ############## File Import #############
  onStop(function() {
    cat("Session stopped cleanning evironment \n")
    unlink(c("beds","fasta","meme_out"), recursive = TRUE)
  })
  
  
  #get hg19 locations of 450k probes
  hg19_loc <- reactive({
    if(input$cpgType == "Infinium Human Methylation 450K"){
      hg19_loc  <-  as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)
    }
    if(input$cpgType == "Infinium MethylationEPIC"){
      IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
      hg19_loc  <-  as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations)
    }
    if(input$cpgType == "Custom annotation"){
        hg19_loc  <-  as.data.frame(customAnn())
        rownames(hg19_loc) <- hg19_loc$cpgId
        print(nrow(hg19_loc))
    }
    hg19_loc
  })
  
  
  observeEvent(input$cpgType, {
    if(input$cpgType == "Custom annotation"){
      shinyjs::enable("CustomAnn")
    } else{
      shinyjs::disable("CustomAnn")
    }
  })  
  
  customAnn <- reactive({
    inFile <- input$CustomAnn
    if (is.null(inFile))
      return(NULL)
    read_delim(inFile$datapath, delim = "\t")
  })
  
  backFile <- reactive({
    inFile <- input$backFile
    
    if (is.null(inFile))
      return(NULL)
    
    read_delim(inFile$datapath, delim = "\t")
  })
  
  
  targFile <- reactive({
    targFiles <- list()
    inFiles <- input$targetFiles
    if (is.null(inFiles))
      return(NULL)
    for (i in 1:nrow(inFiles)) {
      targFiles[[inFiles$name[i]]] <-
        read_delim(inFiles$datapath[i], delim = "\t")
    }
    targFiles
  })
  
  
  
  ########################### input raectives ##################
  
  iDMCs <- eventReactive(targFile(), {
    iDMCs = list()
    for (targ_name in names(targFile())) {
      iDMC = targFile()[[targ_name]]
      iDMC <- iDMC[, c("ProbeID", "status")]
      
      
      #get coordinates for the Cpgs of interest
      matches = match(iDMC$ProbeID, rownames(hg19_loc()))
      iDMC = cbind(iDMC, hg19_loc()[matches, ])
      
      # store position and methylation status only
      iDMCs[[targ_name]] = iDMC
    }
    iDMCs
  })
  
  background <- eventReactive(backFile(), {
    back <-  backFile()
    back <- back[, "ProbeID"]
    #get coordinates for the Cpgs of interest
    matches = match(back$ProbeID, rownames(hg19_loc()))
    back = cbind(back, hg19_loc()[matches, ])
    
    back
  })
  
  observe({
    req(isTruthy(iDMCs()), isTruthy(background()))
    #clean directories
    #unlink(c(beds,fasta,meme_out), recursive = T)
    dir.create(beds_dir, showWarnings = F)
    dir.create(fasta_dir, showWarnings = F)
    dir.create(meme_out_fold, showWarnings = F)
    
    #write down beds for IDMCs
    for (name_idmc in names(iDMCs())) {
      #we need the genomic positions and ids only
      iDMC_pos = iDMCs()[[name_idmc]][, c("ProbeID", "chr", "pos")]
      colnames(iDMC_pos) = c("id", "chr", "pos")
      #save the position file of CpGs in the beds folder
      print(paste0(beds_dir, "/", name_idmc, ".bed"))
      write.table(
        iDMC_pos,
        paste0(beds_dir, "/", name_idmc, ".bed"),
        sep = "\t",
        row.names = F,
        quote = F
      )
    }
    
    #write down beds for background
    #save background locations in a tab delimited file
    back_pos = background()[, c("ProbeID", "chr", "pos")]
    colnames(back_pos) = c("id", "chr", "pos")
    #save the position file of CpGs in the beds folder
    write.table(
      back_pos,
      paste0(beds_dir, "/background.bed"),
      sep = "\t",
      row.names = F,
      quote = F
    )
  })
  
  meme_out <- eventReactive(input$startMEME, {
    req(isTruthy(iDMCs()), isTruthy(background()))
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Motif discovery", value = 0)
    
    print("MEME")
    i <- 1
    for (name_idmc in names(iDMCs())) {
      
      progress$inc(1/length(iDMCs()), detail = paste("Analyzing ", name_idmc))
      
      cmd_meme <-
        paste(
          "Rscript ./motifs_lib.R",
          name_idmc,
          "background",
          beds_dir,
          fasta_dir,
          meme_out_fold,
          input$flankSize
        )
      #print(cmd_meme)
      system(cmd_meme)
      Sys.sleep(1)
      i <- i + 1
    }
    get_mot_sum(iDMCs(), meme_out_fold)
    
  })
  
  
  meme_out_dir <- reactive({
    assign_main_dir(types_mot = meme_out(), input$unbalance_tresh)
  })
  
  mot_similiarity <- reactive({
    motif_similiarity(meme_out_dir())
  })
  
  
  motif_tfbs <- eventReactive(
    input$startTomTom,
    {
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Motif tfbs search", value = 0)
      
      annotate_tfbs(
      types = names(iDMCs()),
      types_mot = meme_out_dir(),
      meme_out_fold = meme_out_fold,
      meme_bin_path = meme_bin_path,
      motif_dbs = motif_dbs,
      tf_db = tf_db,
      progress_ind = progress
    )}
  )
  
  
  ############## outputs laoding #############
  output$back_sum_info <- renderText({
    if (isTruthy(backFile()))
    {
      paste0("Loaded ", nrow(backFile()), " background CpGs")
    } else {
      ""
    }
  })
  
  output$targ_sum_info <- renderUI({
    if (isTruthy(targFile()))
    {
      ncpgs <- sapply(targFile(), nrow)
      HTML(paste0(names(targFile()) , " loaded ",  ncpgs , " target CpGs <br/> "))
    } else {
      HTML("<br/>")
    }
  })
  
  
  output$Ann_sum_info <- renderText({
    if (isTruthy(hg19_loc()))
    {
      paste0("Loaded ", nrow(hg19_loc()), " CpG annotations ")
    } else {
      ""
    }
  })

  
  ############## outputs motif summary #############
  
  output$motifs_summary <- DT::renderDataTable(meme_out_dir(),
                                               options = list(scrollX =
                                                                TRUE))
  output$exportSum <- downloadHandler(
    filename = "MotifSummary.csv",
    contentType = "csv",
    content = function(file) {
      write.csv(meme_out_dir(), file)
    }
  )

  
  output$distHeat <- renderPlot({
    plot <- NULL
    
    if(!is.null(mot_similiarity()))
      plot <- heatmap.2(
      as.matrix(as.dist(mot_similiarity())),
      scale = "none",
      margins = c(12, 12),
      cexCol = 0.9,
      cexRow = 0.9
    )
    
    plot
  })
  
  output$distDend <- renderPlot({
    
    if(!is.null(mot_similiarity())){
      hc_all = hclust(as.dist(mot_similiarity()))
      
      col_dend  <- hc_all %>% as.dendrogram
      dir <-
        as.factor(sapply(strsplit(colnames(
          mot_similiarity()
        ), "-->"), function(x) {
          x[3]
        }))
      dir <- factor(dir, levels = c("neut","hypo","hyper"))
      print(as.dist(mot_similiarity()))
      print(colnames(mot_similiarity()))
      print(dir)
      par(mar = c(15, 5, 2, 2))
      plot(col_dend)
      colored_bars(as.character(c("gray", "green", "red")[dir]), col_dend)
    }
    #add methylation info as a colored bar
  })
  
  
  ############## outputs tfbs summary #############
  
  output$motifs_summary_tfbs <- DT::renderDataTable(motif_tfbs(),
                                                    options = list(scrollX =
                                                                     TRUE))
  output$exportSumTfbs <-
    downloadHandler(
      filename = "MotifSummaryTfbs.csv",
      contentType = "csv",
      content = function(file) {
        write.csv(motif_tfbs(), file)
      }
    )
  
  output$tfbsHeat <- renderPlot({
    ann_tfbs <- motif_tfbs()
    
    rows <- unique(ann_tfbs$tfbs2)
    columns <- unique(ann_tfbs$type)
    
    #create a binary matrix where we mark the presence of the TFs in the experimnts
    mat = matrix(
      0,
      nrow = length(rows),
      ncol = length(columns),
      dimnames = list(rows, columns)
    )
    for (i in 1:nrow(ann_tfbs)) {
      if(ann_tfbs[i,]$main_dir == "neut"){
        mat[ann_tfbs$tfbs2[i], ann_tfbs$type[i]] = 1
      } else if(ann_tfbs[i,]$main_dir == "hypo"){
        mat[ann_tfbs$tfbs2[i], ann_tfbs$type[i]] = 2
      } else if(ann_tfbs[i,]$main_dir == "hyper"){
        mat[ann_tfbs$tfbs2[i], ann_tfbs$type[i]] = 3
      }
    }
    
    heatmap.2(
      t(mat),
      Colv = T,
      Rowv = T,
      dendrogram = "both",
      density.info = "none",
      trace = "none",
      colsep = 1:(nrow(mat) + 1),
      sepwidth = c(0.001, 0.001),
      rowsep = 1:(ncol(mat) + 1),
      sepcolor = "black",
      col = c("white", "gray", "green","red"),
      cexCol = 0.9,
      cexRow = 0.9,
      srtCol = 45,
      key = F
    )
    
    legend("topleft",     
           legend = c("not sig","neutral","hypo","hyper"),
           fill = c("white", "gray", "green","red")
    )
  })
  
}

shinyApp(ui, server)