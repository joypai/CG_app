library(shiny)
library(shinythemes)

ui <- shinyUI(fluidPage(theme=shinytheme("cerulean"), title="CG analysis",
  #shinythemes::themeSelector(),
  
  titlePanel("CG analysis"),
  
  # sidebar with parameters
  fluidRow(
    column(width = 2, 
           wellPanel(
             selectInput("patient", "Patient:", c(601,603,605,610,611,613,616)),
             radioButtons("seq_type", "Sequence Type: ", c("full", "env"), inline=TRUE)
           )
    ),
    
    # main plot
    fluidRow(
      column(width=4, h4("CG plot"),
             plotlyOutput("cg_plot", height=500),
             br(),
             p(strong("Hover"), em("over data point to highlight selected sequence in tree")),
             p(strong("Click"), em("on data point to display CG distribution graph"))
      ),
      column(width=5,
             h4("tree"),
             plotOutput("tree_plot", height=600)
      )
    ),
    fluidRow(
      column(width=8, offset=2,
             hr(),
             h4("CG distribution"),
             plotOutput("dist_plot", height=200)
      )
    )
  )
))

#tabPanel("all dinucleotides", plotOutput("di_plot", height = 1000)),
# tabPanel("tree", plotOutput("tree_plot"))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library(Biostrings)
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
  library(plotly)
  library(ggsci)
  library(ggtree)
  library(ape)
  library(phytools)
  
  read_fasta <- reactive({
    patient <- input$patient
    seq_type <- input$seq_type
    if (seq_type == "env") {
      fasta_file = paste("~/data_alpha/home/jpai/CG_dinucleotide/fastas/",patient,"_env.fasta",sep="")
    } else {
      fasta_file = paste("~/data_alpha/home/jpai/CG_dinucleotide/fastas/",patient,"_chinglan_full.fasta",sep="")
    }
    fasta <- readDNAStringSet(fasta_file)
  })
  
  calculate_frequencies <- reactive({
    fasta <- read_fasta()
    dfreq_all <- data.frame(base1=character(0), base2=character(0), frequency=numeric(0), di=character(0), sample=character(0))
    
    # read fasta and calculate dinucleotide frequencies
    for (i in 1:length(fasta)) {
      hiv_seq <- fasta[[i]]
      obs_df <- dinucleotideFrequency(hiv_seq, as.matrix=TRUE)
      dfreq <- as.data.frame(as.table((obs_df/(length(hiv_seq)-1))))
      colnames(dfreq) <- c("base1","base2","frequency")
      dfreq <- dfreq[order( dfreq[,1], dfreq[,2] ),]
      dfreq$di <- paste(dfreq$base1, dfreq$base2,sep="")
      dfreq$sample <- names(fasta[i])
      
      dfreq_all <- rbind(dfreq_all, dfreq)
    }
    
    if (input$seq_type == "env") {
      dfreq_all$group <- ifelse(grepl("env",dfreq_all$sample),"ching-lan",
                                ifelse(grepl("DNA",dfreq_all$sample),"DNA",
                                       ifelse(grepl("POSSC",dfreq_all$sample),"POSSC",
                                              sub("\\d+[_-](\\w+?)[_-].*","\\1",dfreq_all$sample))))
    } else {
      dfreq_all$group <- "chinglan"
    }
    dfreq_all
  })
  
  
  output$cg_plot <- renderPlotly({ 
    dfreq_all <- calculate_frequencies()
    CGfreq_all <- subset(dfreq_all, di=="CG")
    p <- ggplot(CGfreq_all, aes(x=di, y=frequency, label=sample)) + geom_jitter(aes(color=group)) + 
      xlab(NULL) + ylab("frequency")  + theme_bw() + scale_color_aaas()
    ggplotly(p, tooltip = c("y", "label")) %>% layout(margin=list(l = 80), height=500)
  })
  
  
  output$di_plot <- renderPlot({
    dfreq_all <- calculate_frequencies()
    
    ggplot(dfreq_all, aes(x=di, y=frequency)) + geom_point() + facet_wrap(~sample) + 
      xlab("dinucleotide") + ylab("frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 40)) 
  })
  
  
  output$dist_plot <- renderPlot({
    d <- event_data("plotly_click")

    if (is.null(d)) {
      ggplot() + theme_classic()
    } else {
      fasta <- read_fasta()
      dfreq_all <- calculate_frequencies()
      CGfreq_all <- subset(dfreq_all, di=="CG")
      curve <- sort(unique(CGfreq_all$group))[d$curveNumber+1]
      seq_id <- subset(CGfreq_all, group==curve)[d$pointNumber+1,]$sample
      hiv_seq <- fasta[[seq_id]]

      window_size = 200
      increment = 10
      if (input$seq_type != "env") {
        start_offset = 650
        bp_max = 10000
      } else {
        start_offset = 0 
        bp_max=3000
      }
      starts <- seq(1, length(hiv_seq)-window_size, by = increment)
      n <- length(starts)
      dist_df <- data.frame(bp=starts)
      dist_df$CG_count <- sapply(starts, function(x) countPattern("CG",hiv_seq[x:(x+window_size)]))
      dist_df$bp_offset <- dist_df$bp + start_offset
      cg_pos <- (start(matchPattern("CG",hiv_seq))+start_offset)
      p2 <- ggplot(dist_df, aes(x=bp_offset, y=CG_count)) + geom_line() +  coord_cartesian(xlim=c(0,bp_max), ylim=c(0,20), expand=F) + xlab("bp") + ylab("CG count") +
        theme_classic() + geom_point(data=data.frame(bp=cg_pos,y=rep(max(dist_df$CG_count)+3,length(cg_pos))), aes(x=bp,y=y), color="red", shape=124, size=4) + 
        ggtitle(paste("CG distribution: ",seq_id,sep="")) + theme(plot.title = element_text(hjust = 0.5))
      p2
    }
  })
  
  
  output$tree_plot <- renderPlot({
    tree_file <- paste("~/data_alpha/home/jpai/CG_dinucleotide/trees/",input$patient,"_",input$seq_type,"_tree.newick",sep="")
    tree <- read.newick(tree_file)
    midtree <- midpoint.root(tree)
    midtree$tip.label <- gsub("-","_",midtree$tip.label)
    max_x <- max(midtree$edge.length)*4
    data <- data.frame(tip=midtree$tip.label)
    
    if (input$seq_type == "full") {
      data$tip_group <- "ching-lan"
    } else {
      data$tip_group <- ifelse(grepl("env",data$tip),"ching-lan",
                               ifelse(grepl("DNA",data$tip),"DNA",
                                      ifelse(grepl("POSSC",data$tip),"POSSC",
                                             sub("\\d+[_-](\\w+?)[_-].*","\\1",data$tip))))
    }
  
    h <- event_data("plotly_hover")
    dfreq_all <- calculate_frequencies()
    CGfreq_all <- subset(dfreq_all, di=="CG")
    curve <- sort(unique(CGfreq_all$group))[h$curveNumber+1]
    seq_id <- subset(CGfreq_all, group==curve)[h$pointNumber+1,]$sample
  
    if (is.null(h)) {
      tp <- ggtree(midtree) %<+% data + geom_tiplab(aes(color=tip_group), size=3) + geom_treescale(offset=-1) + scale_color_aaas() + xlim(NA, max_x)
    } else {
      data$alpha = ifelse(data$tip==seq_id,1,0.8)

      tp <- ggtree(midtree) + geom_treescale(offset=-1) + xlim(NA, max_x)
      tp <- tp %<+% data + geom_tiplab(aes(color=tip_group, alpha=alpha), size=3) + scale_color_aaas() + scale_alpha_continuous(range=c(0.2,1))
    }  
    tp
  })
})

# Run the application 
shinyApp(ui = ui, server = server)