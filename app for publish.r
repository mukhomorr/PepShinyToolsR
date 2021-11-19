library(shiny)
library(marker)
library(ampir)
library(shinythemes)
library(stringr)
library(shinydashboard)
library(Peptides)
library(protViz)
library(Biostrings)
library(stringi)
library(tools)
library(msa)
library(shiny)
library(Biostrings)
library(BiocManager)

options(repos = BiocManager::repositories())

wd <- ""
setwd(wd)
###########
###########
all_aa <- c("A","R","N","D","C","E","G","H","I","L","K",
            "M","F","P","O","S","U","T","W","Y","V","B","Z","X","J")
nonpolar <- c("G","V","A","C","P","K","R","H","L","I","M","V","F","W")
polar <- c("S","T","Y","N","Q","D","E")
Arg_C <- paste0("R",all_aa)[-14]
Asp_N <- paste0(all_aa,"D")
CNBr <- paste0("M",all_aa)
Chymotrypsin <- c(paste0("F",all_aa)[-14],paste0("L",all_aa)[-14],paste0("W",all_aa)[-14],paste0("Y",all_aa)[-14])
Glu_C_bicarbonate <- paste0("E",all_aa)[-14]
Glu_C_phosphate <- c(paste0("E",all_aa)[-14],paste0("D",all_aa)[-14])
Lys_C <- paste0("K",all_aa)
Pepsin_pH_1.3 <- c(paste0("F",all_aa),paste0("L",all_aa))
Pepsin_pH_more2 <- c(paste0("F",all_aa),paste0("L",all_aa),paste0("W",all_aa),paste0("Y",all_aa),
                     paste0("A",all_aa),paste0("E",all_aa),paste0("Q",all_aa))
Proteinase_K <- c(paste0("A",all_aa),paste0("C",all_aa),paste0("G",all_aa),paste0("M",all_aa),
                  paste0("F",all_aa),paste0("S",all_aa),paste0("Y",all_aa),paste0("W",all_aa))
Trypsin <- Arg_C <- c(paste0("R",all_aa)[-14],paste0("K",all_aa)[-14])


example_path <-  system.file("examples", "exampleAA.fasta", package="msa")
##############
noPTM <- 0
cys <- 57.0215
glu <- 0.984
asp <- 0.984
meth <- 15.9949
ac <- 42.0106
qq <- -17.0265
sp <- 105.0578
ptm_frag2 <- function(seq,cysteins, glutamine, aspargine,methionine,sp){
  frag_spec <- as.data.frame(fragmentIon(seq))
  nCist <- unlist(lapply(frag_spec$fragment, function(X){
    newSeq <- unlist(strsplit(X,''))
    nCist <- ifelse(test = 'C' %in% newSeq, yes = table(newSeq)[['C']], no = 0)
  }))
  nGlut <- unlist(lapply(frag_spec$fragment, function(X){
    newSeq <- unlist(strsplit(X,''))
    nGlut <- ifelse(test = 'Q' %in% newSeq, yes = table(newSeq)[['Q']], no = 0)
  }))
  nMet <- unlist(lapply(frag_spec$fragment, function(X){
    newSeq <- unlist(strsplit(X,''))
    nMet <- ifelse(test = 'M' %in% newSeq, yes = table(newSeq)[['M']], no = 0)
  }))
  nAsp <- unlist(lapply(frag_spec$fragment, function(X){
    newSeq <- unlist(strsplit(X,''))
    nAsp <- ifelse(test = 'N' %in% newSeq, yes = table(newSeq)[['N']], no = 0)
  }))
  #frag_spec$mass <-   frag_spec$mass + nCist * cysteins+nGlut*glutamine +nAsp*aspargine +nMet*methionine+nCist*sp
  frag_spec$mass <- - nCist*57.0215  + frag_spec$mass + nCist*cysteins+nGlut*glutamine +nAsp*aspargine +nMet*methionine+nCist*sp
  return(frag_spec)  
}
mz_extended <-
  function(seq, charge=2, label="none", aaShift = NULL, cysteins, glutamine, aspargine,methionine, acetylation,qq,sp){
    
    # Check for correct input
    if(!is.numeric(charge) | length(charge) != 1){
      stop("Charge must be given as an integer (typically between 1-4).")
    }
    
    # Calculate mass of uncharged peptide
    mass <- mw(seq=seq, label=label, aaShift=aaShift, monoisotopic = TRUE)
    
    # Add modification at cysteins
    # DO modified the call to the str_count function and vectorized it using native R functions. 
    #mass <- mass + str_count(seq, 'C') * cysteins
    #Cysteins
    nCist <- unlist(lapply(seq, function(X){
      newSeq <- unlist(strsplit(X,''))
      nCist <- ifelse(test = 'C' %in% newSeq, yes = table(newSeq)[['C']], no = 0)
    }))
    #Glutamines
    nGlut <- unlist(lapply(seq, function(X){
      newSeq <- unlist(strsplit(X,''))
      nGlut <- ifelse(test = 'Q' %in% newSeq, yes = table(newSeq)[['Q']], no = 0)
    }))
    #####aspargines
    nAsp <- unlist(lapply(seq, function(X){
      newSeq <- unlist(strsplit(X,''))
      nAsp <- ifelse(test = 'N' %in% newSeq, yes = table(newSeq)[['N']], no = 0)
    }))
    ##### Methionine
    nMet <- unlist(lapply(seq, function(X){
      newSeq <- unlist(strsplit(X,''))
      nMet <- ifelse(test = 'M' %in% newSeq, yes = table(newSeq)[['M']], no = 0)
    }))
    ##### Q -> q
    nQq <- unlist(lapply(seq, function(X){
      newSeq <- unlist(strsplit(X,''))
      nQq <- as.numeric(ifelse(test = 'Q' %in% newSeq[1], yes = table(newSeq[1])[['Q']], no = 0))
      
    }))
    mass <- mass + nCist * cysteins + nGlut * glutamine + nAsp*aspargine +nMet*methionine +acetylation + nQq*qq+nCist*sp
      # Modify for charged peptides.
      if (charge >= 0){
        mass <- mass + charge * 1.007276 # Add weight of H+ ions.
        mass <- mass / charge # Divide by charge state.
      }
      
      return(mass)
    }
######info
type <- c("Tiny","Small","Aliphatic","Aromatic","Non-polar","Polar","Charged","Basic","Acidic")
acids <- c("A + C + G + S + T","A + B + C + D + G + N + P + S + T + V","A + I + L + V","F + H + W + Y","A + C + F + G + I + L + M + P + V + W + Y",
           "D + E + H + K + N + Q + R + S + T + Z","B + D + E + H + K + R + Z","H + K + R","B + D + E + Z")
types <- as.data.frame(cbind(type,acids))

####UI part
ui <- navbarPage("Protpep tools",
                 theme = shinytheme("united"), ### chose theme
  tabPanel("Sequence marker", 
           #img(src = "E:/Rdirect/pep_marker/Examples/amino-acid-structures_med.png"),
               # tags$style(type='text/css', '#test {white-space: pre-wrap;}'),
  mainPanel(h3(textOutput("fasta_nam")),
    h4(id = "text-to-mark",textOutput("test"),align = "wide")), #textOutput("test") ####your sequence for marks viewed here
  
  
  use_marker(),
  tags$head(
    tags$style(type='text/css', '#test {width:800px; word-wrap:break-word;line-height:26pt;}',
      ".red{background-color:#FFB8C3;}.blue{background-color:#6ECFEA;}.orange{background-color:#FFA500;}.lime{background-color:#6495ED;}"
    )
  ),
  #img(src = "aa.png",height = 100, width = 200),
  sidebarPanel(
  fileInput(inputId = "datafile", label = "Upload a data file. Should be .fasta or .txt with first line containing >name. This line will not be highlited" , multiple = FALSE,
            placeholder = "No file selected", accept = "txt"),
            
   #                   actionButton("change", "Load seq", class = "btn btn-success"),
  #actionButton("example", "Load example seq"),
  actionButton(inputId = "update", label = "Calc pI & mw"),
  selectInput('columns', 'Select sequence', ""),
  h2(textOutput("pI")),h2(textOutput("mw")),
  h2(textOutput("x")),
  tableOutput("summary"),

  #tableOutput("class"),
  h4("Mark digestion sites:"),
  actionButton(inputId = "trypsin", label = "Trypsin",class = "btn btn-warning"),
  actionButton(inputId = "arg_c", label = "Arg C"),
  actionButton(inputId = "asp_n", label = "Asp N"),
  actionButton(inputId = "chymotrypsin", label = "Chymotrypsin"),
  actionButton(inputId = "cnbr", label = "CNBr"),
  actionButton(inputId = "glu_c_bi", label = "Glu C (bicarbonate)"),
  actionButton(inputId = "glu_c_ph", label = "Glu C (phosphate)"),
  actionButton(inputId = "lys_c", label = "Lys C"),
  actionButton(inputId = "pepsin13", label = "Pepsin at pH 1.3"),
  actionButton(inputId = "pepsin2", label = "Pepsin at pH > 2 "),
  actionButton(inputId = "proteinase_K", label = "Proteinase K"),
  actionButton(inputId = "clear", label = "No digestion",class = "btn btn-danger"),
  
  
  h4("Mark polar or nonpolar AAs:"),
  actionButton(inputId = "polar", label = "Polar",class = "btn btn-info"),
  actionButton(inputId = "nonpolar", label = "Nonpolar",class = "btn btn-warning"),
  actionButton(inputId = "clear_pol", label = "unmark polar"),
  actionButton(inputId = "clear_nonpol", label = "Unmark nonpolar"),
  h4("Mark specific AAs:"),
  actionButton(inputId = "clear_aa", label = "Unmark AA",class = "btn btn-danger"),
   checkboxGroupInput("icons", "Choose amino acids to highlight:",
                     choiceNames =
                       list("Ala 	A 	Alanine",
                            "Arg 	R 	Arginine",
                            "Asn 	N 	Asparagine",
                            "Asp 	D 	Aspartic acid",
                            "Cys 	C 	Cysteine",
                            "Glu 	E 	Glutamic acid",
                            "Gly 	G 	Glycine",
                            "His 	H 	Histidine",
                            "Ile 	I 	Isoleucine",
                            "Leu 	L 	Leucine",
                            "Lys 	K 	Lysine",
                            "Met 	M 	Methionine",
                            "Phe 	F 	Phenylalanine",
                            "Pro 	P 	Proline",
                            "Pyl 	O 	Pyrrolysine",
                            "Ser 	S 	Serine",
                            "Sec 	U 	Selenocysteine",
                            "Thr 	T 	Threonine",
                            "Trp 	W 	Tryptophan",
                            "Tyr 	Y 	Tyrosine",
                            "Val 	V 	Valine",
                            "Asx 	B 	Aspartic acid or Asparagine",
                            "Glx 	Z 	Glutamic acid or Glutamine",
                            "Xaa 	X 	Any amino acid",
                            "Xle 	J 	Leucine or Isoleucine"
                       ),
                     choiceValues =
                       list("A","R","N","D","C","E","G","H","I","L","K",
                            "M","F","P","O","S","U","T","W","Y","V","B","Z","X","J")
                     
  ))),
  ########
  tabPanel("Peptide Summary",
           mainPanel(  
            
           h4(tableOutput("pep_table")),
           tableOutput("pep_summary")
             ),
           sidebarPanel(
             textInput("peptide","Insert peptide sequence"),
  actionButton(inputId = "pep", label = "Calc", class = "btn btn-success"),


  ##### input stuff for PTM
  #tableOutput("classs"),
  selectInput(inputId = "cysteins", label = "Carbamidomethylation", choices = c(noPTM,cys)),
  selectInput(inputId = "sp", label = "S-pyridylethylation", choices = c(noPTM,sp)),
  selectInput(inputId = "glu_asp", label = "Deamidation", choices = c(noPTM,glu)),
  selectInput(inputId = "methionine", label = "Methionine Oxidation", choices = c(noPTM,meth)),
  selectInput(inputId = "acetylation", label = "Acetylation", choices = c(noPTM,ac)),
  selectInput(inputId = "qq", label = "Pyro-glu from Q", choices = c(noPTM,qq)),
  textOutput("mz_1"),textOutput("mz_2"),textOutput("mz_3"),textOutput("mz_4"),textOutput("mz_5")
  ###########

  #####chose PTM
  #checkboxGroupInput("PTMs", "Choose PTMs:",
   #                  choiceNames =
    #                   list("cysteins","glutamine","aspargine","methionine","acetylation"
     #                  ),
      #               choiceValues =
       #                list(57.02,0.98,0.98,15.99,42.01)
                     
  #),
 )),
 tabPanel("Peptide fragmentation",
          mainPanel(
          h4("CID"),
          tableOutput("fragmentation_by"),
          h4("ETD"),
          tableOutput("fragmentation_zc")
 ),
 sidebarPanel(
   textInput("peptide2","Insert peptide sequence"),
   actionButton(inputId = "frag", label = "frag", class = "btn btn-success"),
   selectInput(inputId = "cysteins1", label = "Carbamidomethylation", choices = c(noPTM,cys)),
   selectInput(inputId = "sp1", label = "S-pyridylethylation", choices = c(noPTM,sp)),
   selectInput(inputId = "glu_asp1", label = "Deamidation", choices = c(noPTM,glu)),
   selectInput(inputId = "methionine1", label = "Methionine Oxidation", choices = c(noPTM,meth))
   #selectInput(inputId = "acetylation1", label = "Acetylation", choices = c(noPTM,ac)),
   #selectInput(inputId = "qq1", label = "Pyro-glu from Q", choices = c(noPTM,qq)),
 )),
 tabPanel("genome sequence translation",
          fileInput(inputId = "datafile_fasta", label = "Upload a data file", multiple = FALSE,
                    placeholder = "No file selected", accept = "fasta"),
          actionButton("translation", "Load fasta and translate"),
          downloadButton("downloadData_3", "Download 3 frames"),
          downloadButton("downloadData_3_2", "Download 3 reversed frames"),
          downloadButton("downloadData_6", "Download 6 frames"),
          downloadButton("downloadData1", "Download first frame"),
          downloadButton("downloadData2", "Download second frame"),
          downloadButton("downloadData3", "Download third frame"),
          downloadButton("downloadData4", "Download fourth frame"),
          downloadButton("downloadData5", "Download fifth frame"),
          downloadButton("downloadData6", "Download sixth frame"),
          textOutput("info")
 ),
 tabPanel("Amino acid sequence allignment",
          h3("1. Chose files bellow"),
          h3("2. Press load button for every sequence"),
          h3("3. Press allign. From this point allignment, TXT and PDF generating starts, wait a few seconds"),
          h3("If pdf does not appear and error message is displayed - one of your sequence is too long"),
          h3("4. Check your working folder"),
          
          fileInput(inputId = "datafile_fasta1", label = "Upload first seq", multiple = FALSE,
                    placeholder = "No file selected", accept = "fasta"),
          fileInput(inputId = "datafile_fasta2", label = "Upload second seq", multiple = FALSE,
                    placeholder = "No file selected", accept = "fasta"),
          actionButton("load1", "Load first seq"),
          actionButton("load2", "Load second seq"),
          actionButton("allign", "allign sequences"),
          h3("Those sequences are going to be alligned"),
          textOutput("info1"),
          textOutput("info2")),
 tabPanel("Remove * from fasta",fileInput(inputId = "zvezda", label = "Upload fasta", multiple = FALSE,
                                          placeholder = "No file selected", accept = "fasta"),
          actionButton("remove", "Load first seq"),
          downloadButton("download_zvezda", "Download fasta with no stop codons"),
          textOutput("zx")
          )
 )
                  


server <- function(input, output,session) {
  #setwd(wd)
  contentsrea <- reactive({
    if (input$example){
      inFile <- input$datafile
      input$datafile <- example_path
      readChar(inFile$datapath,nchars = 10000)
    }
    else{
    #File inpet
    inFile <- input$datafile
    if (is.null(inFile))
      return(NULL)
    readChar(inFile$datapath,nchars = 10000)
  }
  })
 
  #### Show list of seqs
  outVar = reactive({
    req(input$datafile$datapath)
    data <- read_faa(input$datafile$datapath)
    rownames(data) <- data$seq_name
    rownames(data)
  })
  observe({
    updateSelectInput(session, "columns",
                      choices = outVar()
    )})
  
  #### Print example sequence
  title_change <- reactive({
    input$columns
    
    req(input$datafile$datapath) ## your app wont quit if there is no file
    data <- read_faa(input$datafile$datapath) ####your list path
    rownames(data) <- data$seq_name
    data[input$columns,2]
  })

  
     ## chose list of sequences
 
    

  observeEvent(  ### printing your sequence
    input$columns, { output$test <- renderText({
      req(input$datafile$datapath) ## your app wont quit if there is no file
      data <- read_faa(input$datafile$datapath) ####your list path
      rownames(data) <- data$seq_name
      data[input$columns,2]
      
    })
    
    req(input$datafile$datapath) ## your app wont quit if there is no file
    fasta_name <- input$columns
    output$fasta_nam <- renderText({fasta_name})
    })
  
  
  #output$some_text <- renderText({
  #highlight(input$datafile, input$search)

  

  observeEvent(input$update,{
    req(input$datafile$datapath) ## your app wont quit if there is no file
    data <- read_faa(input$datafile$datapath) ####your list path
    rownames(data) <- data$seq_name
    out <- data[input$columns,2]
    output$pI <- renderText({  
      
      c("pI =",as.character(round(pI(out),3)))})})
  
  observeEvent(input$update,{
    req(input$datafile$datapath) ## your app wont quit if there is no file
    data <- read_faa(input$datafile$datapath) ####your list path
    rownames(data) <- data$seq_name
    out <- data[input$columns,2]
    output$mw <- renderText({
   
    c("mw = ",as.character(mw(out,monoisotopic = T)),"Da")})})
  
  observeEvent(input$update,{
    req(input$datafile$datapath) ## your app wont quit if there is no file
    data <- read_faa(input$datafile$datapath) ####your list path
    rownames(data) <- data$seq_name
    out <- data[input$columns,2]
    summary <- as.data.frame(aaComp(out))
    summary <- cbind(summary,types$acids)
    colnames(summary)[1:3] <- c("Number","%","acids")
    output$summary <- renderTable(summary,rownames = T)
                                               })
  #####Peptide calc part
  observeEvent(input$pep,{pep_pI <- round(pI(str_remove_all(input$peptide," ")),3)
  pep_hp <- round(hydrophobicity(str_remove_all(input$peptide," ")),3)
  mw_pep <- mw(str_remove_all(input$peptide," "),monoisotopic = T)
  #output$pep_summary <- renderTable(aaComp(str_remove_all(as.character(input$peptide)," ")),rownames = T)
  summary <- as.data.frame(aaComp(str_remove_all(as.character(input$peptide)," ")))
  summary <- cbind(summary,types$acids)
  colnames(summary)[1:3] <- c("Number","%","acids")
  output$pep_summary <- renderTable(summary,rownames = T)

  



    mz_1 <- mz_extended(input$peptide,1,"none",NULL,as.numeric(input$cysteins),as.numeric(input$glu_asp),
                                                                 as.numeric(input$glu_asp),as.numeric(input$methionine),
                                                                 as.numeric(input$acetylation),as.numeric(input$qq),as.numeric(input$sp))
    mz_2 <- mz_extended(input$peptide,2,"none",NULL,
                                                                 as.numeric(input$cysteins),as.numeric(input$glu_asp),
                                                                 as.numeric(input$glu_asp),as.numeric(input$methionine),as.numeric(input$acetylation),as.numeric(input$qq),as.numeric(input$sp))
    mz_3 <- mz_extended(input$peptide,3,"none",NULL,as.numeric(input$cysteins),
                                                                 as.numeric(input$glu_asp),
                                                                 as.numeric(input$glu_asp),as.numeric(input$methionine),
                                                                 as.numeric(input$acetylation),as.numeric(input$qq),as.numeric(input$sp))
    mz_4 <- mz_extended(input$peptide,4,"none",NULL,as.numeric(input$cysteins),
                                                                 as.numeric(input$glu_asp),
                                                                 as.numeric(input$glu_asp),as.numeric(input$methionine),
                                                                 as.numeric(input$acetylation),as.numeric(input$qq),as.numeric(input$sp))
    mz_5 <- mz_extended(input$peptide,5,"none",NULL,as.numeric(input$cysteins),
                                                                 as.numeric(input$glu_asp),
                                                                 as.numeric(input$glu_asp),as.numeric(input$methionine),
                                                                 as.numeric(input$acetylation),as.numeric(input$qq),as.numeric(input$sp))
    
names_vec <- c("pI","Hydrophobicity","mw","1+","2+","3+","4+","5+") 
value <- c(pep_pI,pep_hp,round(mw_pep,4),round(mz_1,4),round(mz_2,4),round(mz_3,4),round(mz_4,4),round(mz_5,4))
pep_table <- cbind(names_vec,value)
rownames(pep_table) <- names_vec
                   output$pep_table <- renderTable(pep_table[,-1],rownames = T,colnames = F)
    
    
    })
  
  observeEvent(input$frag,{
    #####frag cacl and output!!!!!!
    
    #a <- as.numeric(input$cys)
    #x <- as.numeric(input$glu_asp)
    #y <- as.numeric(input$glu_asp)
    #z <- as.numeric(input$methionine)
    #w <- as.numeric(input$sp)
    frag_spec <- ptm_frag2(input$peptide2,as.numeric(input$cysteins1),as.numeric(input$glu_asp1),as.numeric(input$glu_asp1),as.numeric(input$methionine1),as.numeric(input$sp1))
    #    frag_spec <- ptm_frag2(input$peptide,200,0,0,0,0)
    by <- frag_spec[which(frag_spec$type == "b" | frag_spec$type == "y"),-5]
    by_ions <- cbind(by,
                     cbind(by$mass-18.01528,by$mass-17.031))
    colnames(by_ions)[5:6] <- c("-H2O","-NH3")
    
    output$fragmentation_by <- 
      renderTable(by_ions[,-c(1,3)],
                  rownames = T)
    zc <- frag_spec[which(frag_spec$type == "z" | frag_spec$type == "c"),-5]
    output$fragmentation_zc <- 
      renderTable(zc[,-c(1,3)],
                  rownames = T)
    
    #output$fragmentation_zc <- 
    #        renderTable(frag_spec,
    #            rownames = T)
  })
  #observeEvent(input$pep,{output$mz_2 <- renderText(c("2+ = ",as.character(mz_extended(input$peptide,2,cysteins = cysteins,
   #                                                                                    glutamine = glutamine,
    #                                                                                   aspargine = aspargine,
     #                                                                                  methionine = methionine,acetylation = acetylation
      #                                                                                 ))))})
  #observeEvent(input$pep,{output$mz_3 <- renderText(c("3+ = ",as.character(mz_extended(input$peptide,3,cysteins = cysteins,
   #                                                                                    glutamine = glutamine,
    #                                                                                   aspargine = aspargine,
     #                                                                                  methionine = methionine,acetylation = acetylation
      #                                                                                 ))))})
  #observeEvent(input$pep,{output$mz_4 <- renderText(c("4+ = ",as.character(mz_extended(input$peptide,4,cysteins = cysteins,
   #                                                                                    glutamine = glutamine,
    #                                                                                   aspargine = aspargine,
     #                                                                                  methionine = methionine,acetylation = acetylation
      #                                                                                 ))))})
  #observeEvent(input$pep,{output$mz_5 <- renderText(c("5+ = ",as.character(mz_extended(input$peptide,5,cysteins = cysteins,
   #                                                                                    glutamine = glutamine,
    #                                                                                   aspargine = aspargine,
     #                                                                                  methionine = methionine,acetylation = acetylation
      #                                                                                ))))})
  ##############
  marker <- marker$new("#text-to-mark")
  
  observeEvent(input$icons, {
    marker$
      unmark(className = "blue")$ # unmark red class
      mark(input$icons, className = "blue") # add red class
  })
  
  observeEvent(input$trypsin, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Trypsin, className = "red") # add red class

    
  })
  observeEvent(input$arg_c, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Arg_C, className = "red") # add red class
    
  })
  observeEvent(input$asp_n, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Asp_N, className = "red") # add red class
    
  })
  observeEvent(input$chymotrypsin, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Chymotrypsin, className = "red") # add red class
    
  })
  observeEvent(input$cnbr, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(CNBr, className = "red") # add red class
    
  })
  observeEvent(input$glu_c_bi, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Glu_C_bicarbonate, className = "red") # add red class
    
  })
  observeEvent(input$glu_c_ph, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Glu_C_phosphate, className = "red") # add red class
    
  })
  observeEvent(input$lys_c, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Lys_C, className = "red") # add red class
    
  })
  observeEvent(input$pepsin13, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Pepsin_pH_1.3, className = "red") # add red class
    
  })
  observeEvent(input$pepsin2, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Pepsin_pH_more2, className = "red") # add red class
    
  })
  observeEvent(input$proteinase_K, { 
    marker$
      unmark(className = "red")$ # unmark red class
      mark(Proteinase_K, className = "red") # add red class
    
  })
  observeEvent(input$clear, { 
    marker$
      unmark(className = c("red")) # unmark red class
       # add red class
  })
  observeEvent(input$clear_aa, { 
    marker$
      unmark(className = c("blue")) # unmark red class
    # add red class
  })
  observeEvent(input$clear_pol, { 
    marker$
      unmark(className = c("lime")) # unmark red class
    # add red class
  })
  observeEvent(input$clear_nonpol, { 
    marker$
      unmark(className = c("orange")) # unmark red class
    # add red class
  })
  observeEvent(input$nonpolar, { 
    marker$
      unmark(className = "orange")$ # unmark red class
      mark(polar, className = "orange") # add red class
    
  })
  observeEvent(input$polar, { 
    marker$
      unmark(className = "lime")$ # unmark red class
      mark(nonpolar, className = "lime") # add red class
    
  })
  
  contentsrea <- reactive({
    #File inpet
    inFile_f <- input$datafile_fasta
    if (is.null(inFile_f))
      return(NULL)
    readDNAStringSet(inFile_f$datapath)
  })
  observeEvent(input$translation, {
    if(!is.null(input$datafile_fasta)){
      x <- readDNAStringSet(as.character(input$datafile_fasta$datapath))
      ######
      ####get normal translation
      ####get normal translation
      x1 <- DNAStringSet(x,start = 1)
      xx1 <- translate(x1,)
      
      stri_sub(xx1@ranges@NAMES,1,2) <- "fr_1_"
      
      x2 <- DNAStringSet(x,start = 2)
      xx2 <- translate(x2)
      stri_sub(xx2@ranges@NAMES,1,2) <- "fr_2_"
      
      x3 <- DNAStringSet(x,start = 3)
      xx3 <- translate(x3)
      stri_sub(xx3@ranges@NAMES,1,2) <- "fr_3_"
      
      ##### Get reversed translation
      y <- reverseComplement(x)
      y1 <- DNAStringSet(y,start = 1)
      yy1 <- translate(y1)
      stri_sub(yy1@ranges@NAMES,1,2) <- "fr_4_"
      
      y2 <- DNAStringSet(y,start = 2)
      yy2 <- translate(y2)
      stri_sub(yy2@ranges@NAMES,1,2) <- "fr_5_"
      
      y3 <- DNAStringSet(y,start = 3)
      yy3 <- translate(y3)
      stri_sub(yy3@ranges@NAMES,1,2) <- "fr_6_"
      
      
      zz3 <- c(xx1,xx2,xx3)
      zz3_2 <- c(yy1,yy2,yy3)
      zz6 <- c(xx1,xx2,xx3,yy1,yy2,yy3)
      #######     
  
      #######     
      
      output$info <- renderText(zz3@ranges@NAMES[1])
      output$downloadData_3 <- downloadHandler(
        filename = function() {
          paste("data-", "_3fr_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(zz3,file)
        }
      )
      output$downloadData_6 <- downloadHandler(
        filename = function() {
          paste("data-","_6fr_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(zz6,file)
        }
      )
      
      output$downloadData_3_2 <- downloadHandler(
        filename = function() {
          paste("data-","_3fr_reversed_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(zz3_2,file)
        }
      )
      output$downloadData1 <- downloadHandler(
        filename = function() {
          paste("data-","_fr1_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(xx1,file)
        }
      ) 
      output$downloadData2 <- downloadHandler(
        filename = function() {
          paste("data-","_fr2_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(xx2,file)
        }
      ) 
      output$downloadData3 <- downloadHandler(
        filename = function() {
          paste("data-","_fr3_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(xx3,file)
        }
      )
      output$downloadData4 <- downloadHandler(
        filename = function() {
          paste("data-","_fr4_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(yy1,file)
        }
      )
      output$downloadData5 <- downloadHandler(
        filename = function() {
          paste("data-","_fr5_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(yy2,file)
        }
      )
      output$downloadData6 <- downloadHandler(
        filename = function() {
          paste("data-","_fr6_", ".fasta", sep="")
        },
        content = function(file) {
          writeXStringSet(yy3,file)
        }
      )
      
    }})

contentsrea <- reactive({
  #File inpet
  inFile1 <- input$datafile_fasta1
  if (is.null(inFile1))
    return(NULL)
  readAAStringSet(inFile1$datapath)
})

contentsrea <- reactive({
  #File inpet
  inFile2 <- input$datafile_fasta2
  if (is.null(inFile2))
    return(NULL)
  readAAStringSet(inFile2$datapath)
})
observeEvent(input$load1, {
  if(!is.null(input$datafile_fasta1)){
    x <<- readAAStringSet(as.character(input$datafile_fasta1$datapath))
    
    output$info1 <- renderText(x@ranges@NAMES)
  }})
observeEvent(input$load2, {
  if(!is.null(input$datafile_fasta2)){
    y <<- readAAStringSet(as.character(input$datafile_fasta2$datapath))
    
    output$info2 <- renderText(y@ranges@NAMES)
    
    
    #######     
    
    
  }})
observeEvent(input$allign, {
  #setwd(wd)
  x <<- readAAStringSet(as.character(input$datafile_fasta1$datapath))
  y <<- readAAStringSet(as.character(input$datafile_fasta2$datapath))
  
  xy <- c(x,y)
  allignment <<- msa(xy)
  z <<- format(Sys.time(), "%Y%m%d%H%M%S")
  options(width = 225)
  capture.output(print(allignment,show = "complete",showConsensus = T),
                 file = paste0(gsub(":","-",z),"allignment.txt"))

  #msaPrettyPrint(allignment, output="tex", showNames="none",
  #               showLogo="none", askForOverwrite=FALSE, verbose=FALSE,file = paste0(gsub(":","-",z),"allignment.tex"))

  #texi2pdf(paste0(gsub(":","-",z),"allignment.tex"), clean=T, quiet = T)

  #######     
  
  
})
contentsrea <- reactive({
  #File inpet
  inFile3 <- input$zvezda
  if (is.null(inFile3))
    return(NULL)
  readChar(inFile3$datapath,nchar = 10000000)
})

observeEvent(input$remove, {
  if(!is.null(input$zvezda)){
    sp <- as.character(readChar(input$zvezda$datapath,nchar = 10000000))
    spp <- strsplit(sp, ">")[1]
    s <- spp[[1]][-1]
    sppp <- str_replace_all(s,"\\*","")
    sppp <- str_replace(sppp,"\r\n","_ABZAC_")
    sppp <- paste0(">",sppp)
    sppp <- str_replace(sppp,">","_ABZAC_>")
    sppp <- str_replace_all(sppp,"\r\n","")
    #head(sppp)
    sppp <- str_replace_all(sppp,"_ABZAC_","\n")
    #head(sppp)
    output$zx <- renderText({"Done!"})
    sppp <- paste0(sppp,sep = "",collapse = "")
    output$download_zvezda <- downloadHandler(
      filename = function() {
        paste("data-","no_sym", ".fasta", sep="")
      },
      content = function(file) {
        writeLines(sppp,file)
      }
    )
  }})

}

shinyApp(ui = ui, server = server)