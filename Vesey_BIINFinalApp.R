library(shiny)
library(shinythemes)
library(Biostrings)
library(stringr)

options(shiny.maxRequestSize = 30*1024^2)

# Define UI ----
ui <- fluidPage(
  theme = shinytheme("readable"), 
  titlePanel("Yeast Gene Sequence Search"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta", "Upload your Yeast FASTA or FNA genome file"),
      textInput("gene", "Enter your Nucleic Acid sequence for your gene of interest:"),
      numericInput("max_mismatches", "Max mismatches allowed:", value = 3, min = 0, max = 10),
      checkboxInput("allow_indels", "Allow Insertions and Deletions (indels)?", value = FALSE),
      actionButton("submit", "Submit", class = "btn-primary")  
    ),
    mainPanel(
      htmlOutput("HitsOutput", class = "non-scrolling-box")  
    )
  ),
  tags$div(style = "margin-top: 50px;"),  
  tags$head(
    tags$style(
      HTML(
        "
      .non-scrolling-box {
        max-width: 800px; 
        word-wrap: break-word;
        background-color: #f2f2f2; 
        padding: 10px; 
        border-radius: 5px; 
      }
      .results-text {
        font-size: 20px;
        color: #333; 
        margin-bottom: 10px; 
      }
      .results-count {
        font-size: 16px; 
        color: #555; 
        margin-top: 10px;
      }
      .hit-info {
        font-size: 14px; 
      }
      .btn-primary {
        color: #fff; 
        background-color: #007bff; 
        border-color: #007bff; 
      }
      "
      )
    )
  )
)


# Important Functions --

extract_name <- function(description) {
  name <- str_extract(description, "(chromosome|contig|chr|ctg)\\s*\\S+")
  if (!is.na(name)) {
    location <- gsub("^\\s+|,$", "", name)
    return(location)
  } else {
    return(NULL)
  }
}

check_file_size <- function(file_path, max_size_mb) {
  file_size_mb <- file.info(file_path)$size / 1024^2
  return(file_size_mb <= max_size_mb)
}

# Define server logic ----
server <- function(input, output) {
  hits_data <- reactiveVal(NULL)  
  
  observeEvent(input$submit, {
    if (!grepl("\\.(fasta|fna)$", input$fasta$name, ignore.case = TRUE)) {
      showNotification("File type not accepted.", type = "error")
      return()
    }
    
    if (!grepl("^[ACGTN]+$", input$gene, ignore.case = TRUE)) {
      showNotification("Invalid input: Please enter a valid nucleotide sequence.", type = "error")
      return()
    }
    
    max_file_size_mb <- 30
    if (!check_file_size(input$fasta$datapath, max_file_size_mb)) {
      showNotification(paste("File size exceeds the maximum limit of", max_file_size_mb, "MB."), type = "error")
      return()
    }
    
    # Read FASTA file
    fasta_file <- input$fasta$datapath
    genome_sequences <- tryCatch(
      readDNAStringSet(fasta_file),
      error = function(e) {
        showNotification("Error reading the FASTA file.", type = "error")
        return(NULL)
      }
    )
    
    # Return if genome sequences not read successfully
    if (is.null(genome_sequences)) return()
    
    # Process query sequence
    query_sequence <- DNAString(input$gene)
    
    # Check for query length and max mismatches
    if (input$max_mismatches >= nchar(input$gene)) {
      showNotification("Maximum mismatches cannot be greater than or equal to the length of the input query.", type = "error")
      return()
    }
    
    query_sequence_revcomp <- reverseComplement(query_sequence)
    
    # Lists to store hits and hit sequences
    hits_list <- vector("list", length(genome_sequences))
    hit_sequences <- vector("list", length(genome_sequences))
    
    # Loop through genome sequences
    for (i in seq_along(genome_sequences)) {
      chromosome_sequence <- genome_sequences[[i]]
      
      # Find matches and reverse complement matches
      matches <- matchPattern(query_sequence, chromosome_sequence, max.mismatch = input$max_mismatches, min.mismatch = 0, with.indels = input$allow_indels)
      matches_revcomp <- matchPattern(query_sequence_revcomp, chromosome_sequence, max.mismatch = input$max_mismatches, min.mismatch = 0, with.indels = input$allow_indels)
      
      # Combine matches
      all_matches <- c(matches, matches_revcomp)
      
      # Get start and end positions of matches
      start_positions <- start(all_matches)
      end_positions <- end(all_matches)
      
      # Check for invalid genome positions
      tryCatch({
        valid_indices <- start_positions <= length(chromosome_sequence) & 
          end_positions <= length(chromosome_sequence) & 
          start_positions <= end_positions & 
          start_positions > 0
        
        start_positions <- start_positions[valid_indices]
        end_positions <- end_positions[valid_indices]
        
        if (length(start_positions) > 0) {
          hits_list[[i]] <- data.frame(
            location = extract_name(names(genome_sequences)[i]),
            start_position = start_positions,
            end_position = end_positions,
            sequence = paste(subseq(chromosome_sequence, start_positions, end_positions), collapse = "")
          )
          
          # Extract hit sequences
          hit_seqs <- unlist(subseq(chromosome_sequence, start_positions, end_positions))
          hit_sequences[[i]] <- hit_seqs
        } else {
          hits_list[[i]] <- NULL
          hit_sequences[[i]] <- NULL
        }
      }, error = function(e) {
        showNotification(paste("Error processing hits for search query. Please make sure the supplied 'start', 'end' and 'width' arguments are defining a region that is within the limits of the sequence.", e$message), type = "error")
        hits_list[[i]] <- NULL
        hit_sequences[[i]] <- NULL
      })
    }
    
    # Combine hits and hit sequences
    hits_list <- hits_list[!sapply(hits_list, is.null)]
    output_data <- lapply(seq_along(hits_list), function(i) {
      hit_df <- hits_list[[i]]
      if (nrow(hit_df) == 1) {
        hit_df <- hit_df[1, ]
      } else if (length(hit_sequences[[i]]) > 0) {
        merged_hits <- cbind(hit_df, hit_sequence = paste(hit_sequences[[i]], collapse = ""))
        merged_hits
      } else {
        hit_df
      }
    })
    
    # Update hits data
    hits_data(output_data)  
    
    output$HitsOutput <- renderText({
     
       # Check if any hits found...
      if (length(output_data) == 0 || all(sapply(output_data, is.null))) {
        return("No results found")
      }
      
      # Prepare the output HTML
      output_html <-   # Prepare the output HTML
        output_html <- "<div class='results-text'><strong>Results:</strong>"
      
      # Initialize count of results
      result_count <- 0
      
      # Print hit information
      for (data in output_data) {
        if (!is.null(data)) {
          result_count <- result_count + 1
          
          hit_info <- paste(
            "<div class='hit-info'>",
            "<ul>",
            paste("<li><strong>Location:</strong> ", data$location, "</li>"),
            paste("<li><strong>Start Position:</strong> ", data$start_position, "</li>"),
            paste("<li><strong>End Position:</strong> ", data$end_position, "</li>"),
            paste("<li><strong>Sequence:</strong> ", data$sequence, "</li>"),
            "</ul>",
            "</div>",
            "<hr>",
            sep = ""
          )
          output_html <- paste(output_html, hit_info)
        }
      }
      
      # Add result count above hits
      result_count_html <- paste("<div class='results-count'>", result_count, " results found.</div>")
      output_html <- paste(result_count_html, output_html)
      
      HTML(output_html)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)










