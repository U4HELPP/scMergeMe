# Global Variables and Libraries (to be placed at the top of app.R or in global.R)
library(shiny)
library(shinydashboard) # For a more structured UI
library(DT)           # For interactive tables
library(shinyWidgets) # For enhanced select input (pickerInput)
library(readr)        # For efficient CSV reading
library(stringr)      # For string manipulation
library(dplyr)

# --- Paths ---
MASTER_SAMPLE_LIST_PATH <- "/lustre/home/harrell_lab/scRNASeq/config_slurm/Master_Sample_List.csv"
MERGE_SCRIPT_PATH <- "/lustre/home/harrell_lab/src/WCCTR_RNASeq_Pipeline/SingleCell/SeuratMerge_100322.R"

# --- Helper Function: Get current date in YYMMDD format ---
get_yymmdd <- function() {
  format(Sys.Date(), "%y%m%d")
}

# --- Helper Function: Check for required 10X files ---
# This function will be called during validation.
# It checks if a given path is a directory and contains the 3 core 10X files.
check_10x_output_dir <- function(path) {
  if (!dir.exists(path)) {
    return(paste0("Path does not exist or is not a directory: ", path))
  }
  required_files <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
  # Handle cases where features.tsv.gz might be genes.tsv.gz
  if (!file.exists(file.path(path, "features.tsv.gz")) && !file.exists(file.path(path, "genes.tsv.gz"))) {
    return(paste0("Missing 'features.tsv.gz' or 'genes.tsv.gz' in ", path))
  }
  for (file in required_files[1:2]) { # Check matrix and barcodes
    if (!file.exists(file.path(path, file))) {
      return(paste0("Missing '", file, "' in ", path))
    }
  }
  return(TRUE) # All checks passed
}


ui <- dashboardPage(
  dashboardHeader(title = "scMergeMe: Self-Service scRNASeq Sample Merge App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Select Samples", tabName = "sample_selection", icon = icon("mouse-pointer")),
      menuItem("2. Configure Merge", tabName = "merge_config", icon = icon("cogs")),
      menuItem("3. Run Merge", tabName = "run_merge", icon = icon("play-circle"))
    )
  ),
  dashboardBody(
    tabItems(
      # --- Tab 1: Sample Selection ---
      tabItem(
        tabName = "sample_selection",
        h2("1. Select Samples for Merging"),
        fluidRow(
          box(
            title = "Available Samples",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            # Use pickerInput for multi-selection with search
            pickerInput(
              inputId = "available_samples_picker",
              label = "Select samples from the list:",
              choices = NULL, # Populated by server
              multiple = TRUE,
              options = list(
                `actions-box` = TRUE,
                `live-search` = TRUE,
                size = 10,
                `selected-text-format` = "count > 3"
              )
            ),
            actionButton("add_selected_samples", "Add Selected Samples to Merge List", icon = icon("arrow-down")),
            hr(),
            # Display all available samples for context (optional, but good for large lists)
            h4("Full List of Available Samples (for reference):"),
            DTOutput("full_sample_list_dt")
          ),
          box(
            title = "Samples to be Merged",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            helpText("Samples added here will be included in the merge. Click 'Remove' to deselect."),
            DTOutput("selected_samples_dt"),
            actionButton("clear_selected_samples", "Clear All Selected Samples", icon = icon("trash"))
          )
        ),
        fluidRow(
          box(
            title = "Sample Selection Status",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("sample_selection_status")
          )
        )
      ),
      
      # --- Tab 2: Configure Merge Parameters ---
      tabItem(
        tabName = "merge_config",
        h2("2. Configure Merge Parameters"),
        fluidRow(
          box(
            title = "Basic Run Details",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            textInput("run_id_input", "Unique Run ID (e.g., MyProjectMerge)", value = ""),
            helpText(textOutput("run_id_preview")), # Display generated run ID
            textInput("output_dir_input", "Output Directory (e.g., ~/my_merge_results)", value = "~/", placeholder = "~/"),
            numericInput("mem_limit_input", "Memory Limit (GB) for SLURM Job", value = 200, min = 1, max = 1000),
            numericInput("num_cores_input", "Number of Cores for R Script and SLURM Job", value = 12, min = 1, max = 64) # Default as per example
          ),
          box(
            title = "Merge & Normalization Options",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            selectInput("normalization_type", "Normalization Type",
                        choices = c("LogNormalize", "SCT"), selected = "LogNormalize"), # Default changed
            selectInput("merge_type", "Merge Type",
                        choices = c("simple", "integration"), selected = "simple"),
            checkboxInput("parallel_checkbox", "Enable Parallelization", value = TRUE),
            selectInput("species_select", "Species for Cell Cycle Scoring",
                        choices = c("human", "mouse"), selected = "human"),
            checkboxInput("save_h5_checkbox", "Save Merged Seurat Object as .h5Seurat (instead of .RData)", value = FALSE)
          )
        ),
        fluidRow(
          box(
            title = "Filtering & Regression Options",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            checkboxInput("filter_cells_checkbox", "Filter to only selected cells (requires Cells2Keep column in sample sheet)", value = TRUE),
            checkboxInput("ambient_rna_adjust_checkbox", "Adjust for Ambient RNA Contamination (requires RunSoupX and raw10Xdata columns, only adjusts samples where RunSoupX is 1)", value = TRUE),
            sliderInput("downsample_percentage", "Downsample Percentage (Keep % of cells from each sample)", min = 0, max = 100, value = 100),
            numericInput("num_anchors_input", "Number of Anchor Genes (for integration merge type)", value = 2000, min = 1),
            checkboxInput("regress_cell_cycle_checkbox", "Regress out Cell Cycle Differences", value = FALSE),
            checkboxInput("regress_umi_checkbox", "Regress out UMI Counts", value = FALSE)
          ),
          box(
            title = "Optional Input Files",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            textInput("features_file_path", "Path to Features List File (Optional)", placeholder = "/path/to/features.txt"),
            textInput("exclude_barcodes_file_path", "Path to Barcodes to Exclude File (Optional)", placeholder = "/path/to/exclude_barcodes.tsv"),
            textInput("keep_barcodes_file_path", "Path to Barcodes to Keep File (Optional)", placeholder = "/path/to/keep_barcodes.tsv")
          )
        )
      ),
      
      # --- Tab 3: Run Merge & Output ---
      tabItem(
        tabName = "run_merge",
        h2("3. Run Merge Job"),
        fluidRow(
          box(
            title = "Review & Submit",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("final_review_summary"),
            actionButton("submit_job", "Submit Merge Job to SLURM", icon = icon("cloud-upload-alt"), class = "btn-success btn-lg"),
            hr(),
            h4("Job Submission Output:"),
            verbatimTextOutput("job_submission_output"),
            h4("Monitoring Instructions:"),
            uiOutput("monitoring_instructions")
          )
        )
      )
    )
  )
)












server <- function(input, output, session) {
  
  # --- Reactive Values ---
  # Stores the full master sample list
  master_sample_data <- reactiveVal(NULL)
  # Stores samples currently selected by the user for merging
  selected_samples_for_merge <- reactiveVal(data.frame())
  
  # --- 1. Load Master Sample List on Startup ---
  observeEvent(TRUE, { # Runs once on app start
    req(file.exists(MASTER_SAMPLE_LIST_PATH))
    tryCatch({
      data <- read_csv(MASTER_SAMPLE_LIST_PATH, show_col_types = FALSE)
      master_sample_data(data)
      
      # Populate pickerInput choices
      # Use SampleName for display, but keep original data for internal use
      updatePickerInput(
        session = session,
        inputId = "available_samples_picker",
        choices = setNames(data$SampleName, data$SampleName) # Both value and label are SampleName for simplicity
      )
    }, error = function(e) {
      showNotification(paste("Error loading master sample list:", e$message), type = "error")
      master_sample_data(NULL) # Set to NULL if loading fails
    })
  }, once = TRUE)
  
  # --- Display Full Sample List (for reference) ---
  output$full_sample_list_dt <- renderDT({
    req(master_sample_data())
    # Only display relevant columns
    display_cols <- c("PDXSource", "SampleName", "Sex", "PrimaryCancerType", "AcquiredResistance", "MetastaticSampleTissue", "RunSoupX")
    DT::datatable(
      master_sample_data()[, display_cols],
      options = list(pageLength = 5, scrollX = TRUE),
      selection = 'none',
      rownames = FALSE
    )
  }, server = FALSE) # Use server=FALSE for smaller datasets, otherwise server=TRUE
  
  # --- 2. Sample Selection Logic ---
  observeEvent(input$add_selected_samples, {
    req(input$available_samples_picker) # Ensure some samples are selected in the picker
    current_selected_names <- selected_samples_for_merge()$SampleName # Names already in the merge list
    
    # Filter master data to get details of newly selected samples
    new_selections <- master_sample_data() %>%
      filter(SampleName %in% input$available_samples_picker, !(SampleName %in% current_selected_names))
    
    if (nrow(new_selections) > 0) {
      updated_selected <- bind_rows(selected_samples_for_merge(), new_selections)
      selected_samples_for_merge(updated_selected)
      
      # Remove added samples from the available picker list
      remaining_choices <- master_sample_data() %>%
        filter(!(SampleName %in% updated_selected$SampleName))
      updatePickerInput(
        session = session,
        inputId = "available_samples_picker",
        choices = setNames(remaining_choices$SampleName, remaining_choices$SampleName),
        selected = character(0) # Clear selected items in picker
      )
    }
  })
  
  # Render selected samples in a DT table
  output$selected_samples_dt <- renderDT({
    req(selected_samples_for_merge())
    if (nrow(selected_samples_for_merge()) == 0) {
      return(NULL)
    }
    # Display selected columns for review, add a 'Remove' button column
    display_cols <- c("PDXSource", "SampleName", "Sex", "PrimaryCancerType", "AcquiredResistance", "MetastaticSampleTissue", "RunSoupX")
    df <- selected_samples_for_merge()[, display_cols, drop = FALSE]
    # Add a 'Remove' button column
    df$Remove <- paste0('<button class="btn btn-danger btn-sm remove_btn" data-id="', df$SampleName, '">Remove</button>')
    DT::datatable(
      df,
      escape = FALSE, # Allow HTML in 'Remove' button
      options = list(dom = 't', # Only show table, no search/pagination
                     columnDefs = list(list(className = 'dt-center', targets = '_all'))),
      selection = 'none',
      rownames = FALSE
    )
  }, server = FALSE) # Use server=FALSE for reactivity with buttons
  
  # Observe clicks on 'Remove' buttons in the selected samples table
  observeEvent(input$selected_samples_dt_cell_clicked, {
    info <- input$selected_samples_dt_cell_clicked
    if (is.null(info$value)) return() # No button clicked
    if (str_detect(info$value, "remove_btn")) {
      sample_name_to_remove <- str_extract(info$value, 'data-id="([^"]+)"', group = 1)
      current_selected <- selected_samples_for_merge()
      updated_selected <- current_selected %>%
        filter(SampleName != sample_name_to_remove)
      selected_samples_for_merge(updated_selected)
      
      # Add the removed sample back to the available picker list
      master_data <- master_sample_data()
      remaining_choices <- master_data %>%
        filter(!(SampleName %in% updated_selected$SampleName))
      updatePickerInput(
        session = session,
        inputId = "available_samples_picker",
        choices = setNames(remaining_choices$SampleName, remaining_choices$SampleName),
        selected = character(0) # Clear selected items in picker
      )
    }
  })
  
  # Clear all selected samples
  observeEvent(input$clear_selected_samples, {
    selected_samples_for_merge(data.frame())
    # Restore all samples to the available picker list
    all_master_samples <- master_sample_data()
    updatePickerInput(
      session = session,
      inputId = "available_samples_picker",
      choices = setNames(all_master_samples$SampleName, all_master_samples$SampleName),
      selected = character(0) # Clear selected items in picker
    )
  })
  
  # Status message for sample selection
  output$sample_selection_status <- renderText({
    n_selected <- nrow(selected_samples_for_merge())
    if (n_selected == 0) {
      "Please select at least two samples to merge."
    } else if (n_selected == 1) {
      "You have selected 1 sample. Please select at least two for merging."
    } else {
      paste0("Ready to merge ", n_selected, " samples.")
    }
  })
  
  # --- 3. Parameter Input Logic ---
  # Reactive for Unique Run ID preview
  output$run_id_preview <- renderText({
    req(input$run_id_input)
    paste0("Final Run ID will be: ", input$run_id_input, "_", get_yymmdd())
  })
  
  # --- 4. Validation and Submission ---
  # Reactive expression to gather all parameters and validate
  all_params <- reactive({
    params <- list(
      runid_user = input$run_id_input,
      outdir = paste0(input$output_dir_input,"/",input$run_id_input),
      mem_limit_gb = input$mem_limit_input,
      num_cores = input$num_cores_input,
      normalization = input$normalization_type,
      merge_type = input$merge_type,
      parallel = input$parallel_checkbox,
      species = input$species_select,
      saveH5 = input$save_h5_checkbox,
      filter = input$filter_cells_checkbox,
      ambientRNAadjust = input$ambient_rna_adjust_checkbox,
      downsample = input$downsample_percentage,
      numAnchors = input$num_anchors_input,
      regressCellCycle = input$regress_cell_cycle_checkbox,
      regressUMI = input$regress_umi_checkbox,
      features_file = input$features_file_path,
      exclude_file = input$exclude_barcodes_file_path,
      keep_file = input$keep_barcodes_file_path
    )
    params
  })
  
  # Final review summary before submission
  output$final_review_summary <- renderText({
    params <- all_params()
    n_samples <- nrow(selected_samples_for_merge())
    if (n_samples < 2) {
      return("Please select at least two samples in the 'Select Samples' tab before reviewing.")
    }
    if (is.null(params$runid_user) || params$runid_user == "") {
      return("Please enter a Unique Run ID in the 'Configure Merge' tab.")
    }
    if (is.null(params$outdir) || params$outdir == "") {
      return("Please enter an Output Directory in the 'Configure Merge' tab.")
    }
    
    # Build summary string
    summary_text <- paste0(
      "--- Merge Job Summary ---\n",
      "Number of Samples Selected: ", n_samples, "\n",
      "Proposed Final Run ID: ", params$runid_user, "_", get_yymmdd(), "\n",
      "Output Directory: ", params$outdir, "\n",
      "SLURM Memory Limit: ", params$mem_limit_gb, "G\n",
      "Cores for Job: ", params$num_cores, "\n",
      "Normalization: ", params$normalization, "\n",
      "Merge Type: ", params$merge_type, "\n",
      "Parallelization: ", params$parallel, "\n",
      "Species: ", params$species, "\n",
      "Save as .h5Seurat: ", params$saveH5, "\n",
      "Filter Cells: ", params$filter, "\n",
      "Ambient RNA Adjustment: ", params$ambientRNAadjust, "\n",
      "Downsample Percentage: ", params$downsample, "%\n",
      "Number of Anchors: ", params$numAnchors, "\n",
      "Regress Cell Cycle: ", params$regressCellCycle, "\n",
      "Regress UMI: ", params$regressUMI, "\n",
      "Optional Features File: ", ifelse(params$features_file == "", "None", params$features_file), "\n",
      "Optional Exclude Barcodes File: ", ifelse(params$exclude_file == "", "None", params$exclude_file), "\n",
      "Optional Keep Barcodes File: ", ifelse(params$keep_file == "", "None", params$keep_file), "\n",
      "\n--- Ready to Submit! ---"
    )
    summary_text
  })
  
  
  # --- Submit Job Button Logic ---
  observeEvent(input$submit_job, {
    # 1. Basic Parameter Validation
    params <- all_params()
    n_samples <- nrow(selected_samples_for_merge())
    if (n_samples < 2) {
      showNotification("Error: Please select at least two samples for merging.", type = "error")
      return()
    }
    if (is.null(params$runid_user) || params$runid_user == "") {
      showNotification("Error: Unique Run ID is required.", type = "error")
      return()
    }
    if (is.null(params$outdir) || params$outdir == "") {
      showNotification("Error: Output Directory is required.", type = "error")
      return()
    }
    
    # Expand tilde for output directory
    final_outdir <- path.expand(params$outdir)
    if (!dir.exists(final_outdir)) {
      tryCatch({
        dir.create(final_outdir, recursive = TRUE)
        showNotification(paste0("Created output directory: ", final_outdir), type = "message")
      }, error = function(e) {
        showNotification(paste0("Error creating output directory: ", e$message), type = "error")
        return()
      })
    }
    
    # 2. File Existence Validation for Samples
    validation_errors <- c()
    for (i in 1:nrow(selected_samples_for_merge())) {
      sample_row <- selected_samples_for_merge()[i, ]
      sample_name <- sample_row$SampleName
      sample_path <- sample_row$SamplePath
      data_type <- sample_row$DataType
      
      if (data_type == "10X") {
        check_result <- check_10x_output_dir(sample_path)
        if (check_result != TRUE) {
          validation_errors <- c(validation_errors, paste0("Sample '", sample_name, "' (10X): ", check_result))
        }
      } else if (data_type == "Seurat") {
        if (!file.exists(sample_path) || !str_detect(sample_path, "\\.h5Seurat$")) {
          validation_errors <- c(validation_errors, paste0("Sample '", sample_name, "' (Seurat): Path is not an existing .h5Seurat file: ", sample_path))
        }
      } else {
        validation_errors <- c(validation_errors, paste0("Sample '", sample_name, "': Unknown DataType '", data_type, "'."))
      }
      
      # Validate Cells2Keep if it exists in the sample sheet (even if filter option is off, it's there)
      if (!is.na(sample_row$Cells2Keep) && sample_row$Cells2Keep != "") {
        if (!file.exists(sample_row$Cells2Keep)) {
          validation_errors <- c(validation_errors, paste0("Sample '", sample_name, "': Cells2Keep file not found: ", sample_row$Cells2Keep))
        }
      }
      # Validate raw10Xdata if it exists in the sample sheet
      if (!is.na(sample_row$raw10Xdata) && sample_row$raw10Xdata != "") {
        if (!file.exists(sample_row$raw10Xdata)) {
          validation_errors <- c(validation_errors, paste0("Sample '", sample_name, "': raw10Xdata file not found: ", sample_row$raw10Xdata))
        }
      }
    }
    
    # Validate optional input files from UI
    if (params$features_file != "" && !file.exists(params$features_file)) {
      validation_errors <- c(validation_errors, paste0("Optional Features File not found: ", params$features_file))
    }
    if (params$exclude_file != "" && !file.exists(params$exclude_file)) {
      validation_errors <- c(validation_errors, paste0("Optional Exclude Barcodes File not found: ", params$exclude_file))
    }
    if (params$keep_file != "" && !file.exists(params$keep_file)) {
      validation_errors <- c(validation_errors, paste0("Optional Keep Barcodes File not found: ", params$keep_file))
    }
    
    if (length(validation_errors) > 0) {
      showNotification(
        HTML(paste("Validation Errors:", paste(validation_errors, collapse = "<br>"))),
        type = "error", duration = NULL
      )
      return()
    }
    
    # 3. Generate Sample Sheet
    sample_sheet_filename <- paste0(params$runid_user, "_", get_yymmdd(), "_sample_sheet.csv")
    sample_sheet_path <- file.path(final_outdir, sample_sheet_filename)
    
    # Prepare data frame for sample sheet
    # Ensure column names match exactly and values are correctly mapped
    sample_sheet_df <- selected_samples_for_merge() %>%
      select(
        SampleName,
        DataType,
        SamplePath,
        PDXSource,
        MouseDepletion,
        AcquiredResistance,
        Sex,
        PrimaryCancerType,
        MetastaticSampleTissue,
        Pipeline, 
        BatchID,
        Cells2Keep,
        raw10Xdata,
        RunSoupX
      )
    write_csv(sample_sheet_df, sample_sheet_path)
    showNotification(paste0("Sample sheet generated at: ", sample_sheet_path), type = "message")
    
    # 4. Generate SLURM Bash Script
    slurm_script_filename <- paste0("06_slurm_", params$runid_user, "_", params$species, "_", get_yymmdd(), ".sh")
    slurm_script_path <- file.path(final_outdir, slurm_script_filename)
    
    # Build Rscript arguments
    rscript_args <- c(
      paste0("-r ", params$runid_user, "_", get_yymmdd()),
      paste0("-o ", final_outdir),
      paste0("-c ", sample_sheet_path),
      paste0("-i ", params$normalization),
      paste0("-t ", params$merge_type)
    )
    if (params$parallel) { rscript_args <- c(rscript_args, "--parallel") }
    rscript_args <- c(rscript_args, paste0("-n ", params$num_cores)) # Use the user-defined cores
    if (params$filter) { rscript_args <- c(rscript_args, "--filter") }
    if (params$ambientRNAadjust) { rscript_args <- c(rscript_args, "--ambientRNAadjust") }
    if (params$downsample != 100) { rscript_args <- c(rscript_args, paste0("-s ", params$downsample)) }
    if (params$numAnchors != 2000) { rscript_args <- c(rscript_args, paste0("-a ", params$numAnchors)) }
    if (params$species != "human") { rscript_args <- c(rscript_args, paste0("-p ", params$species)) }
    if (params$saveH5) { rscript_args <- c(rscript_args, "--saveH5") }
    if (params$regressCellCycle) { rscript_args <- c(rscript_args, "--regressCellCycle") }
    if (params$regressUMI) { rscript_args <- c(rscript_args, "--regressUMI") }
    if (params$features_file != "") { rscript_args <- c(rscript_args, paste0("-f ", params$features_file)) }
    if (params$exclude_file != "") { rscript_args <- c(rscript_args, paste0("-e ", params$exclude_file)) }
    if (params$keep_file != "") { rscript_args <- c(rscript_args, paste0("-k ", params$keep_file)) }
    
    # Construct bash script content
    bash_script_content <- c(
      "#!/bin/bash",
      paste0("#SBATCH --job-name Merge_", params$runid_user), # Make job name more specific
      paste0("#SBATCH --output ", file.path(final_outdir, paste0("06_", params$runid_user, "_", get_yymmdd(), "_output.log"))),
      paste0("#SBATCH --error ", file.path(final_outdir, paste0("06_", params$runid_user, "_", get_yymmdd(), "_error.log"))),
      paste0("#SBATCH --cpus-per-task ", params$num_cores),
      paste0("#SBATCH --mem ", params$mem_limit_gb, "G"),
      "",
      "module load R/4.4.1",
      "module load hdf5/1.14.3",
      "",
      paste0("echo \"Starting Seurat Merge job for ", params$runid_user, "...\""),
      paste0("Rscript ", MERGE_SCRIPT_PATH, " ", paste(rscript_args, collapse = " \\\n  ")),
      "",
      paste0("echo \"Seurat Merge job for ", params$runid_user, " completed.\"")
    )
    writeLines(bash_script_content, slurm_script_path)
    Sys.chmod(slurm_script_path, mode = "0755") # Make executable
    showNotification(paste0("SLURM script generated at: ", slurm_script_path), type = "message")
    
    # 5. Submit Job
    submit_command <- paste0("sbatch ", slurm_script_path)
    submission_result <- system(submit_command, intern = TRUE, ignore.stderr = FALSE) # Capture output
    
    output$job_submission_output <- renderText({
      if (length(submission_result) > 0) {
        paste(submission_result, collapse = "\n")
      } else {
        "No output from sbatch command."
      }
    })
    
    # Extract Job ID (example: "Submitted batch job 12345")
    job_id_match <- str_extract(submission_result[1], "Submitted batch job (\\d+)")
    job_id <- ifelse(!is.na(job_id_match), str_extract(job_id_match, "\\d+"), "N/A")
    
    output$monitoring_instructions <- renderUI({
      div(
        p(HTML(paste0("Your SLURM job has been submitted!"))),
        p(HTML(paste0("<strong>Job ID:</strong> ", job_id))),
        p(HTML(paste0("You can check job status in your terminal using: <code>squeue -u $USER</code> (replace $USER with your username)"))),
        p(HTML(paste0("Output Log File: <code>", file.path(final_outdir, paste0("06_", params$runid_user, "_", get_yymmdd(), "_output.log")), "</code>"))),
        p(HTML(paste0("Error Log File: <code>", file.path(final_outdir, paste0("06_", params$runid_user, "_", get_yymmdd(), "_error.log")), "</code>"))),
        p("You can monitor the logs using `tail -f` (e.g., `tail -f <output_log_file>`)."),
        p("The Shiny App does not need to remain open to monitor the job progress on the cluster.")
      )
    })
    
    showNotification("Job submitted successfully!", type = "message")
    
  }) # End observeEvent submit_job
  
} # End server function

# Run the application
# shinyApp(ui = ui, server = server) # This line would be at the very end of app.R