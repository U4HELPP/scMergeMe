# scMergeMe: A Single Cell RNA-Seq Data Merger App

**Author:** Amy Olex  
**Date:** June 6, 2025  
**Generative AI Statement:** This application was designed in collaboration with large language models (Gemini 2.5 Flash and Claude Opus 4.6) to assist in bioinformatics software engineering tasks.

## Click here for [USER Documentation Website](https://u4helpp.github.io/scMergeMe/)

---

## Table of Contents

- [Purpose of the Application](#purpose-of-the-application)
- [Features](#features)
- [Installation](#installation)
   - [Prerequisites](#prerequisites)
   - [R Package Dependencies](#r-package-dependencies)
   - [Configuration Files](#configuration-files)
- [Usage](#usage)
   - [Running the Application](#running-the-application)
   - [Step 1. Select Samples Tab](#step-1-select-samples-tab)
   - [Step 2. Configure Merge Tab](#step-2-configure-merge-tab)
   - [Step 3. Run Merge Tab](#step-3-run-merge-tab)
- [Key File Structure](#key-file-structure)
- [Troubleshooting / Known Issues](#troubleshooting--known-issues)
- [Contributing](#contributing)
- [License](#license)

---

## Purpose of the Application

This Shiny R application provides a user-friendly graphical interface to streamline the process of merging single-cell RNA-seq datasets. It simplifies the setup and submission of complex Seurat merge jobs on a SLURM-managed cluster, allowing users to select samples from a master list, configure various bioinformatics parameters, and generate/submit the necessary scripts without direct command-line interaction for each step.

The application automates:
* Loading and displaying available samples from a master list.
* Managing selected samples for a merge operation.
* Validation of input file paths and user parameters.
* Generation of a formatted sample sheet CSV.
* Generation of a SLURM batch submission script (`.sh`).
* Submission of the SLURM job to the cluster.

---

## Features

* **Interactive Sample Selection:** Browse and select multiple samples from a pre-defined master list using a searchable picker interface.
* **Dynamic Sample Management:** Add and remove samples from the merge list on the fly.
* **Comprehensive Parameter Configuration:** Adjust key parameters for the Seurat merge script, including:
    * Run ID and Output Directory
    * SLURM resource allocation (Memory, Cores)
    * Normalization and Merge Type (LogNormalize/SCT, simple/integration)
    * Parallelization options
    * Species for cell cycle scoring
    * Output file format (RData/h5Seurat)
    * Advanced filtering and regression options (Cells2Keep, Ambient RNA Adjustment, Downsampling, Cell Cycle Regression, UMI Regression)
    * Optional input files (features lists, barcode filters)
* **Pre-submission Validation:** The app performs checks on required inputs and file paths before job submission.
* **Automated Script Generation:** Automatically creates a structured sample sheet and an executable SLURM submission script based on user inputs.
* **SLURM Job Submission:** Integrates directly with `sbatch` to submit the merge job to your cluster.
* **Real-time Feedback:** Provides status updates and job submission output directly within the app.

---

## Installation

### Prerequisites

* **R (>= 4.0):** The core programming language.
* **RStudio (Recommended):** An integrated development environment for R.
* **SLURM Workload Manager:** This application is designed to submit jobs to a SLURM cluster. Ensure `sbatch` and other SLURM commands are available in your environment.
* **Git:** For cloning the repository.

### R Package Dependencies

Install the necessary R packages using the following commands in your R console:

```R
install.packages(c("shiny", "shinydashboard", "DT", "shinyWidgets", "readr", "stringr", "dplyr"))
```

### Configuration Files

This application relies on two external files at specific paths. You **must** ensure these files exist and are correctly formatted, or update their paths within the `app.R` script if they are located elsewhere in your environment.

1.  **`Master_Sample_List.csv`**
    * **Default Path:** `/lustre/home/harrell_lab/scRNASeq/config_slurm/Master_Sample_List.csv`
    * **Purpose:** Contains metadata for all available single-cell RNA-seq samples.
    * **Required Columns:**
        * `SampleName`: Unique identifier for each sample. Includes the PDXSource ID and the mouse number.
        * `DataType`: Type of data (e.g., "10X", "Seurat").
        * `SamplePath`: Absolute path to the sample data (10X folder or .h5Seurat file).
        * `PDXSource`: The base PDX identifier without a mouse number.
        * `MouseDepletion`: Indicates if mouse depletion was performed for the sample.
        * `AcquiredResistance`: Indicates if the PDX has any aquired resistance to therapies.
        * `Sex`: Male or Female.
        * `PrimaryCancerType`: The primary cancer types the PDX was derived from.
        * `MetastaticSampleTissue`: The tissue type a metastasis was removed from. The primary cancer type can be breast while the metastatic tissue can be lung.
        * `Pipeline`: Barnyard or Human Only.  Barnyard pipeline include aligning to the human/mouse merged reference and first performing species deconvolution.  Human Only are samples directly aligned to the human reference genome with no species deconvolution.
        * `BatchID`: The ID from the sequencing facility.
        * `Cells2Keep`: (Optional) Path to a file specifying cells to keep for filtering.
        * `raw10Xdata`: (Optional) Path to raw 10X data for ambient RNA adjustment.
        * `RunSoupX`: (Optional) Indicator for SoupX processing. If set to 1 AND AmbientRNA  adjustment is selected, then only these samples will be adjusted.
    * **Note:** If your CSV has different column names or is located elsewhere, update `MASTER_SAMPLE_LIST_PATH` and the `display_cols` and `select` statements in `server.R` accordingly.

2.  **`SeuratMerge_100322.R`**
    * **Default Path:** `/lustre/home/harrell_lab/src/WCCTR_RNASeq_Pipeline/SingleCell/SeuratMerge_100322.R`
    * **Purpose:** This is the backend R script that performs the actual Seurat merging operations. The Shiny app generates a sample sheet and SLURM script to call this R script.
    * **Note:** Ensure this script is executable and its path is correctly set in `MERGE_SCRIPT_PATH` in `app.R`.

---

## Usage

### Running the Application

1.  Clone this repository to your local machine:
    ```bash
    git clone [repository-url]
    cd [repository-directory]
    ```
2.  Open the `app.R` file in RStudio.
3.  Click the "Run App" button in RStudio, or execute from your R console:
    ```R
    shiny::runApp()
    ```

### Step 1) Select Samples Tab

This is your starting point for choosing the samples you wish to merge.

* **Available Samples:**
    * The large table at the top displays all samples found in your `Master_Sample_List.csv`.
    * Use the "Select samples from the list:" picker input to search for and select samples. You can select multiple samples.
    * Click "Add Selected Samples to Merge List" to move them to the lower table.
* **Samples to be Merged:**
    * This table shows the samples currently slated for the merge job.
    * Each row has a "Remove" button to deselect a sample from the merge list.
    * The "Clear All Selected Samples" button will empty the entire merge list and return all samples to the "Available Samples" picker.
* **Sample Selection Status:** Provides a quick count of selected samples and prompts if more are needed.

### Step 2) Configure Merge Tab

In this tab, you define all the parameters for your Seurat merge job.

* **Basic Run Details:**
    * `Unique Run ID`: A short, descriptive identifier for your merge run (e.g., `MyProjectMerge`). This will be incorporated into output file names and SLURM job names.
    * `Output Directory`: The absolute path where all outputs (sample sheet, SLURM script, logs, merged Seurat object) will be saved. You can use `~/` for your home directory.
    * `Memory Limit (GB) for SLURM Job`: Requested RAM for your SLURM job.
    * `Number of Cores for R Script and SLURM Job`: Number of CPU cores requested for the job.
* **Merge & Normalization Options:**
    * `Normalization Type`: Choose between `LogNormalize` and `SCT` (Seurat's SCTransform).
    * `Merge Type`: Select `simple` (standard merge) or `integration` (Seurat's integration workflow for batch correction).
    * `Enable Parallelization`: Check this box to enable parallel processing within the R script (recommended for performance).
    * `Species for Cell Cycle Scoring`: Select `human` or `mouse` for cell cycle regression.
    * `Save Merged Seurat Object as .h5Seurat`: Choose to save the final Seurat object in `.h5Seurat` format (compatible with `SeuratDisk`) instead of the default `.RData`.
* **Filtering & Regression Options:**
    * `Filter to only selected cells`: If checked, the app will use the paths specified in the `Cells2Keep` column of your sample sheet to filter cells.
    * `Adjust for Ambient RNA Contamination`: Check to enable SoupX-based ambient RNA adjustment for those samples with a 1 in the RunSoupX column. Requires `RunSoupX` and `raw10Xdata` columns in your sample sheet.
    * `Downsample Percentage`: Downsample each sample to a percentage of its cells (e.g., 50 for 50%).
    * `Number of Anchor Genes`: Relevant for `integration` merge type.
    * `Regress out Cell Cycle Differences`: Apply regression for cell cycle effects.
    * `Regress out UMI Counts`: Apply regression for UMI count differences.
* **Optional Input Files:**
    * Paths to external files for `Features List`, `Barcodes to Exclude`, or `Barcodes to Keep`. Leave blank if not needed.

### Step 3) Run Merge Tab

This tab provides a final summary and allows you to submit your job.

* **Review & Submit:**
    * A summary of all your selected samples and configured parameters is displayed. Review this carefully.
    * **Click "Submit Merge Job to SLURM"**:
        * The app will perform final validation checks.
        * It will then generate the sample sheet (`[RunID]_[YYMMDD]_sample_sheet.csv`) in your specified output directory.
        * It will generate the SLURM submission script (`06_slurm_[RunID]_[Species]_[YYMMDD].sh`) in the same directory.
        * Finally, it will execute the `sbatch` command to submit the job to the cluster.
* **Job Submission Output:** Displays the output from the `sbatch` command, including your job ID.
* **Monitoring Instructions:** Provides commands to monitor your job's status (`squeue`) and view its output/error logs (`tail -f`).

---

## Key File Structure

.

├── app.R                       # The main Shiny application code

├── Master_Sample_List.csv      # (External) Your master list of samples

└── SeuratMerge_100322.R        # (External) The backend R script for merging

**Note:** `Master_Sample_List.csv` and `SeuratMerge_100322.R` are expected at specific absolute paths in the current configuration. Adjust `MASTER_SAMPLE_LIST_PATH` and `MERGE_SCRIPT_PATH` in `app.R` if your setup differs.

---

## Troubleshooting / Known Issues

* **"object 'dplyr::filter' not found" or similar for `%>%`:** Ensure `library(dplyr)` is included at the top of your `app.R` file.
* **"object 'SampleName' not found" or other column errors:** Double-check that the column names in your `Master_Sample_List.csv` exactly match what's expected by the app (case-sensitive!).
* **`showNotification` error about `type` argument:** The `type` argument for `showNotification` must be one of "default", "message", "warning", or "error". Change any `type = "success"` to `type = "message"`.

---

## Contributing

Feel free to open issues or submit pull requests for any improvements or bug fixes.

---

## License

This application is licensed under the **GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007**.

You can find the full text of the GNU General Public License v3 at:  
[https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html)
