# scMergeMe: A Single-Cell RNA-Seq Data Merger App

**Author:** Amy Olex

**Date:** April 8, 2026

**Generative AI Statement:** This application was designed in collaboration with large language models (Gemini Pro 2025 and Claude Opus 4.6) to assist in bioinformatics software engineering tasks.

---

## Table of Contents

- [Purpose](#purpose)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Getting Started on Apollo](#getting-started-on-apollo)
- [Using the App](#using-the-app)
  - [Step 1: Settings](#step-1-settings)
  - [Step 2: Select Samples](#step-2-select-samples)
  - [Step 3: Configure Merge](#step-3-configure-merge)
  - [Step 4: Run Merge](#step-4-run-merge)
- [What the App Generates](#what-the-app-generates)
- [Monitoring Your Job](#monitoring-your-job)
- [Configuration Reference](#configuration-reference)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Purpose

scMergeMe provides a graphical interface for merging single-cell RNA-seq datasets on the VCU HPRC Apollo cluster. Instead of manually writing sample sheets and SLURM scripts, you select your samples, configure your parameters, and the app generates everything and submits the job for you.

The app automates loading samples from a master list, validating file paths, generating a formatted sample sheet CSV, generating a SLURM batch script, and submitting the job to the cluster.

---

## Features

- Searchable, interactive sample selection from the master sample list
- Support for simple merge, Harmony integration, and Seurat CCA integration workflows
- LogNormalize and SCTransform normalization options
- Cell filtering, ambient RNA correction (SoupX), downsampling, and cell cycle/UMI regression
- Configurable SLURM resource allocation (memory and cores)
- Pre-submission validation of all file paths and parameters
- Automatic SLURM script generation and job submission
- Configurable paths so the app works with different master sample lists or merge scripts

---

## Prerequisites

- Access to the VCU HPRC Apollo cluster (`apollo.hprc.vcu.edu`)
- Access to the Harrell Lab shared directories on `/lustre`
- An RStudio session running on the cluster (via Open OnDemand or similar)

---

## Getting Started on Apollo

### 1. Start an RStudio Session

Log in to the HPRC Apollo Open OnDemand portal at **apollo.hprc.vcu.edu** and launch an interactive RStudio session. Request enough resources for the RStudio session itself (the actual merge job will be submitted separately to SLURM with its own resource allocation).

### 2. Navigate to the App Directory

In the RStudio **Files** pane (bottom-right), navigate to:

```
/lustre/home/harrell_lab/src/scMergeMe
```

Alternatively, you can set your working directory in the R console:

```r
setwd("/lustre/home/harrell_lab/src/scMergeMe")
```

### 3. Open the App

Click on the `app.R` file in the Files pane to open it in the RStudio editor.

### 4. Run the App

Click the **"Run App"** button in the top-right corner of the editor pane. RStudio will launch the Shiny app in a new window or in the built-in viewer.

### 5. Dismiss the Initial Error

The first time the app launches, you will see an error notification. This happens because the app attempts to load the master sample list from the default path on startup. Simply **dismiss the error by clicking the "x"** on the notification. You will load the sample list properly in the next step using the Settings tab.

---

## Using the App

The app is organized into four tabs in the left sidebar. Work through them in order.

### Step 1: Settings

This tab lets you configure where the app looks for its key input files.

**Path to Master Sample List CSV** — This defaults to the standard lab location (`/lustre/home/harrell_lab/scRNASeq/config_slurm/Master_Sample_List.csv`). If you are working with a different sample list, update this path to point to your own CSV file. Your CSV must follow the same column format as the master list (see [Configuration Reference](#configuration-reference) below).

**Path to Seurat Merge R Script** — This defaults to the lab's merge script (`/lustre/home/harrell_lab/src/WCCTR_RNASeq_Pipeline/SingleCell/SeuratMerge_100322.R`). Only change this if you are using a modified version of the merge script.

Once both paths are set, click the **"Load/Reload Sample List"** button. You should see a green notification confirming the sample list loaded successfully. You can return to this tab and reload at any time if you update your master sample list.

### Step 2: Select Samples

This tab is where you choose which samples to include in the merge.

**Selecting samples:**

1. Use the **"Select samples from the list"** dropdown at the top. It supports live search — start typing a sample name to filter the list. You can also click "Select All" or "Deselect All" using the buttons inside the dropdown.
2. After highlighting your desired samples in the dropdown, click **"Add Selected Samples to Merge List"**. The selected samples will move from the dropdown into the **"Samples to be Merged"** table below.
3. The reference table at the bottom of the "Available Samples" box shows all samples with key metadata columns (PDX source, sex, cancer type, treatment, etc.) so you can look up details before selecting.

**Managing your selection:**

- To remove a single sample, click the red **"Remove"** button in its row in the "Samples to be Merged" table.
- To start over, click **"Clear All Selected Samples"** to return all samples to the available pool.

**Status check:** The "Sample Selection Status" box at the bottom tells you how many samples are selected. You need **at least two samples** to proceed with a merge.

### Step 3: Configure Merge

This tab contains all the parameters that control how the merge is performed and how the SLURM job is configured.

**Basic Run Details:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Unique Run ID | A short name for this merge run (e.g., `TNBC_Comparison`). The app appends the current date in YYMMDD format automatically. | (empty — you must provide one) |
| Output Directory | Where all outputs will be saved. The app creates a subdirectory using your Run ID and date. | `~/` |
| Memory Limit (GB) | RAM requested for the SLURM job. Increase for large merges. | 200 |
| Number of Cores | CPU cores for both the SLURM job and R script parallelization. | 12 |

**Merge and Normalization Options:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Normalization Type | `LogNormalize` (standard log-normalization) or `SCT` (SCTransform). | LogNormalize |
| Merge Type | `simple` (concatenate samples), `harmony` (Harmony integration for batch correction), or `integration` (Seurat CCA-based integration). | simple |
| Enable Parallelization | Run computations in parallel using the specified number of cores. | Enabled |
| Species | `human` or `mouse` — used for cell cycle gene lists during scoring. | human |
| Save as .h5Seurat | Save output in h5Seurat format instead of .RData. | Disabled |

**Filtering and Regression Options:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Filter to only selected cells | Use the `Cells2Keep` column from the sample sheet to subset cells. | Enabled |
| Adjust for Ambient RNA | Run SoupX correction on samples where `RunSoupX` is set to 1. Requires the `raw10Xdata` column. | Enabled |
| Downsample Percentage | Keep only a percentage of cells from each sample (100 = keep all). | 100 |
| Number of Anchor Genes | Number of anchor features for the Seurat `integration` merge type. Not used for `simple` or `harmony`. | 2000 |
| Regress Cell Cycle | Regress out cell cycle score differences during scaling. | Disabled |
| Regress UMI | Regress out UMI count differences during scaling. | Disabled |

**Optional Input Files:**

These are paths to external files and can be left blank if not needed.

- **Features List File** — A text file with a list of genes/features to use instead of the default variable features.
- **Barcodes to Exclude File** — A TSV of cell barcodes to remove from the analysis.
- **Barcodes to Keep File** — A TSV of cell barcodes to retain (all others are removed).

### Step 4: Run Merge

This final tab shows a summary of everything you have configured and lets you submit the job.

**Review your settings:** The summary box displays the number of selected samples, run ID, output directory, and all parameter values. Read through this carefully before submitting.

**Submit the job:** Click the **"Submit Merge Job to SLURM"** button. The app will then:

1. Validate all file paths for the selected samples (checks that 10X directories contain the expected `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` files, or that Seurat `.h5Seurat` files exist).
2. Validate optional file paths if any were provided.
3. Create the output directory if it does not already exist.
4. Generate a sample sheet CSV in the output directory.
5. Generate a SLURM batch script in the output directory.
6. Run `sbatch` to submit the job to the cluster.

If validation fails, you will see error notifications describing exactly which files or paths are missing. Fix the issues and try again.

On success, the app displays the **Job ID** and instructions for monitoring.

---

## What the App Generates

After a successful submission, your output directory will contain:

| File | Description |
|------|-------------|
| `[RunID]_[YYMMDD]_sample_sheet.csv` | The sample sheet CSV passed to the merge script |
| `06_slurm_[RunID]_[Species]_[YYMMDD].sh` | The SLURM batch submission script |
| `06_[RunID]_output.log` | Standard output log from the job (created when job runs) |
| `06_[RunID]_error.log` | Standard error log from the job (created when job runs) |
| Merged Seurat object | The final `.RData` or `.h5Seurat` file (created by the merge script) |

---

## Monitoring Your Job

Once your job is submitted, **you can close the Shiny app** — the job runs independently on the cluster.

### Checking if the Job is Still Running

If you are comfortable with the command line, open a terminal on Apollo and run:

```bash
squeue -u $USER
```

If your job appears in the list, it is still running. If it does not appear, it has either completed or failed.

### Watching the Logs

To watch the output log in real time:

```bash
tail -f /path/to/your/output/06_[RunID]_output.log
```

To watch the error log:

```bash
tail -f /path/to/your/output/06_[RunID]_error.log
```

### Knowing When the Job is Done

The job is complete when the **very last line** of the output log file reads:

```
Seurat Merge job for [RunID] completed.
```

where `[RunID]` is your run ID (e.g., `Seurat Merge job for TNBC_Comparison_250408 completed.`). You can check this at any time by opening the output log file or running:

```bash
tail -1 /path/to/your/output/06_[RunID]_output.log
```

If the job failed or was killed before finishing, this line will not be present. In that case, check the error log for details on what went wrong.

---

## Configuration Reference

### Master Sample List CSV

The master sample list must be a CSV file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `SampleName` | Yes | Unique sample identifier (includes PDX source and mouse number) |
| `DataType` | Yes | `10X` or `Seurat` |
| `SamplePath` | Yes | Absolute path to the 10X output directory or `.h5Seurat` file |
| `PDXSource` | Yes | Base PDX identifier without mouse number |
| `MouseDepletion` | Yes | Whether mouse depletion was performed |
| `AcquiredResistance` | Yes | Whether the PDX has acquired resistance to therapies |
| `Treatment` | Yes | Treatment information for the sample |
| `Sex` | Yes | Male or Female |
| `PrimaryCancerType` | Yes | Cancer type the PDX was derived from |
| `MetastaticSampleTissue` | Yes | Tissue type of metastatic sample |
| `Pipeline` | Yes | `Barnyard` or `Human Only` |
| `CellRangerVersion` | Yes | Cell Ranger version used for alignment |
| `BatchID` | Yes | Sequencing facility batch ID |
| `Cells2Keep` | Optional | Path to file specifying cells to retain |
| `raw10Xdata` | Optional | Path to raw 10X data for SoupX ambient RNA correction |
| `RunSoupX` | Optional | Set to `1` to enable SoupX correction for that sample |

### R Package Dependencies

The app will attempt to install missing packages automatically on startup. The required packages are: `shiny`, `shinydashboard`, `DT`, `shinyWidgets`, `readr`, `stringr`, and `dplyr`.

---

## Troubleshooting

**The app shows an error on first launch** — This is expected. The app tries to load the default master sample list on startup. Go to the Settings tab, verify the path, and click "Load/Reload Sample List."

**"Column not found" errors in the reference table** — Your master sample list CSV may be missing expected columns or have different column names. Column names are case-sensitive and must match exactly (e.g., `SampleName`, not `samplename`).

**Validation errors on submit** — The app checks that every sample's data files exist on disk. If you see path-related errors, verify that the `SamplePath` entries in your master CSV are correct and that the directories/files are accessible from your current session.

**Job submitted but fails immediately** — Check the error log file listed in the submission output. Common causes include the R module not being available (`module load R/4.4.1`), missing R packages on the cluster, or insufficient memory for the number of samples.

**The app is slow to load the sample list** — Large CSV files may take a moment. Be patient after clicking "Load/Reload Sample List."

---

## Contributing

Feel free to open issues or submit pull requests for improvements or bug fixes.

---

## License

This application is licensed under the **GNU General Public License v3.0**.
See [https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html) for details.
