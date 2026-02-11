# 🧬 MultiomicMenu

**Lead Scientist:** Dr. William G Ryan V
**Last Updated:** February 11, 2026 
**Status:** [Draft]

This repository provides a standardized, configuration-driven framework for integrating high-dimensional multi-modal data (e.g., RNA-seq, Proteomics, etc.). The purpose of utilizing **Quarto (`.qmd`)** in this pipeline is to transform raw data into a fully documented, publication-ready HTML or PDF reports.

---

## 🚀 Quick Start

You do not need to rewrite the analysis code. Follow these steps to generate your report:

1.  **Use the Template, Download, and Navigate to Folder:**
    a. Navigate to [https://github.com/CogDisResLab/MultiomicMenu](https://github.com/CogDisResLab/MultiomicMenu)
    b. Click "Use this template" on the upper right corner of the page.
    c. Choose a name for your repository. Remember to make it private for sensitive data.
    d. Use git clone on your new repository:
    ```bash
    git clone `git clone [https://github.com/](https://github.com/)[username]/[repo-name].git`
    cd [repo-name]
    ```

2.  **Environment Setup:**
    Open R in the project working directory and restore the required package versions:
    ```r
    if (!require("renv")) install.packages("renv")
    renv::restore()
    ```
3.  **Configure:**
    Edit `index.qmd` with the following *important details*:
    a. title
    b. species
    c. gmt (MUST be changed if not using human data)
    d. data -> value -> data files
4.  **Render:**
    Execute the pipeline via the terminal:
    ```bash
    quarto render .
    ```
---

## 📂 Repository Structure

* `index.qmd`: **The Single Source of Truth.** Define all hyperparameters and file paths here.
* `data/`: Put your DEGs, DPPs, DAPs, or DAKs here.
* `learn.qmd`: Nothing is altered here: this is for informational purposes.
* `data.qmd`: All data processing is done in this step, including integration and PPI generation using Kinograte.
* `pathways.qmd`: All pathway data is processed here using PAVER.
* `networks.qmd`: Network diagrams are generated here using igraph.

---

## ⚙️ Configuration (`config.yaml`)

The pipeline logic is controlled entirely by the `config.yaml` file. Common parameters include:

| Parameter | Description | Default/Example |
| :--- | :--- | :--- |
| `title` | Title used in the final report headers. | `"MultiomicMenu Report"` |
| `species` | Species used in the experiment. | `[human, rat, or mouse]` |
| `gmt` | GMT file to use for pathway information. | `gmt: "https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Human_GO_AllPathways_noPFOCR_with_GO_iea_February_03_2026_symbol.gmt"` |
| `data` | Data files to use. | `["kinase_stk: 'data/kinase.csv'"]` |

---

## 🧪 Statistical Framework

The integration primarily utilizes **Prize Collecting Steiner Forest (PCSF)** analysis with a PPI graph.

---

## 🛠 Requirements

* **R:** ≥ 4.5.0
* **Quarto CLI:** System-level installation required ([Download](https://quarto.org/docs/get-started/)).
* **Memory:** 8GB+ RAM required. 16GB+ RAM recommended.

---

## 📊 Outputs

Upon completion, the `results/` directory will contain:
1.  **`index.html`**: An interactive report featuring all of the below pages.
2.  **`learn.html`**: A page dedicated to describing the statistics and bioinformatics of MultiomicMenu.
3.  **`data.html`**: A page for describing the data, including a "prize plot" which shows the cumulative distribution of "prizes".
4.  **`pathways.html`**: A page for pathways that are hits in the analysis.
5.  **`networks.html`**: A concatenation of all networks generated in the analysis, typically including a sample network and a seeded network.

