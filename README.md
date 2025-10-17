# 🍄 Magic Mushroom Proteomics

**Bioinformatics workflow for extracting psilocybin biosynthetic gene cluster (BGC) protein sequences from fungal genomic DNA.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

This workflow replicates the protein extraction pipeline from [Bradshaw et al. (2024)](https://doi.org/10.1073/pnas.2311245121) for ancestral sequence reconstruction of psilocybin biosynthesis genes in *Psilocybe* mushrooms.

---

## 🎯 Purpose

Extract clean, intron-free protein sequences for the four core psilocybin biosynthetic genes:
- **PsiD** - Tryptophan decarboxylase (tryptophan → tryptamine)
- **PsiK** - Kinase (phosphorylation)
- **PsiM** - Methyltransferase (methylation)
- **PsiH** - P450 monooxygenase (C-4 hydroxylation, critical for psychoactivity)

## 📊 Background

Based on the phylogenomics study of psychoactive *Psilocybe* mushrooms by Bradshaw et al. (2024), which revealed:
- Psilocybin biosynthesis arose ~67 million years ago
- Two distinct gene cluster organizations (Clade I vs II)
- Evidence of horizontal gene transfer to other mushroom genera
- 71 specimens analyzed (23 type specimens from museum collections)

## 🚀 Quick Start

### 1. Clone and Setup
```bash
git clone https://github.com/YOUR_USERNAME/magic-mushroom-proteomics.git
cd magic-mushroom-proteomics

# Create conda environment
conda env create -f environment.yml
conda activate AncSeq
```

### 2. Run the Workflow
```bash
jupyter lab
# Open: 01_extract_protein_sequences.ipynb
# Run cells sequentially
```

### 3. Verify Results
The notebook will extract protein sequences and compare them to validated examples.

Expected output:
```
✅ PROTEIN: Perfect match!
✅ CDS: Perfect match!
✅ GFF: Perfect match! (or minor date differences)
```

---

## 📁 Repository Structure

```
magic-mushroom-proteomics/
├── 01_extract_protein_sequences.ipynb  # Main workflow notebook
├── psilocybin_utils.py                 # Python utility functions
├── environment.yml                     # Conda environment specification
├── results/                            # Output directory (gitignored)
└── README.md                           # This file
```

## 📓 Notebooks

### `01_extract_protein_sequences.ipynb`
Step-by-step workflow for protein extraction:
1. **Check input files** - Verify reference proteins and genome scaffolds
2. **Run Exonerate** - Extract coding sequences (CDS) using protein2genome alignment
3. **Translate** - Convert CDS to protein using EMBOSS transeq
4. **Quality Control** - Check for stop codons and length validation
5. **Verify** - Compare outputs against validated examples

## 🧰 Utility Functions (`psilocybin_utils.py`)

- **`ExonerateRunner`** - Wrapper for exonerate protein2genome extraction
- **`translate_sequence()`** - EMBOSS transeq wrapper
- **`quality_control()`** - QC checks (stop codons, length)
- **`run_mafft()`** - Multiple sequence alignment
- **`run_iqtree()`** - Phylogenetic tree construction
- **`combine_fasta_files()`** - Merge FASTA files for batch processing

## 🔬 Workflow Details

### Pipeline Overview
```
Genomic DNA (with introns)
         ↓
    [Exonerate]  ← Reference Protein
         ↓
   CDS (no introns)
         ↓
    [Transeq]
         ↓
  Protein Sequence
         ↓
   [Quality Control]
         ↓
  ✅ Clean Protein
```

### Key Parameters
- **Exonerate**: `--model protein2genome --bestn 1`
- **Transeq**: `-frame 1`
- **QC**: No internal stop codons, length validation

---

## 🧬 Psilocybin Biosynthetic Gene Cluster

### The Four Core Enzymes

| Gene | Enzyme | Function | Color Code |
|------|--------|----------|------------|
| **PsiD** | Tryptophan decarboxylase | Tryptophan → Tryptamine | 🔴 Red |
| **PsiK** | Kinase | Phosphorylation | 🟢 Green |
| **PsiM** | Methyltransferase | Methylation | 🟡 Gold |
| **PsiH** | P450 monooxygenase | C-4 hydroxylation (critical!) | 🟣 Purple |

### Gene Cluster Organization

**Clade I:** `PsiD → PsiK → PsiH → PsiM`
**Clade II:** `PsiD → PsiM → PsiH → PsiK` *(canonical)*

Two distinct arrangements corresponding to an ancient phylogenetic split ~67 million years ago!

---

## 🛠️ Installation

### Requirements
- Python 3.9+
- Conda or Mamba
- Exonerate 2.4.0
- EMBOSS 6.6.0
- MAFFT 7.526
- IQ-TREE 3.0.1

### Create Environment
```bash
conda create -n AncSeq python=3.9
conda activate AncSeq
conda install -c bioconda -c conda-forge exonerate emboss mafft iqtree
pip install biopython jupyter pandas matplotlib seaborn
```

Or use the provided `environment.yml`:
```bash
conda env create -f environment.yml
```

---

## 📚 References

**Primary Publication:**
- Bradshaw AJ, Ramírez-Cruz V, Awan AR, Furci G, Guzmán-Dávalos L, Dentinger BTM (2024). Phylogenomics of the psychoactive mushroom genus *Psilocybe* and evolution of the psilocybin biosynthetic gene cluster. *PNAS* 121(3):e2311245121. [DOI: 10.1073/pnas.2311245121](https://doi.org/10.1073/pnas.2311245121)

**Psilocybin Biosynthesis:**
- Fricke J, Blei F, Hoffmeister D (2017). Enzymatic synthesis of psilocybin. *Angew. Chem. Int. Ed.* 56:12352–12355.

---

## 📝 License

MIT License - See LICENSE file for details

## 🤝 Contributing

Contributions welcome! Please open an issue or pull request.

## 📧 Contact

For questions about this workflow or ancestral sequence reconstruction, please open an issue.

---

**⚠️ Legal Notice:** This repository is for research purposes only. Psilocybin is a controlled substance in many jurisdictions. This workflow is designed for academic study of fungal evolution and biosynthesis.
