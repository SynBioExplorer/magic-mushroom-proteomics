# Psilocybe Protein Extraction Pipeline

Extraction and analysis of psilocybin biosynthesis genes (PsiD, PsiK, PsiM, PsiH) from 71 *Psilocybe* species genome assemblies.

## Overview

This pipeline extracts protein-coding sequences for the four core psilocybin biosynthesis enzymes from fungal genome scaffolds using exonerate's protein-to-genome alignment. The output is suitable for AlphaFold structure prediction and phylogenetic analysis.

### Psilocybin Biosynthesis Pathway

| Enzyme | Function | Typical Length |
|--------|----------|----------------|
| **PsiD** | Tryptophan decarboxylase (tryptophan → tryptamine) | ~430-440 aa |
| **PsiK** | Kinase (phosphorylation) | ~350-360 aa |
| **PsiM** | Methyltransferase (methylation) | ~295-310 aa |
| **PsiH** | P450 monooxygenase (C-4 hydroxylation) | ~500-510 aa |

## Method

### Extraction Approach

We use **exonerate** with the `--ryo %tcs` (Roll-Your-Own) output format to extract pre-spliced coding sequences directly. This approach:

1. Aligns reference P. cubensis proteins to target genome scaffolds
2. Identifies exon boundaries accounting for introns
3. Returns the correctly spliced CDS with proper phase handling
4. Avoids manual exon stitching errors

```bash
exonerate --model protein2genome \
  --ryo ">%ti|%s|%tS\n%tcs\n" \
  --bestn 1 --percent 30 \
  reference.faa genome.fasta
```

### Translation

Proteins are translated with `to_stop=True` to produce clean sequences for AlphaFold (no internal stop codons). Internal stops are detected separately for quality control flagging.

## Results Summary

### Extraction Statistics

| Enzyme | Sequences Extracted | Species with Hits | CHECK_ORF | NO_HIT |
|--------|---------------------|-------------------|-----------|--------|
| **PsiD** | 51 | 51 | 0 | 20 |
| **PsiK** | 51 | 51 | 1 | 20 |
| **PsiM** | 45 | 45 | 8 | 26 |
| **PsiH** | 51 | 51 | 1 | 20 |

- **71 total species** analyzed
- **51 species** contain psilocybin biosynthesis genes
- **20 species** lack all four genes (non-psilocybin producers)

### Species with Pseudogenes (CHECK_ORF)

These species have internal stop codons indicating pseudogenes or assembly errors. The output proteins are truncated at the first stop codon.

| Species | Affected Gene | Output Length | Expected Length | Notes |
|---------|---------------|---------------|-----------------|-------|
| P. aztecorum var. bonetii | PsiK | 205 aa | ~360 aa | Truncated |
| P. cubensis | PsiM | 150 aa | ~310 aa | Truncated |
| P. fagicola var. mesocystidiata | PsiM | 39 aa | ~310 aa | Severely truncated |
| P. acutissima | PsiM | 25 aa | ~310 aa | Severely truncated |
| P. hoogshageni | PsiM | 25 aa | ~310 aa | Severely truncated |
| P. subhoogshagenii | PsiM | 25 aa | ~310 aa | Severely truncated |
| P. subyungensis | PsiM | 25 aa | ~310 aa | Severely truncated |
| P. xalapensis | PsiM | 25 aa | ~310 aa | Severely truncated |
| P. yungensis | PsiM | 25 aa | ~310 aa | Severely truncated |
| P. yungensis | PsiH | 291 aa | ~510 aa | Truncated |

### Non-Psilocybin Producing Species (NO_HIT for all genes)

These 20 species lack the psilocybin biosynthesis gene cluster entirely:

- P. argentina
- P. apelliculosa
- P. chionophila
- P. clavata
- P. fuliginosa
- P. laticystis
- P. lazoi
- P. magica
- P. montana
- P. phyllogena
- P. polytrichophila
- P. rhomboidospora
- P. sabulosa
- P. strictipes
- P. subcoprophila
- P. subpsilocybioides
- P. subviscida
- P. tuberosa
- P. washingtonensis

## Directory Structure

```
05_Felix_code/
├── second_try.ipynb          # Main extraction notebook
├── README.md                 # This file
│
├── reference_panels/         # P. cubensis reference proteins
│   ├── Proteins_PsiD/
│   │   └── PsiD_refs.faa
│   ├── Proteins_PsiH/
│   │   └── PsiH_refs.faa
│   ├── Proteins_PsiK/
│   │   └── PsiK_refs.faa
│   └── Proteins_PsiM/
│       └── PsiM_refs.faa
│
├── cds_extraction_outputs/   # Per-species extraction results
│   ├── Psilocybe_*/
│   │   ├── *_PsiD.cds.fa     # Nucleotide CDS
│   │   ├── *_PsiD.prot.fa    # Protein sequence
│   │   ├── *_PsiK.cds.fa
│   │   ├── *_PsiK.prot.fa
│   │   ├── *_PsiM.cds.fa
│   │   ├── *_PsiM.prot.fa
│   │   ├── *_PsiH.cds.fa
│   │   ├── *_PsiH.prot.fa
│   │   └── *_summary.json    # Extraction metadata
│   ├── extraction_log.json   # Full run log
│   ├── problems_report.tsv   # All issues flagged
│   └── species_problem_summary.tsv
│
├── alphafold_input/          # Combined FASTA files for AlphaFold
│   ├── PsiD_sequences.faa    # 51 sequences
│   ├── PsiK_sequences.faa    # 51 sequences
│   ├── PsiM_sequences.faa    # 45 sequences
│   └── PsiH_sequences.faa    # 51 sequences
│
├── mafft_alignments/         # Multiple sequence alignments
│   ├── PsiD_mafft_aligned.faa
│   ├── PsiH_mafft_aligned.faa
│   ├── PsiK_mafft_aligned.faa
│   ├── PsiM_mafft_aligned.faa
│   ├── alignment_metadata.tsv
│   ├── alignment_viewer.html  # Interactive Plotly viewer
│   └── nucleotide_alignments/
│       ├── PsiD_cds_aligned.fna
│       ├── PsiH_cds_aligned.fna
│       ├── PsiK_cds_aligned.fna
│       └── PsiM_cds_aligned.fna
│
└── protein_structure_prediction/
    └── 01_ColabFold_Psilocybin.ipynb  # ColabFold structure prediction
```

## Protein Structure Prediction

The `protein_structure_prediction/01_ColabFold_Psilocybin.ipynb` notebook predicts 3D protein structures for all extracted psilocybin enzymes using **ColabFold** (AlphaFold2 with MMseqs2). It is designed to run on **Google Colab** with a GPU runtime (A100 or T4).

### What it does

1. **Environment setup** - Installs ColabFold, JAX with CUDA support, and HHsuite for template search
2. **Input** - Reads the combined FASTA files from `alphafold_input/` (PsiD, PsiK, PsiM, PsiH) plus ancestral reconstructions (`ASR1_PsiD`)
3. **Structure prediction** - Runs ColabFold batch predictions with configurable parameters:
   - 5 models per sequence (AlphaFold2 model ensemble)
   - 6 recycles for improved accuracy
   - Template-based modeling enabled
4. **Quality assessment** - Extracts pLDDT and pTM confidence scores from prediction outputs
5. **Results** - Saves PDB structures and quality summary CSV to Google Drive

### Current results

| Enzyme | Sequences Predicted | PDB Files |
|--------|---------------------|-----------|
| PsiD | 53 | 265 |
| PsiK | 24 | 120 |
| PsiM | 0 | 0 |
| PsiH | 0 | 0 |
| ASR1_PsiD | 66 | 330 |

PsiM and PsiH predictions are pending.

## Output Files

### AlphaFold Input (`alphafold_input/`)

Combined FASTA files ready for structure prediction:

| File | Sequences | Description |
|------|-----------|-------------|
| `PsiD_sequences.faa` | 51 | Tryptophan decarboxylase |
| `PsiK_sequences.faa` | 51 | Kinase |
| `PsiM_sequences.faa` | 45 | Methyltransferase |
| `PsiH_sequences.faa` | 51 | P450 monooxygenase |

**All files are clean** - no internal stop codons (`*`).

### MAFFT Alignments (`mafft_alignments/`)

Multiple sequence alignments for phylogenetic analysis:

- **Protein alignments**: `*_mafft_aligned.faa`
- **Nucleotide alignments**: `nucleotide_alignments/*_cds_aligned.fna`
- **Interactive viewer**: `alignment_viewer.html` (open in browser)

### Per-Species Results (`cds_extraction_outputs/`)

Each species folder contains:
- `.cds.fa` - Nucleotide coding sequence
- `.prot.fa` - Translated protein
- `_summary.json` - Extraction metadata (score, strand, status)

## Notebook Cells

| Cell | Description |
|------|-------------|
| **Cell 1** | Setup & config - loads references from `03_Sequences_Paper/References/` |
| **Cell 2** | Core extraction functions using exonerate `--ryo` |
| **Cell 3** | Batch parallel processing (71 species, ~6 min) |
| **Cell 3b** | Re-run failed species (serial mode) |
| **Cell 3c** | Export combined FASTA for AlphaFold |
| **Cell 4** | MAFFT protein alignments |
| **Cell 5-6** | Plotly alignment visualization |
| **Cell 7** | Nucleotide alignments |

## Requirements

```bash
# System tools
apt-get install exonerate  # or brew install exonerate

# Python packages
pip install biopython pandas tqdm plotly ipywidgets
```

## Usage

1. **Run Cell 1** to set up paths and build reference panels
2. **Run Cell 2** to define extraction functions (test with P. baeocystis)
3. **Run Cell 3** to batch process all 71 species
4. **Run Cell 3b** if any species failed with errors
5. **Run Cell 3c** to export combined FASTA files for AlphaFold
6. **Run Cells 4-7** for alignments and visualization

## References

- Bradshaw et al. (2024) "Phylogenomics of the psychoactive mushroom genus Psilocybe"
- Reference proteins: *Psilocybe cubensis* (GenBank: ASU62239.1 PsiD)

## Citation

If you use this pipeline, please cite:
- The Bradshaw et al. (2024) paper for the genome assemblies
- Exonerate: Slater & Birney (2005) BMC Bioinformatics

## License

MIT License
