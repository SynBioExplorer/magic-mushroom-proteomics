# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics pipeline for extracting psilocybin biosynthesis genes (PsiD, PsiK, PsiM, PsiH) from fungal genome assemblies. The project uses Exonerate for protein-to-genome alignment to identify and extract coding sequences (CDS) from scaffold FASTAs, based on reference protein panels derived from the Bradshaw et al. (2024) paper on Psilocybe phylogenomics.

## Key Architecture

### Data Flow Pipeline

1. **Reference Panel Creation** (Cell 1):
   - Builds four enzyme-specific reference panels from DNA sequences in `03_Sequences_Paper/DNA_All/`
   - Falls back to protein predictions from `03_Sequences_Paper/Protein_predictions/` if available
   - Outputs: `reference_panels/Proteins_{PsiD,PsiK,PsiM,PsiH}/{enzyme}_refs.faa`

2. **CDS Extraction** (Cell 2):
   - For each genome scaffold in `03_Sequences_Paper/Assembly_scaffolds/`:
     - Runs Exonerate (protein2genome model) against each enzyme reference panel
     - Parses GFF output to identify best hit based on scoring criteria
     - Extracts CDS, handles reverse-complement for minus-strand hits
     - Translates CDS to protein and validates reading frame
   - Outputs per species: `cds_extraction_outputs/{species}/{species}_{enzyme}.{cds,prot}.fa`

3. **Batch Processing** (Cell 3):
   - Parallel processing of all species using ProcessPoolExecutor
   - Generates per-species JSON summaries and global extraction log

### Critical Dependencies

- **Exonerate**: Must be installed and on PATH (`brew install exonerate` on macOS)
- **BioPython**: For FASTA parsing, sequence manipulation, translation
- **Data Directory**: The notebook auto-discovers `03_Sequences_Paper/` relative to the notebook location
- **tqdm** (optional): For progress bars (`pip install tqdm`)

### Optimization Mode (Default)

**The notebook is optimized for speed** using:
- **Single reference**: P. cubensis only (not multi-species panel)
- **BESTN = 1**: Best hit only (vs. top 3)
- **Result**: ~3-4x faster than comprehensive mode

This is ideal for the Bradshaw dataset where all 71 species are closely related Psilocybe.

### Scoring Logic

The `score_candidate()` function prioritizes hits by:
1. No internal stop codons (critical)
2. In-frame translation (length divisible by 3)
3. Within expected amino acid length bands (enzyme-specific soft constraints)
4. Higher Exonerate alignment score
5. Longer CDS length

Length bands (amino acids):
- PsiD: 350-600 (TDC, ~480 typical)
- PsiK: 250-420 (kinase)
- PsiM: 180-320 (methyltransferase)
- PsiH: 400-650 (P450, ~510 typical)

## Running the Pipeline

### Prerequisites
```bash
# Ensure exonerate is installed
exonerate --version

# Verify Python environment has biopython
python -c "import Bio; print(Bio.__version__)"
```

### Execution Modes

**Single Species Test:**
```python
# Cell 2: Modify the glob pattern to target specific species
one_species = sorted(SCAFF_DIR.glob("Psilocybe_baeocystis*.fa*"))[0]
_ = run_species(one_species, OUT_DIR)
```

**Batch Processing:**
```python
# Cell 3: Run all species in parallel with progress tracking
# Adjust max_workers based on CPU cores (default: min(12, cpu_count))
max_workers = min(12, os.cpu_count() or 8)
# Progress bars (via tqdm) show:
# - Overall species progress (e.g., "Species 15/71")
# - Per-species enzyme progress (e.g., "PsiD, PsiK, PsiM, PsiH")
```

### Tunable Parameters

In Cell 1:
- `REFERENCE_SPECIES = "Psilocybe_cubensis"`: Which species to use as reference
- Change this to use a different reference (must exist in DNA_All)

In Cell 2:
- `BESTN = 1`: Number of top hits (1 = fastest, 3 = more thorough)
- `PERCENT = 30`: Minimum percent identity for Exonerate alignment
- `AA_BANDS`: Expected amino acid length ranges (can be set to `None` to disable)

## Output Structure

```
cds_extraction_outputs/
├── {species_name}/
│   ├── {species}_{enzyme}.cds.fa      # nucleotide CDS
│   ├── {species}_{enzyme}.prot.fa     # translated protein
│   └── {species}_summary.json         # extraction metadata
├── extraction_log.json                # global batch log with all results
├── problems_report.tsv                # detailed table of all problematic genes
└── species_problem_summary.tsv        # per-species problem counts
```

### Status Codes in Summaries
- `OK`: Clean extraction (in-frame, no stop codons)
- `CHECK_ORF`: Out-of-frame or has internal stops (manual review needed)
- `NO_HIT`: No alignment met minimum thresholds
- `ERROR`: Exception during processing (see error field)

### QC Reports

The batch processing generates two TSV files for quality control:

1. **problems_report.tsv**: Lists all genes with issues (CHECK_ORF, NO_HIT, ERROR)
   - Columns: species, enzyme, status, details
   - Use for manual curation and debugging

2. **species_problem_summary.tsv**: Per-species summary of problem counts
   - Columns: species, total_problems, CHECK_ORF, NO_HIT, ERROR
   - Helps identify species requiring attention

## Common Development Tasks

### Adjusting Alignment Sensitivity
To increase hit rate for divergent sequences, lower the identity threshold:
```python
PERCENT = 25  # from 30
```

### Debugging Failed Extractions
Check individual species logs:
```python
import json
with open("cds_extraction_outputs/{species}/{species}_summary.json") as f:
    print(json.dumps(json.load(f), indent=2))
```

### Re-running Single Enzyme
Modify `REFS` dict to include only target enzyme before calling `run_species()`.

## Data Dependencies

The pipeline expects `03_Sequences_Paper/` at the parent directory level, containing:
- `Assembly_scaffolds/`: Genome FASTAs (`.fa`, `.fasta`, `.fa.gz`)
- `DNA_All/{PsiD,PsiK,PsiM,PsiH}/`: Reference DNA sequences with strand info
- `Protein_predictions/`: Optional `.faa` files with enzyme annotations

## Known Limitations

- Exonerate alignments with low score may miss true positives (tune `PERCENT` threshold)
- CDS extraction assumes standard genetic code (table 1)
- Reference panel quality directly affects hit sensitivity
- Parallelization limited by Exonerate CPU usage (I/O-bound on many cores)

## Performance Notes

### Exonerate Runtime

Exonerate runtime is primarily determined by:
1. **Number of scaffolds** in the target genome (not total genome size)
2. **BESTN parameter** (checking top N hits)
3. **PERCENT identity threshold** (lower = more candidates to evaluate)

Example runtimes per enzyme on different genome assemblies (OPTIMIZED MODE):
- Well-assembled genome (~100 scaffolds): 2-5 seconds per enzyme
- Fragmented genome (~100,000 scaffolds): 60-120 seconds per enzyme

**P. baeocystis example**: With 111,256 scaffolds:
- **Optimized mode** (bestn=1, single ref): ~60-90s per enzyme (~4-6 min total)
- Old mode (bestn=3, 60 refs): ~180-240s per enzyme (~12-16 min total)
- **Speedup: 3-4x faster!**

### Switching Between Modes

**Current: Optimized Mode** (default, recommended for Psilocybe species)
- Uses P. cubensis single reference
- BESTN = 1
- ~3-4x faster

**To Enable Comprehensive Mode** (for divergent species):

In Cell 1, change:
```python
USE_SINGLE_REFERENCE = False  # Use multi-species panel
```

In Cell 2, change:
```python
BESTN = 3  # Check top 3 hits
```

Trade-off: Slower but better for distantly related species or detecting duplications.

The verbose logging shows per-enzyme timing to help identify problematic genomes.
