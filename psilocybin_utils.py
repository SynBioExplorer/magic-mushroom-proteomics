"""
Utility functions for Psilocybin BGC Ancestral Sequence Analysis

This module provides helper functions for:
- Running exonerate for gene extraction
- Quality control of sequences
- Running bioinformatics tools (MAFFT, IQ-TREE, etc.)
- Sequence manipulation and analysis

Based on Bradshaw et al. (2024) workflow
"""

import subprocess
import os
from pathlib import Path
from typing import Union, List, Dict, Optional
from Bio import SeqIO, Seq, AlignIO
from Bio.SeqRecord import SeqRecord


class ExonerateRunner:
    """Wrapper for running exonerate protein2genome alignments."""

    def __init__(self, exonerate_path: str = "exonerate"):
        """Initialize ExonerateRunner.

        Parameters
        ----------
        exonerate_path : str
            Path to exonerate executable (default: assumes in PATH)
        """
        self.exonerate_path = exonerate_path

    def run_protein2genome(
        self,
        query_protein: Union[str, Path],
        target_genome: Union[str, Path],
        output_cds: Union[str, Path],
        output_gff: Optional[Union[str, Path]] = None,
        bestn: int = 1,
        verbose: bool = True
    ) -> bool:
        """
        Run exonerate protein2genome to extract CDS from genome.

        Parameters
        ----------
        query_protein : str or Path
            Reference protein FASTA file
        target_genome : str or Path
            Target genome scaffold FASTA file
        output_cds : str or Path
            Output CDS FASTA file
        output_gff : str or Path, optional
            Output GFF annotation file
        bestn : int
            Number of best hits to report (default: 1)
        verbose : bool
            Print command and output

        Returns
        -------
        bool
            True if successful, False otherwise
        """
        # Extract CDS
        cmd_cds = [
            self.exonerate_path,
            "--model", "protein2genome",
            str(query_protein),
            str(target_genome),
            "--bestn", str(bestn),
            "--showalignment", "no",
            "--showvulgar", "no",
            "--verbose", "0",
            "--ryo", ">%ti:%tab-%tae(%tS)\\n%tcs\\n"
        ]

        if verbose:
            print(f"Running exonerate CDS extraction: {' '.join(cmd_cds)}")

        with open(output_cds, 'w') as outfile:
            result = subprocess.run(cmd_cds, stdout=outfile, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(f"✗ Exonerate CDS extraction failed: {result.stderr}")
            return False

        if verbose:
            print(f"✓ CDS extracted: {output_cds}")

        # Extract GFF if requested
        if output_gff:
            cmd_gff = [
                self.exonerate_path,
                "--model", "protein2genome",
                str(query_protein),
                str(target_genome),
                "--bestn", str(bestn),
                "--showalignment", "no",
                "--showvulgar", "no",
                "--verbose", "0",
                "--showtargetgff", "yes"
            ]

            with open(output_gff, 'w') as outfile:
                result = subprocess.run(cmd_gff, stdout=outfile, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                print(f"✗ Exonerate GFF extraction failed: {result.stderr}")
                return False

            if verbose:
                print(f"✓ GFF extracted: {output_gff}")

        return True


def translate_sequence(
    input_cds: Union[str, Path],
    output_protein: Union[str, Path],
    frame: int = 1,
    transeq_path: str = "transeq"
) -> bool:
    """
    Translate CDS to protein using EMBOSS transeq.

    Parameters
    ----------
    input_cds : str or Path
        Input CDS FASTA file
    output_protein : str or Path
        Output protein FASTA file
    frame : int
        Reading frame (default: 1)
    transeq_path : str
        Path to transeq executable

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    cmd = [
        transeq_path,
        "-sequence", str(input_cds),
        "-outseq", str(output_protein),
        "-frame", str(frame)
    ]

    print(f"Translating CDS: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"✓ Translation complete: {output_protein}")
        return True
    else:
        print(f"✗ Translation failed: {result.stderr}")
        return False


def quality_control(
    protein_file: Union[str, Path],
    check_stops: bool = True,
    check_length: bool = True,
    min_length: int = 100,
    max_stops: int = 0
) -> Dict[str, any]:
    """
    Perform quality control checks on protein sequences.

    Parameters
    ----------
    protein_file : str or Path
        Protein FASTA file
    check_stops : bool
        Check for internal stop codons
    check_length : bool
        Check sequence length
    min_length : int
        Minimum acceptable protein length (aa)
    max_stops : int
        Maximum acceptable internal stop codons

    Returns
    -------
    dict
        QC results with keys: 'passes_qc', 'stop_codons', 'length', 'warnings'
    """
    results = {
        'passes_qc': True,
        'stop_codons': 0,
        'length': 0,
        'warnings': []
    }

    try:
        record = SeqIO.read(protein_file, "fasta")
        sequence = str(record.seq)

        # Check length
        results['length'] = len(sequence)
        if check_length and results['length'] < min_length:
            results['passes_qc'] = False
            results['warnings'].append(f"Sequence too short: {results['length']} aa < {min_length} aa")

        # Check for stop codons (excluding terminal)
        if check_stops:
            internal_seq = sequence[:-1] if sequence.endswith('*') else sequence
            results['stop_codons'] = internal_seq.count('*')

            if results['stop_codons'] > max_stops:
                results['passes_qc'] = False
                results['warnings'].append(
                    f"Too many internal stop codons: {results['stop_codons']} > {max_stops}"
                )

        if results['passes_qc']:
            print(f"✅ QC passed: {protein_file.name if hasattr(protein_file, 'name') else protein_file}")
        else:
            print(f"⚠️  QC warnings for {protein_file.name if hasattr(protein_file, 'name') else protein_file}:")
            for warning in results['warnings']:
                print(f"   - {warning}")

    except Exception as e:
        results['passes_qc'] = False
        results['warnings'].append(f"Error reading file: {str(e)}")
        print(f"✗ QC failed: {str(e)}")

    return results


def run_mafft(
    input_fasta: Union[str, Path],
    output_fasta: Union[str, Path],
    algorithm: str = "auto",
    threads: int = -1,
    mafft_path: str = "mafft"
) -> bool:
    """
    Run MAFFT multiple sequence alignment.

    Parameters
    ----------
    input_fasta : str or Path
        Input unaligned FASTA file
    output_fasta : str or Path
        Output aligned FASTA file
    algorithm : str
        MAFFT algorithm: 'auto', 'linsi', 'ginsi', 'einsi' (default: auto)
    threads : int
        Number of threads (-1 for auto, default: -1)
    mafft_path : str
        Path to MAFFT executable

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    # Build command
    cmd = [mafft_path]

    if algorithm != "auto":
        cmd.append(f"--{algorithm}")
    else:
        cmd.append("--auto")

    cmd.extend([
        "--reorder",
        "--thread", str(threads),
        str(input_fasta)
    ])

    print(f"Running MAFFT: {' '.join(cmd)}")

    with open(output_fasta, 'w') as outfile:
        result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)

    if result.returncode == 0:
        print(f"✓ Alignment complete: {output_fasta}")
        return True
    else:
        print(f"✗ MAFFT failed: {result.stderr}")
        return False


def run_iqtree(
    alignment_file: Union[str, Path],
    output_prefix: Union[str, Path],
    model: str = "MFP",
    bootstrap: int = 1000,
    threads: str = "AUTO",
    iqtree_path: str = "iqtree"
) -> bool:
    """
    Run IQ-TREE for phylogenetic tree construction.

    Parameters
    ----------
    alignment_file : str or Path
        Input alignment file
    output_prefix : str or Path
        Prefix for output files
    model : str
        Substitution model (default: MFP for ModelFinder Plus)
    bootstrap : int
        Number of ultrafast bootstrap replicates (default: 1000)
    threads : str
        Number of threads (default: AUTO)
    iqtree_path : str
        Path to IQ-TREE executable

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    cmd = [
        iqtree_path,
        "-s", str(alignment_file),
        "-pre", str(output_prefix),
        "-m", model,
        "-bb", str(bootstrap),
        "-nt", threads
    ]

    print(f"Running IQ-TREE: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"✓ IQ-TREE analysis complete")
        print(f"  Tree file: {output_prefix}.treefile")
        return True
    else:
        print(f"✗ IQ-TREE failed: {result.stderr}")
        return False


def combine_fasta_files(
    input_files: List[Union[str, Path]],
    output_file: Union[str, Path],
    add_species_prefix: bool = True
) -> int:
    """
    Combine multiple FASTA files into one.

    Parameters
    ----------
    input_files : list of str or Path
        List of input FASTA files
    output_file : str or Path
        Output combined FASTA file
    add_species_prefix : bool
        Add species name as prefix to sequence IDs (default: True)

    Returns
    -------
    int
        Number of sequences combined
    """
    combined_records = []

    for fasta_file in input_files:
        fasta_path = Path(fasta_file)
        species_name = fasta_path.stem.replace('_PsiD', '').replace('.', '_')

        for record in SeqIO.parse(fasta_file, "fasta"):
            if add_species_prefix:
                record.id = f"{species_name}|{record.id}"
                record.description = f"{species_name} {record.description}"
            combined_records.append(record)

    SeqIO.write(combined_records, output_file, "fasta")
    print(f"✓ Combined {len(combined_records)} sequences into {output_file}")

    return len(combined_records)


# BGC gene information
BGC_GENES = {
    'PsiD': {
        'name': 'Tryptophan decarboxylase',
        'function': 'Converts tryptophan to tryptamine',
        'color': 'red'
    },
    'PsiK': {
        'name': 'Kinase',
        'function': 'Phosphorylates intermediate',
        'color': 'green'
    },
    'PsiM': {
        'name': 'Methyltransferase',
        'function': 'Methylates tryptamine',
        'color': 'gold'
    },
    'PsiH': {
        'name': 'P450 monooxygenase',
        'function': 'Hydroxylates at C-4 position (critical for activity)',
        'color': 'purple'
    }
}

# BGC organization by clade
BGC_ORGANIZATION = {
    'Clade_I': ['PsiD', 'PsiK', 'PsiH', 'PsiM'],
    'Clade_II': ['PsiD', 'PsiM', 'PsiH', 'PsiK']
}


if __name__ == "__main__":
    print("Psilocybin BGC Analysis Utilities")
    print("=" * 50)
    print("\nBGC Genes:")
    for gene, info in BGC_GENES.items():
        print(f"  {gene}: {info['name']}")
        print(f"      {info['function']}")

    print("\nBGC Organization:")
    for clade, order in BGC_ORGANIZATION.items():
        print(f"  {clade}: {' > '.join(order)}")
