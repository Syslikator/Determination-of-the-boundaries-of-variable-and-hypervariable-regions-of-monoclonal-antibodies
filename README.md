# Determination-of-the-boundaries-of-variable-and-hypervariable-regions-of-monoclonal-antibodies
Vabs.py - Protocol for searching for variable regions of mouse antibodies. Translates and identifies variable regions of kappa, lambda, heavy chains by searching for common amino acid motifs.<br>
<br>
CDRs.py - Protocol for searching for hypervariable regions in mouse antibodies. Boundaries are determined by searching for specific sequences flanking hypervariable regions.

Algorithms (scripts) VAbs.py and CDRs.py are written using python. <br>
The script structure accepts an input file in FASTA format containing the cDNA sequence of the heavy or light chain of the antibody.
This cDNA is a nucleotide sequence containing information about the primary structure of the antibody and includes all C, V, (D) and J or only V (D) and J segments of the somatic recombination of heavy or light chain immunoglobulin genes. <br><br>
The main part of the algorithm includes information about motifs and combinations of amino acids defining conserved regions of the variable domain (VAbs.py) and amino acid combinations flanking hypervariable regions (CDRs.py). Each strand type has its own argument name for object methods (cDNA sequences). <br>
Algorithms are run through the terminal. The program path is specified (VAbs.py or CDRs.py). Then, using the assignment function (<), the name of the file containing the cDNA sequence is entered. The results indicate the chain type (kappa, lambda, heavy), the amino acid motifs of the conserved region of the variable domain, and the translated sequence of the studied cDNA.

