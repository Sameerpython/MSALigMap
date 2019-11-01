MSALigMap (Multiple Sequence Alignment Ligand Mapping) is a tool for mapping active site amino acid residues that bind selected ligands on to target protein sequences of interest.

MSALigMap requires Python 2.7. The running dependencies are:
1. Python 2.7 or later
2. Installation of ClustalO

Example:

A basic example to execute MSALigMap is:

Usage: MSALigMap <input.fasta> <ligand.txt>

Input: Fasta file containing sequences that includes sequences of PDB structures and non-PDB sequences
Ligand: is the file that contains PDB code and bound ligand code from PDB database.

The format of PDB ids and ligand names in ligand.txt should be as follows:

3wxb:A	NDP
3o26:A	NDP

A pdb code is followed by ':' and the ligand name (code as provided in PDB database) is separated by a tab.

In the sequence file, the sequence name for the PDB structures should also be similar as provided in ligand.txt file. For eg:

>3WXB:A
MGSSHHHHHHSSGLVPRGSHMGELRVRSVLVTGANRGIGLGFVQHLLALSNPPEWVFATCRDPKGQRAQELQKLASKHPN
LVIVPLEVTDPASIKAAAASVGERLKGSGLNLLINNAGIARANTIDNETLKDMSEVYTTNTIAPLLLSQAFLPMLKKAAQ
ENPGSGLSCSKAAIINISSTAGSIQDLYLWQYGQALSYRCSKAALNMLTRCQSMGYREHGIFCVALHPGWVKTDMGGTLE
DKSRVTVDESVGGMLKVLSNLSEKDSGAFLNWEGKVMAW

Sample file can be found in example folder.

