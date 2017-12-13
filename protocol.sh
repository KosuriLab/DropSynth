#!/bin/bash
#protocol.sh

echo "Hi, $USER. This is the main DropSynth oligo generation script."
echo "Please read the README.md before running"

echo "Screening assembly and amplification primers."
#screen primers
python primer_screen.py

echo "Downloading DHFR records from NCBI Entrez. (slow)"
#get XML files for each homolog (can take long if download required)
python entrez_get_XML.py

echo "Sort records by sequence length and finding antibiotic resistant homologs."
#sort sequences by length into those needed 4 oligo and 5 oligo for assembly
#find DHFR sequences which are resistant (put into a separate library later)
python sort_by_length_and_find_resistant_seq.py

echo "The resulting file:"
echo "db/DHFR_IP_targetlist_parsed_4_5_res_oligo.csv"
echo "should be manually edited to add or remove sequences as desired."
echo "The edited file should be saved to:"
echo "db/DHFR_IP_targetlist_parsed_4_5_res_oligo_toFASTA.csv"

read -n1 -r -p "Press any key to continue..." key

#(optional) generate identity distribution (very slow)
#python identity_dist.py

#(optional) generate identity to E. coli and Human sequences (very slow)
#python identity_Ecoli.py

echo "Generating library sequence input files."
#take sequence files and split them into files for 4 oligo and 5 oligo
python generate_libraries.py

echo "Generating genes and splitting into oligo payload."
#splits genes to use a payload in DropSynth oligos (without amp primers, microbead barcodes, or nt.BspQI nicking sites)
python split_genes_for_oligos.py

echo "Generating CSV file with all library info."
#make a file with all info needed downstream
python make_master_CSV.py

#(optional) screen primers for off-target hybridization
#./BLAT_primers.sh
#python BLAT_hits_parse.py

echo "Finalizing DropSyth oligos. Adding amplification primers, microbead barcodes, and nicking sites."
#finalize DropSynth oligos by adding amp primers, microbead barcodes, and nt.BspQI nicking sites
python oligobuffergen_DHFR.py

echo "Final oligo check of design rules."
#check as many design rules as possible
#do not resuse code in this script to avoid propagating errors
python final_oligo_check.py

echo "Making CSV file with protein info."
#make a csv file with all protein seqs (also done in make_master_CSV.py)
python FASTA_to_csv_proteins.py

echo "Making CSV file with oligos."
#make a csv file with all oligo seqs to order microarray
python FASTA_to_csv_for_Chip.py

echo "Checking oligos for homopolymers."
#check homopolymers
python homopolymer_stats.py

echo "Choosing which oligo strand to synthesize (fewest Adenines)."
#choose the strand with fewest Adenines for synthesis (depurination)
python choose_lowest_A.py chip_out/DHFR_chip_3_raw.csv chip_out/DHFR_chip_3_raw_RC_lowest_A.csv

echo "Done."