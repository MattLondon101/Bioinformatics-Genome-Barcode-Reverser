# Genome-Barcoder-Reverser
Command line application to get reverse barcodes and their frequencies in genome sequences.

1. Run barcode_reverse_pair.py in Linux Terminal

2. You will be given three user input prompts:   
"Enter path to .fastq file:"  
"Enter start position of sequence:"  
"Enter stop position of sequence:"  

3. The terminal will output the time in seconds (s) it took for the program to execute and the memory usage of the program in megabytes (MB):  
"Program execution time = * s"  
"Memory usage of program = * MB"

4. Two CSV files will be written to the wkdir:  
* paired_barcodes_revcomps.csv : Barcode and Reverse Barcode pairs and their frequency counts.  
* unpaired_barcodes.csv. : Unpaired barcodes and their frequency counts.
