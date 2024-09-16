# lRM001_pegRNA_screen
scripts related to the analysis of lRM001

note: only ~3700 out of the approximately 6000 pegRNA designs in the TWIST library are related to the variant tiling screen, see _____ for details including the relevant barcodes

**to run demux_pegRNA_screen_2024.py:**

python3 demux_pegRNA_screen_2024.py ${library}-read-1.fastq.gz ${library}-read-3.fastq.gz

e.g. python3 demux_pegRNA_screen_2024.py 3018__0__HTS_RAM_062-read-1.fastq.gz 3018__0__HTS_RAM_062-read-3.fastq.gz


**to run analyze_pegRNA_screen_2024_v2.py:**

python3 analyze_pegRNA_screen_2024_v2.py pegRNA_index library_name

e.g. python3 analyze_pegRNA_screen_2024_v2.py 0 3018__0__HTS_RAM_062


