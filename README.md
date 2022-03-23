# genlen

Call genes with prodigal, calculate gene length distribution and optionally compare to the nearest BLAST hits.

positional arguments:  
  contigs_file          File containing contigs in fasta format.  

optional arguments:  
  -h, --help                        show this help message and exit  
  --blast                           Compare gene lengths to nearest BLAST hits in db.  
  --db ref_db                       Diamond reference database for gene alignment.  
  --cutoff id_cutoff                Identity cutoff for a valid hit. Default: 97.0.  
  --plot                            Plot histograms of gene length and ratio (if --blast).  
  -t, threads, --threads threads    Number of compute threads. Default: 16.  

