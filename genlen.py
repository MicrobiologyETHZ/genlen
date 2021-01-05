import argparse, os, subprocess, sys
from Bio import Seq, SeqIO
from io import StringIO

def blast_genes(faa, db, threads, id):
    arg = "diamond blastp --threads {} --max-target-seqs 1 --db {} --query {} --outfmt 6 qseqid sseqid qlen slen pident length".format(threads, db, faa)
    results = subprocess.check_output(arg, shell=True).decode()

    # Parse results
    hits = results.splitlines()
    sys.stderr.write("{} genes with blast hits\n".format(len(hits)))
    hits = [line.split('\t') for line in hits]

    # Filter for valid hits
    valid_hits = [hit for hit in hits if float(hit[4]) >= float(id)]
    sys.stderr.write("{} genes with blast hits > {}% identity\n".format(len(valid_hits), id))

    return(valid_hits)

if True:
#def __main__():
    parser = argparse.ArgumentParser(description='Call genes with prodigal, calculate gene length distribution and compare to the nearest BLAST hits')
    parser.add_argument('contigs', metavar='contigs_file', help='File containing contigs in fasta format')
    parser.add_argument('--db', metavar='ref_db', default='/nfs/cds/Databases/DIAMOND/nr', help='Reference database for gene alignment')
    parser.add_argument('--id', metavar='id_cutoff', default=80.0, help='Identity threshold')
    parser.add_argument('-t,', '--threads', metavar='threads', default=16, help='Number of compute threads')

    args = parser.parse_args()

    # Determine file name prefix
    prefix = os.path.splitext(args.contigs)[0]

    # Call genes with prodigal
    prodigal_output = subprocess.check_output('prodigal -c -i {} -o {}.gff -a {}.faa -d {}.fna -f gff'.format(args.contigs, prefix, prefix, prefix),shell=True).decode()

    # Blast protein sequences against reference database
    valid_hits = blast_genes("{}.faa".format(prefix), args.db, args.threads, args.id)

    # Calculate length ratios
    len_ratios = [float(hit[2])/float(hit[3]) for hit in valid_hits]
    len_ratios = sorted(len_ratios)

    for x in len_ratios:
        sys.stdout.write("{}\n".format(str(x)))

#if __name__ == "__main__":
#    __main__()

