import argparse, os, subprocess
from Bio import SeqIO
import matplotlib.pyplot as plt

def blast_genes(faa, db, threads):
    arg = f'diamond blastp --threads {threads} --max-target-seqs 1 --db {db} --query {faa} --outfmt 6 qseqid sseqid slen pident'
    results = subprocess.check_output(arg, shell=True).decode()

    hits = [line.split("\t") for line in results.splitlines()]
    hits = {hit[0]:hit[1:] for hit in hits}

    return(hits)

def plot_hist(outfix, mode, values):
    plt.style.use('seaborn')
    plt.figure()
    plt.ylabel("Count")
    if mode=="len":
        plt.hist(values, 21)
        plt.title("Gene Lengths")
        plt.xlabel("Gene length (bp)")
        plt.savefig(f'{outfix}_gene_lengths.png')
    elif mode=="ratio":
        plt.hist(values, bins=[-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1])
        plt.title("Gene Length Ratios")
        plt.xlabel("Gene length ratio")
        plt.savefig(f'{outfix}_gene_ratios.png')

def __main__():
    parser = argparse.ArgumentParser(description='Call genes with prodigal, calculate gene length distribution and optionally compare to the nearest BLAST hits.')
    parser.add_argument('contigs', metavar='contigs_file', help='File containing contigs in fasta format.')
    parser.add_argument('--blast', action='store_true', help='Compare gene lengths to nearest BLAST hits in db.')
    parser.add_argument('--db', metavar='ref_db', default='/nfs/cds/Databases/UNIPROT/uniprot_trembl_bacteria.dmnd', help='Diamond reference database for gene alignment.')
    parser.add_argument('--cutoff', metavar='id_cutoff', default=97.0, help='Identity cutoff for a valid hit. Default: 97.0.')
    parser.add_argument('--plot', action='store_true', help='Plot histograms of gene length and ratio (if --blast).')
    parser.add_argument('-t,', '--threads', metavar='threads', default=16, help='Number of compute threads. Default: 16.')

    args = parser.parse_args()

    # Determine file names
    dir = os.path.split(args.contigs)[0]
    outdir = os.path.join(dir, "genlen_results")

    prefix = os.path.splitext(os.path.split(args.contigs)[1])[0]
    outfix = os.path.join(outdir, prefix)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Call genes with prodigal
    prodigal_output = subprocess.check_output(f'prodigal -c -i {args.contigs} -o {outfix}.gff -a {outfix}.faa -d {outfix}.fna -f gff',shell=True).decode()

    # Read in protein sequences
    proteins = SeqIO.parse(f'{outfix}.faa', 'fasta')

    # Determine gene lengths
    gene_lengths = {protein.id:len(protein)-1 for protein in proteins}

    if args.blast:
        # Blast protein sequences against reference database
        hits = blast_genes(f'{outfix}.faa', args.db, args.threads)

        # Filter hits
        good_hits = {k:v for k,v in hits.items() if float(v[2]) >= args.cutoff}

        # Determine hit lengths
        hit_lengths = {k:int(v[1]) for k,v in good_hits.items()}

        # Calculate ratios
        length_ratios = {k:gene_lengths[k]/hit_lengths[k] for k in hit_lengths.keys()}

    # Output
    with open(f'{outfix}_results.txt', 'w') as fo:
        for k in gene_lengths.keys():
            try:
                fo.write(f'{k}\t{gene_lengths[k]}\t{hit_lengths[k]}\t{length_ratios[k]}\n')
            except:
                fo.write(f'{k}\t{gene_lengths[k]}\tNA\tNA\n')

    # Graphs
    if args.plot:
        plot_hist(outfix, "len", gene_lengths.values())
        if args.blast:
            plot_hist(outfix, "ratio", length_ratios.values())

if __name__ == "__main__":
    __main__()

