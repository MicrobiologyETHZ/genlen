import argparse, os, subprocess
from Bio import SeqIO
import matplotlib.pyplot as plt

def align_genes(faa, db, threads):
    arg = f'diamond blastp --threads {threads} --db {db} --query {faa} --unal 1 --outfmt 6 qseqid sseqid qlen slen pident'
    results = subprocess.check_output(arg, shell=True).decode()

    hits = [line.split("\t") for line in results.splitlines()]
    genes = set(hit[0] for hit in hits)
    hits = [[hit for hit in hits if hit[0]==y] for y in genes]
    hits = {hit[0][0]:hit[0] for hit in hits}

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

if True:
#ef __main__():
    parser = argparse.ArgumentParser(description='Call genes with prodigal, calculate gene length distribution and optionally compare to the nearest BLAST hits.')
    parser.add_argument('contigs', metavar='contigs_file', help='File containing contigs in fasta format.')
    parser.add_argument('--align', action='store_true', help='Find best hits in db with DIAMOND.')
    parser.add_argument('--db', metavar='ref_db', default='/nfs/cds/Databases/UNIPROT/uniprot_trembl_bacteria.dmnd', help='DIAMOND reference database for gene alignment.')
    parser.add_argument('--cutoff', metavar='id_cutoff', default=97.0, help='Identity cutoff for a valid hit. Default: 97.0.')
    parser.add_argument('--plot', action='store_true', help='Plot histograms of gene length and ratio (if --align).')
    parser.add_argument('-s', '--stop', action='store_true', default=False, help='Database contains stop codons. Default: False.')
    parser.add_argument('-o', '--out', metavar='outdir', help='Output directory. Defaults to genlen_results in the same directory as the contigs file')
    parser.add_argument('-t,', '--threads', metavar='threads', default=16, help='Number of compute threads. Default: 16.')

    args = parser.parse_args()

    # Determine file names
    if not args.out:
        dir = os.path.split(args.contigs)[0]
        args.out = os.path.join(dir, "genlen_results")

    prefix = os.path.splitext(os.path.split(args.contigs)[1])[0]
    outfix = os.path.join(args.out, prefix)

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    # Call genes with prodigal
    prodigal_output = subprocess.check_output(f'prodigal -c -i {args.contigs} -o {outfix}.gff -a {outfix}.faa -d {outfix}.fna -f gff',shell=True).decode()

    # Read in protein sequences
    proteins = SeqIO.parse(f'{outfix}.faa', 'fasta')

    if args.align:
        # Blast protein sequences against reference database
        hits = align_genes(f'{outfix}.faa', args.db, args.threads)

        # Filter hits
        good_hits = {k:v for k,v in hits.items() if float(v[4]) >= args.cutoff}

        # Subtract stop codons from lengths then calculate length ratio
        for hit in hits.values():
            hit[2] = int(hit[2])-1
            if args.stop:
                hit[3] = int(hit[3])-1
            else:
                hit[3] = int(hit[3])
            hit.append(hit[2]/hit[3])
    else:
        hits = {x.id:[x.id, x.id, len(x)] for x in proteins}

    # Output
    with open(f'{outfix}_results.txt', 'w') as fo:
        for protein in proteins:
            try:
                towrite = '\t'.join(map(str, good_hits[protein.id]))
                fo.write(f'{towrite}\n')
            except:
                fo.write(f'{protein.id}\tNA\t{len(protein)}\tNA\tNA\tNA\n')

    # Graphs
    if args.plot:
        plot_hist(outfix, "len", [x[2] for x in hits.values()])
        if args.align:
            plot_hist(outfix, "ratio", [x[2]/x[3] for x in hits.values()])

#if __name__ == "__main__":
#    __main__()

