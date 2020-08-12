import argparse, subprocess
from Bio import Seq, SeqIO
from BCBio import GFF
from io import StringIO

def __main__():
    parser = argparse.ArgumentParser(description='Call genes with prodigal, calculate gene length distribution and compare to the nearest BLAST hits')
    parser.add_argument('contigs',metavar='contigs_file',help='File containing contigs in fasta format')
    parser.add_argument('--gff',metavar='gff_file',help='Save predicted genes to given gff file')

    args = parser.parse_args()

    # Call genes with prodigal
    if args.gff is None:
        prodigal_output = subprocess.check_output('prodigal -i {} -f gff'.format(args.contigs),shell=True).decode()
        prodigal_handle = StringIO(prodigal_output)
        records = GFF.parse(prodigal_handle)
    else:
        prodigal_output = subprocess.check_output('prodigal -i {} -o {} -f gff'.format(args.contigs,args.gff),shell=True).decode()
        records = GFF.parse(args.gff)

    # Determine gene lengths
    gene_lengths = [len(feature) for record in records for feature in record.features]
    gene_lengths = sorted(gene_lengths,reverse=False)

    # Output gene lengths to file
    with open('gene_lengths.txt','w') as fo:
        for gl in gene_lengths:
            fo.write('{}\n'.format(gl))

if __name__ == "__main__":
    __main__()

