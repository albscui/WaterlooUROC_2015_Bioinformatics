from Bio import SeqIO
aacidMassTable = {"A":}
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        new_sequence = some_function(sequence)
        write_fasta(out_file)

def pepCut(string):
    peptides = []
    positions
