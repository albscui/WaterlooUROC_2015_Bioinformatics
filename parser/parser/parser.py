from bio import SeqIO

fasta_sequences = SeqIO.parse(open('../../data/ups.fasta'),'fasta')
for protein in fasta_sequences:
    allPeptides = trypsinCut(protein.seq)
    print allPeptides
