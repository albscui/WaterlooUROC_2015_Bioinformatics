from Bio import SeqIO

def trypsinCut(sequence):
    peptides = []
    nTerm = 0
    cTerm = 0
    for aa in sequence:
        if (aa == "K" or aa == "R") and sequence[cTerm+1] != "P":
            pep = sequence[nTerm:cTerm+1]
            # print pep
            peptides.append(pep)
            nTerm = cTerm+1
        cTerm += 1

    return peptides


fasta_sequences = SeqIO.parse(open('../../data/ups.fasta'),'fasta')
for sequence in fasta_sequences:
    print repr(sequence.seq[0])
