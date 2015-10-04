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
