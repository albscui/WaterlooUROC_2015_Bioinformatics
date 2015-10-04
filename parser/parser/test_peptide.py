from peptide import *
import csv
from read_in_data import *

massTable = readAminoAcidTable('../../data/aaMasses.csv')


def createAminoAcidSeq(sequence):
    amino_acids = []
    for aa in sequence:
        thisMass = calcMass(aa, massTable)
        amino_acid = Amino_Acid(aa, thisMass)
        amino_acids.append(amino_acid)

    return amino_acids

sequence1 = "AHKSEVAHRFKDLGEENFKALVL"
aminoacids1 = createAminoAcidSeq(sequence1)
peptide1 = Peptide(aminoacids1, "someProtein")
print peptide1.sequence()
