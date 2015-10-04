import csv

def readAminoAcidTable(csvFile):
    aa_mass = {}
    with open(csvFile, 'r') as aaCSVFile:
        reader = csv.DictReader(aaCSVFile)
        for row in reader:
            aa_mass[row['Amino-acid']] = row['Monoisotopic']

    return aa_mass

def calcMass(peptide, mass_table):
    mass = 0
    for aa in peptide:
        mass += float(mass_table[aa])

    return mass
