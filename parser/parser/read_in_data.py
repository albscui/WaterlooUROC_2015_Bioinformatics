import csv


with open('../../data/aaMasses.csv', 'r') as aaCSVFile:
    reader = csv.DictReader(aaCSVFile)
    for row in reader:
        print(row['Amino-acid'], row['Monoisotopic'])

def calcMass(peptide, mass_table):
    mass = 0
    for aa in peptide:
        mass += mass_table[aa]
