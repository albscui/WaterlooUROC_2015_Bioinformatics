from Bio import SeqIO
import fileinput
import csv
from peptide import *
import scipy.stats
from collections import Counter

def get_spectra():

	mass_intensity = [ ]
	header_data    = {}
	reading_header = False

	for line in fileinput.input('../../data/test.mgf'):

		line = line[:-1]

		# Spectra start
		if (line == "BEGIN IONS"):
			reading_header = True

		# Spectra end
		elif (line == "END IONS"):

			yield header_data, mass_intensity

			# Empty header data
			mass_intensity = [ ]
			header_data    = {}
			reading_header = False

		# heading in header data
		elif reading_header:

			for i in xrange(len(line)):
				if line[i] == '=':

					key = line[:i]
					val = line[i+1:]
					header_data[key] = val
					break

			if key == "RTINSECONDS":
				reading_header = False

		elif line != '':

			split = line.split(' ')
			mass = float(split[0])
			intensity = float(split[1])
			mass_intensity.append( (mass, intensity) )

def rescale(mass_intensity):

	# Get maximum value
	max_value = float(mass_intensity[len(mass_intensity)-1][0])

	for i in xrange(len(mass_intensity)):
		mass, intensity = mass_intensity[i]
		mass_intensity[i] = (mass, intensity/max_value )

def trypsinCut(sequence):
    peptides = []
    nTerm = 0
    cTerm = 0
    for aa in sequence:
        if (aa == "K" or aa == "R") and (cTerm+1) < len(sequence) and sequence[cTerm+1] != "P":
            pep = sequence[nTerm:cTerm+1]
            # print pep
            peptides.append(pep)
            nTerm = cTerm+1
        cTerm += 1

    return peptides

def createAminoAcidSeq(sequence, masstable):
    amino_acids = []
    for aa in sequence:
        thisMass = calcMass(aa, masstable)
        amino_acid = Amino_Acid(aa, thisMass)
        amino_acids.append(amino_acid)

    return amino_acids

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

def getMeasuredMassArray(massIntensities):
	masses = []
	for i in xrange(len(massIntensities)):
		masses.append( massIntensities[i][0] )

	return masses

if __name__ == '__main__':

	amino_acid_mass_table = readAminoAcidTable('../../data/aaMasses.csv')

	fasta_sequences = SeqIO.parse(open('../../data/ups.fasta'),'fasta')
	peptideDict = {}
	for protein in fasta_sequences:
		allPeptides = trypsinCut(protein.seq)
		for pepseq in allPeptides:
			if pepseq not in peptideDict:
				aminos = createAminoAcidSeq(pepseq, amino_acid_mass_table)
				peptide = Peptide(aminos)
				peptideDict[pepseq] = peptide
			else:
				peptide = peptideDict[pepseq]
			if protein.id not in peptide.proteins:
				peptide.proteins.append(protein.id)


	bestPeptides = {}
	for header, mass_intensity in get_spectra():

		peptide_mass = float(header['PEPMASS'])
		z = float(header['CHARGE'].split('+')[0])
		m = (peptide_mass - 1.007)*z - 18.01

		y_masses = []
		measuredMasses = getMeasuredMassArray(mass_intensity)

		spectra_peptide = None
		biggest_p_val = 0

		for key in peptideDict:
			peptide = peptideDict[key]
			if peptide not in bestPeptides:
				peptide_mass = peptide.mass
				peptide_str = peptide.sequence()

				if abs(m - peptide_mass) <= 0.5:

					for i in xrange(len(peptide_str)):

						prefix = peptide_str[:i]
						suffix = peptide_str[i+1:]

						m_i = calcMass( suffix, amino_acid_mass_table )
						y_masses.append(m_i)

					mannwhitneyu = scipy.stats.mannwhitneyu(measuredMasses, y_masses, use_continuity=False)

					p_val = mannwhitneyu[1]
					if p_val > biggest_p_val:
						biggest_p_val = p_val
						spectra_peptide = peptide

			if spectra_peptide != None:
				bestPeptides[spectra_peptide] = 1

	bestProteins = Counter()
	for key in bestPeptides:
		proteins = key.proteins
		for protein in proteins:
			bestProteins[protein] += 1
	print bestProteins
