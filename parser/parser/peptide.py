class Peptide:

    def __init__(self,amino_acids):
        self.mass = sum([ a.mass for a in amino_acids ])
        self.amino_acids = amino_acids
        self.proteins = []

    def sequence(self):
        sequence = ''.join([ a.letter for a in self.amino_acids ])
        return sequence


class Amino_Acid:
    def __init__(self, letter, mass):
        self.letter = letter
        self.mass = mass
