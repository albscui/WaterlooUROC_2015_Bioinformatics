import fileinput
import numpy


def get_spectra():

	mass_intensity = [ ]
	header_data    = {}
	reading_header = False

	for line in fileinput.input('test.mgf.txt'):
		
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

if __name__ == '__main__':
	

	generator = get_spectra()
	for header, mass_intensity in generator:
		print header

	
