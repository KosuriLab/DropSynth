import argparse

def csv_reader(filename):
	'''
	Read library from .csv file. Assumes the first column are the sequence names
	and the second column is the sequence
	'''

	# dictionary to hold library, where names are keys and sequences are values
	library = {}

	# open input file from given input name
	infile = open(filename, 'rU')

	# for every line
	for line in infile.readlines():
		# strip the newline character and split by comma. The first item will
		# be the name and the second item the sequence
		name, seq = line.strip().split(',')
		# store in dictionary
		library[name] = seq.upper()

	return library

def csv_writer(library, output):
	'''
	Write a dictionary file to output file. Assumes the keys are sequence names
	and the values are the sequences
	'''

	# open output file to write with given output name
	outfile = open(output, 'w')

	# for each sequence
	for name in library:
		# grab corresponding sequence
		seq = library[name]
		# join the two fields with a comma, write to output
		outfile.write(','.join([name, seq]) +'\n')

	outfile.close()

def reverse_complement(seq):
	'''
	Returns the reverse complement of a sequence
	'''

	# reverse the sequence
	rev = seq[::-1]
	# create blank string
	rc = ''

	# for each letter in the reverse sequence
	for nt in rev:
		# add the complement
		if nt == 'A':
			rc += 'T'
		elif nt == 'T':
			rc += 'A'
		elif nt == 'C':
			rc += 'G'
		else: # nt  == 'G'
			rc += 'C'
	
	return rc

def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize.
	'''

	# get the reverse compliment of the oligo
	rc_oligo = reverse_complement(oligo)

	# count the number of A's for both strands
	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	# save the oligo with lower number of A's
	if oligo_As < rc_As:
		final_oligo = oligo
	else:
		final_oligo = rc_oligo

	return final_oligo


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Select strand with lowest A content\
									 from a .csv library file')
	parser.add_argument('library_file', help='.csv file of library, first column\
						 sequence name and second column the sequence')
	parser.add_argument('output_name', help='Name of output file of new library')

	# parse arguments from command line
	args = parser.parse_args()

	# read in library file
	library = csv_reader(args.library_file)

	# create blank dictionary to save new library
	new_library = {}

	for name in library:
		# for each sequence, grab the sequence and choose the strand with the
		# best (lowest) A content
		seq = library[name]
		best_oligo = best_A_content(seq)
		new_library[name] = best_oligo

	# write to file
	csv_writer(new_library, args.output_name)



