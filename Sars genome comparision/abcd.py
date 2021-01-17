from Bio.Blast import NCBIWWW , NCBIXML

fasta_string = open("n-covid genome.fna").read()
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)
blast_record = NCBIXML.read(result_handle)

len(blast_record.alignments)
E_VALUE_THRESH = 0.01
for alignment in blast_record_alignments:
	for hsp in alignment.hsp:
		if hsp.expect < E_VALUE_THRESH:
			print("***Alignment***")
			print("sequence:", alignment.title)
			print("length:",alignment.length)
			print("e value:",hsp.expect)
			print(hsp.query)
			print(hsp.match)
			print(hsp.sbjct)