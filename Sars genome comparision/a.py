import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna,generic_rna,generic_protein
from Bio import SeqIO
import pandas as pd
from collections import Counter
from Bio.PDB import PDBParser
import nglview as nv

for record in SeqIO.parse("sequence.fasta","fasta"):
	ncov_record=SeqIO.read("sequence.fasta","fasta")
	ncov_dna=ncov_record.seq

	# print(ncov_dna)
	print(len(ncov_dna))
#Transcription to mrna
ncov_mrna=ncov_dna.transcribe()


#protein translation
ncov_protein=ncov_dna.translate()
print(len(ncov_protein))

ncov_aa = ncov_protein.split("*")
ncov_clean=[str(i)for i in ncov_aa ]
# print(ncov_aa)
print(len(ncov_aa))
print(ncov_clean)
# print(len(ncov_clean))

df = pd.DataFrame({'amino_acids':ncov_clean})
df['count']=df['amino_acids'].str.len()
df.head()

#largest seq before *
df.nlargest(10,"count")

#count frequency of AA
Counter(ncov_protein).most_common(10)
print(df['count'])
print(Counter(ncov_protein).most_common(10))


parser=PDBParser()
structure= parser.get_structure("6LU7","6LU7.pdb")
model=structure[0]
for chain in model:
	print(chain)

view = nv.show_biopython(structure)
view
