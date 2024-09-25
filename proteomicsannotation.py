import pandas as pd
from Bio import Entrez
from Bio import SeqIO

genebank = pd.read_table("C:/Users/ginlo/Downloads/genebank_failed.txt", sep="\t", header=0)

genebank['ID_Numbers'] = genebank['ID_Numbers'].astype(str)

Entrez.email = '[emailhere]'
Entrez.api_key = '[api key here]'

for i in genebank["ID_Numbers"]:
    handle = Entrez.efetch(db="protein", id=i, rettype="gb", retmode="text")
    x = SeqIO.read(handle, 'genbank')
    print(x.description)
    with open("C:/Users/ginlo/Downloads/gene_function2.txt", "a") as newfile:
        newfile.write(f'\n{x.description}')
