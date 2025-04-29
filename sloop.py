import os
import subprocess

if __name__ == "__main__":
    for n in ['774','59','112','685','1033','48','414']:
        queryfile = f'query/MI_KPC_{n}/'
        subjectfile = f'reference/kpc_assemblies/MI_KPC_{n}_flye_medaka_polypolish.fasta'
        highlightfile = 'highlight_data_KPC.tsv'
        command = ['python','dotplot.py','--query',queryfile,'--subject',subjectfile,'--name','kpc_v1','--highlight',highlightfile]
        subprocess.call(command)







