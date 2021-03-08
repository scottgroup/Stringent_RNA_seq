import pandas as pd
import os

inputs = snakemake.input.counts
output = snakemake.output.merge_tpm

def get_name(file):
    name = file.split('/')[-1].split('.')[0]
    return name


df = pd.read_csv(inputs[0], sep="\t")
df = df[['gene_id','gene_name']]

names = [get_name(x) for x in inputs]

for file, name in zip(inputs, names):
    tmp = pd.read_csv(file, sep="\t")
    df[name] = df.gene_id.map(dict(zip(tmp.gene_id, tmp.tpm)))

df.to_csv(output,sep="\t",index=False)
