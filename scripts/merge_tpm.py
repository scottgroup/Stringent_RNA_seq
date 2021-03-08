import pandas as pd 
import os 


def string(i):
    name,tsv = os.path.splitext(i)
    name = os.path.basename(str(name))
   
    return name
name_list = {
    "brain_A802070_S5": "Liver_1",
    
    "brain_B702085_S4": "Liver_2",
    
    "brain_C210018_S3": "Liver_3",
    
    "BRN117": "Breast_1",
   
    "BRN491": "Breast_2",
    
    "BRN_574": "Breast_3",
    
    
    "liver_B603125_S2": "Brain_1",
   
    "liver_B802097_S1": "Brain_2",
    
    "NPR026": "Prostate_1",
    
    "NPR029": "Prostate_2",
    
    "NPR036": "Prostate_3",
    
    "OVN218": "Ovary_1",
    
    "OVN377": "Ovary_2",
    
    "OVN85": "OvARY_3",
   
    "Skeletal_muscle_1_S3": "Skeletal_muscle_1",
    
    "Skeletal_muscle_2_S4": "Skeletal_muscle_2",
    
    "Skeletal_muscle_3_S5": "Skeletal_muscle_3",
    
    "Testis1_S1": "Testis_3",

    "testis_B606080_S7": "Brain_1",
    
    "testis_B703033_S6": "Brain_2",
    
    "testis_B703034_S8": "Brain_3",
    



  }




inputs = snakemake.input.counts
#passing = snakemake.params.passing
#inputs = [x for x in inputs if f"/{passing}/" in x]
print(inputs)

output = snakemake.output.merge_tpm

merge_file_test = pd.read_csv(inputs[0],sep="\t")
Empty_merge = merge_file_test[['gene_id','gene_name']]

liste_nom = ['gene_id','gene_name']
name_list = list(name_list.values())
for i in range(len(name_list)):    
    liste_nom.append(name_list[i])

Merge_final = []
Merge_final.append(Empty_merge)
for i in inputs:
    merge_file = pd.read_csv(i,sep="\t")
    
    transfer_data = merge_file['tpm']
    transfer_data = transfer_data.rename(columns={'tpm':str(string(i))})
    


    Merge_final.append(transfer_data)

Merge_final = pd.concat(Merge_final,axis=1)

Merge_final.columns= liste_nom



Merge_final.to_csv(output,sep="\t",index=False)






