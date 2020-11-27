# -*- coding: utf-8 -*-
# @Time    : 2020-10-19 16:52
# @Author  : xiaorui su
# @Email   :  suxiaorui19@mails.ucas.edu.cn
# @File    : preprocess_datad.py
# @Software : PyCharm

#data process for drug bank

# In[1]:
# !/usr/bin/python3

import untangle
import pandas as pd
import numpy as np
import os
from yoctol_utils.hash import consistent_hash


from config import RAW_DATA_DIR,PROCESSED_DATA_DIR

save_path_ddi = os.path.join(PROCESSED_DATA_DIR, 'drugbank_ddi')

def load_drugbank():
    filename = os.path.join(RAW_DATA_DIR, 'drugbank', 'drugbank.xml')

    # filename = "drugbank_all_full_database.xml"
    obj = untangle.parse(filename)

    return obj



def extract_smile(obj):

    # Building dataframe of chemical descriptors
    # Data Frame of DrugBank Small Molecule Type Drugs
    df_drugbank_sm = pd.DataFrame(
        columns=["drugbank_id", "name", "cas", "smiles", "logP ALOGPS", "logP ChemAxon", "solubility ALOGPS",
                 "pKa (strongest acidic)", "pKa (strongest basic)"])

    print(df_drugbank_sm)

    i = -1

    # iterate over drug entries to extract information
    for drug in obj.drugbank.drug:
        drug_type = str(drug["type"])

        # select for small molecule drugs
        if drug_type in ["small molecule", "Small Molecule", "Small molecule"]:
            i = i + 1

            # Get drugbank_id
            for id in drug.drugbank_id:
                if str(id["primary"]) == "true":
                    df_drugbank_sm.loc[i, "drugbank_id"] = id.cdata
            # Drug name
            df_drugbank_sm.loc[i, "name"] = drug.name.cdata

            # Drug CAS
            df_drugbank_sm.loc[i, "cas"] = drug.cas_number.cdata

            print(drug.drug_interactions.cdata)

            # Get SMILES, logP, Solubility
            # Skip drugs with no structure. ("DB00386","DB00407","DB00702","DB00785","DB00840",
            #                                            "DB00893","DB00930","DB00965", "DB01109","DB01266",
            #                                           "DB01323", "DB01341"...)
            if len(drug.calculated_properties.cdata) == 0:  # If there is no calculated properties
                continue
            else:
                for property in drug.calculated_properties.property:
                    if property.kind.cdata == "SMILES":
                        df_drugbank_sm.loc[i, "smiles"] = property.value.cdata
                        print(property.value.cdata)

                    if property.kind.cdata == "logP":
                        if property.source.cdata == "ALOGPS":
                            df_drugbank_sm.loc[i, "logP ALOGPS"] = property.value.cdata
                        if property.source.cdata == "ChemAxon":
                            df_drugbank_sm.loc[i, "logP ChemAxon"] = property.value.cdata

                    if property.kind.cdata == "Water Solubility":
                        df_drugbank_sm.loc[i, "solubility ALOGPS"] = property.value.cdata

                    if property.kind.cdata == "pKa (strongest acidic)":
                        df_drugbank_sm.loc[i, "pKa (strongest acidic)"] = property.value.cdata

                    if property.kind.cdata == "pKa (strongest basic)":
                        df_drugbank_sm.loc[i, "pKa (strongest basic)"] = property.value.cdata



    print(df_drugbank_sm.head(10))
    print(df_drugbank_sm.shape)

    # Drop drugs without SMILES from the dataframe
    df_drugbank_smiles = df_drugbank_sm.dropna()
    df_drugbank_smiles = df_drugbank_smiles.reset_index(drop=True)
    print(df_drugbank_smiles.shape)


    df_drugbank_smiles.head()
    df_drugbank_smiles.to_csv("dataset/kegg/drug_smile.csv")

# In[9]:

# write to csv
#df_drugbank_smiles.to_csv("drugbank_smiles.csv", encoding='utf-8', index=False)


#extract drug-drug interactions
def extract_ddi(obj):
    print("ddi")
    ddis = []

    for drug in obj.drugbank.drug:

        drug_type = str(drug["type"])
        drug_group = drug.groups
        flag = 0
        #select approved interactions and drug
        for groups in drug_group:
            for i in groups.group:
                if i.cdata == "approved":
                    print(i.cdata)
                    flag = 1
                    break


        # select for small molecule drugs and approved
        if drug_type in ["small molecule", "Small Molecule", "Small molecule"]:
            #i = i + 1
            # Get drugbank_id
            for id in drug.drugbank_id:
                if str(id["primary"]) == "true" and flag == 1:
                    print(id.cdata)
                    for inters in drug.drug_interactions:
                        print("interactions")
                        if len(inters) :
                            for i in inters.drug_interaction:
                                ddis.append([id.cdata,i.drugbank_id.cdata])

    #print(ddis)
    ddis = np.array(ddis)
    print(ddis.shape)
    np.save(save_path_ddi,ddis)
    #drug_all = np.concatenate((np.array(ddis[:,0]),np.array(ddis[:,1])),axis=0)
    #print(len(np.unique(drug_all)))

def match_drugid_smile():
    drugbank_ddi = np.load("dataset/kegg/drugid.npy")
    #drug_smile = np.loadtxt("dataset/kegg/drug_smile.csv",str,delimiter=",",skiprows=1)
    drugbank_ddi = pd.DataFrame(drugbank_ddi)
    print(drugbank_ddi)
    #print(drug_smile)

    f = open("dataset/kegg/drug_smile.csv", 'r')
    #drug_smile = f.readlines()
    #drug_smile = np.array(drug_smile)
    match_id_smile = []
    for line in f.readlines():
        row = line.strip().split(",")
        drugid = row[1]
        drugsmile = row[4]
        #print(row[1],row[4])
        samples = np.array(drugbank_ddi[drugbank_ddi[0].isin([str(drugid)])])
        if len(samples) == 0:
            continue
        else:
            for i in samples:
                print(i)
                match = [int(i[1]),row[4]]
                #print("match")
                #print(match)
                match_id_smile.append(match)

    results = np.array(match_id_smile)
    print(results[:,0])

    smiles = []
    hashes = []

    f = open("raw_data/kegg/entity2id.txt",'r')
    flag = 1
    for line in f.readlines():
        if flag == 1:
            flag = 0
            continue
        #print(line)
        row = line.strip().split("\t")
        drugid = str(row[1])
        itemindex = [index for index in range(0,len(results)) if results[index,0] == drugid]
        print(itemindex)
        if len(itemindex) == 0:
            smiles.append([" "])
            hash_vec = np.zeros(512)
            hashes.append(hash_vec)
        else:
            strings = []
            for chars in range(0,len(results[itemindex[0], 1])):
                strings.append(results[itemindex[0], 1][chars])
            print(strings)
            hashes.append(hash_seq(strings,512))
            smiles.append(strings)


    print(len(smiles))
    #print(results.shape) #(2014,2) finally, 2014 smiles for 2014 drugs
    #print(results)
    #print(smiles) #smiles for each entity
    print(hashes)
    hashes = np.array(hashes)
    hashes.reshape(len(smiles),512)



    np.save("data/drugbank_id_smile",smiles)
    np.save("data/drugbank_smile_hash", hashes)

    return results


def hash_seq(sequence, max_index):
    print(sequence)
    return np.array([consistent_hash(word) % max_index + 1 for word in sequence])




#extract_ddi(obj)
#extract_smile(obj)
#match_drugid_smile()
hashes = np.load("data/drugbank_smile_hash.npy",allow_pickle=True)
print(hashes.shape)
#hashes.reshape(len(hashes),len(hashes[0]))
new_hash = []
for i in hashes:
    if len(i) == 512:
        print("yes")
        new_hash.append(i)
    else:
        print(i)
        print(len(i))
        padding = np.zeros(512-len(i))
        new_hashes = np.concatenate((padding,i),axis=0)
        new_hash.append(new_hashes)

new_hash = np.array(new_hash)
print(np.array(new_hash).shape)
print(new_hash.reshape(len(hashes),512))
np.save("data/drugbank_smile_hash", new_hash)





#print(np.load("data/drugbank_smile_hash.npy",allow_pickle=True))

