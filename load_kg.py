# -*- coding: utf-8 -*-
# @Time    : 2020-10-20 15:39
# @Author  : xiaorui su
# @Email   :  suxiaorui19@mails.ucas.edu.cn
# @File    : load_kg.py
# @Software : PyCharm

d_name = "ogbl-biokg"
import dgl
import numpy as np
import torch as th
from ogb.linkproppred import DglLinkPropPredDataset
import random


dataset = DglLinkPropPredDataset(name = d_name)
print(dataset[0])
split_edge = dataset.get_edge_split()
#print(split_edge)
#five kinds of nodes: drug, protein,function,disease
train_edge, valid_edge, test_edge = split_edge["train"], split_edge["valid"], split_edge["test"]
graph = dataset[0] # dgl graph object containing only training edges

#head_type,head,relation,tail_type,tail
#ids: head,relation,tail


def construct_entity():
    ##construct entity2id.txt for biokg
    # 10687
    disease_entity = "dataset/ogbl_biokg/mapping/disease_entidx2name.csv"
    # 10533
    drug_entity = "dataset/ogbl_biokg/mapping/drug_entidx2name.csv"
    # 45085
    function_entity = "dataset/ogbl_biokg/mapping/function_entidx2name.csv"
    # 17499
    protein_entity = "dataset/ogbl_biokg/mapping/protein_entidx2name.csv"
    # 9969
    sideeffect_entity = "dataset/ogbl_biokg/mapping/sideeffect_entidx2name.csv"

    entity = [disease_entity, drug_entity, function_entity, protein_entity, sideeffect_entity]

    entityID = []

    for i in entity:
        with open(i, encoding='utf-8') as f:
            data = np.loadtxt(f, delimiter=",", skiprows=1, dtype='str')
            entityID.append(data)

    print(entityID)
    print(len(entityID))
    entity_number = 0
    entity2id = []
    for i in range(0, len(entityID)):
        for j in entityID[i]:
            samples = [j[1], entity_number]
            entity2id.append(samples)
            entity_number = entity_number + 1

    print(entity2id)
    np.savetxt("raw_data/ogb/entity2id.txt", entity2id, fmt='%s')  ##entity2d number 93773
    return ;






##constrcut train2id.txt for biokg [drugid,drugid,relationid]
##constrcut approved.txt for biokg [entity,entity,interact]

def construct_example():
    all_edge_name = [train_edge, valid_edge, test_edge]
    examples = []
    trains = []
    drug_base = 10687
    function_base = 10687 + 10533
    protein_base = 10687 + 10533 + 45085
    side_base = 10687 + 10533 + 45085 + 17499

    print(train_edge['head'][0].item() + 1000)

    for i in all_edge_name:

        for j in range(0, len(i['head_type'])):
            if i['head_type'][j] == 'drug' and i['tail_type'][j] == 'drug':
                samples = [i['head'][j].item() + drug_base, i['tail'][j].item() + drug_base, 1]
                sample = [i['head'][j].item() + drug_base, i['tail'][j].item() + drug_base, i['relation'][j].item()]
                examples.append(samples)
                trains.append(sample)
            elif i['head_type'][j] == 'disease' and i['tail_type'][j] == 'protein':
                samples = [i['head'][j].item(), i['tail'][j].item() + protein_base, i['relation'][j].item()]
                trains.append(samples)
            elif i['head_type'][j] == 'drug' and i['tail_type'][j] == 'disease':
                samples = [i['head'][j].item() + drug_base, i['tail'][j].item(), i['relation'][j].item()]
                trains.append(samples)
            elif i['head_type'][j] == 'drug' and i['tail_type'][j] == 'protein':
                samples = [i['head'][j].item() + drug_base, i['tail'][j].item() + protein_base, i['relation'][j].item()]
                trains.append(samples)
            elif i['head_type'][j] == 'drug' and i['tail_type'][j] == 'sideeffect':
                samples = [i['head'][j].item() + drug_base, i['tail'][j].item() + side_base, i['relation'][j].item()]
                trains.append(samples)
            elif i['head_type'][j] == 'function' and i['tail_type'][j] == 'function':
                samples = [i['head'][j].item() + function_base, i['tail'][j].item() + function_base,
                           i['relation'][j].item()]
                trains.append(samples)
            elif i['head_type'][j] == 'protein' and i['tail_type'][j] == 'function':
                samples = [i['head'][j].item() + protein_base, i['tail'][j].item() + function_base,
                           i['relation'][j].item()]
                trains.append(samples)
            elif i['head_type'][j] == 'protein' and i['tail_type'][j] == 'protein':
                samples = [i['head'][j].item() + protein_base, i['tail'][j].item() + protein_base,
                           i['relation'][j].item()]
                trains.append(samples)

    np.savetxt("raw_data/ogb/train2id.txt", np.array(trains), fmt='%s')
    np.savetxt("raw_data/ogb/approved_example.txt", np.array(examples), fmt='%s')
    print(len(trains))
    print(len(examples))

    return;

def generate_negative_samples():
    import pandas as pd
    examples_path = "raw_data/ogb/approved_example.txt"
    examples = np.loadtxt(examples_path,dtype='int64')
    examples_df = pd.DataFrame(examples)
    print(examples)
    drug_id_all = np.unique(np.concatenate((examples[:,0],examples[:,1]),axis=0))
    print(len(drug_id_all))#9131
    nagative_samples = []
    while len(nagative_samples) < len(examples):
        print(len(nagative_samples))
        random_drug = random.randint(10687,10687+10533)
        print(random_drug)

        random_drug1 = random.randint(10687,10687+10533)
        print(random_drug1)

        if random_drug1 in np.array(examples_df[examples_df[0].isin([random_drug])][1]):
            continue
        else:
            sample = [random_drug,random_drug1,0]
            nagative_samples.append(sample)

    print(nagative_samples)
    print(len(nagative_samples))
    all_example = np.concatenate((examples,np.array(nagative_samples)),axis=0)
    print(all_example)
    np.savetxt("raw_data/ogb/approved_example.txt", all_example, fmt='%d')




generate_negative_samples()









'''
import numpy as np

f = open('raw_data/kegg/entity2id.txt', 'r')
drugbank_id_kegg = []
flag = 0
for line in f.readlines():
    flag = flag + 1
    row = line.strip().split("\t")

    print(row)

    if "drugbank:" in str(row[0]):
        split_result = row[0].split(":")
        #print(row[0].split(":"))
        id = split_result[2][0:7]
        drugbank_id_kegg.append([id,row[1]])

    if "drugbank/" in str(row[0]):
        #print(row[0])
        split_result = row[0].split("/")
        print(split_result)
        id = split_result[4][0:7]
        drugbank_id_kegg.append([id, row[1]])


print(drugbank_id_kegg)
print(len(np.unique(np.array(drugbank_id_kegg)[:,0])))
np.save("dataset/kegg/drugid",drugbank_id_kegg)
#finally there are 1463 unique drugs in kegg dataset'''



