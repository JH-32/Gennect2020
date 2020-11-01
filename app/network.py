import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

"""
    1. Crawl 결과에서 positive 텀과 negative 텀이 포함된 문장을 다시 골라내어 이들을 network 형태로 저장 
    2. Jacard distance를 계산하여 relation을 prediction 하고, 
    3. 특정 threshold를 넘는 relation만 선정하여 pickle 형태로 저장 
"""

tt = pd.read_excel('static/data/Crawl.xlsx')

positives = ["prolong", "regulator", "accelerate", "modulate", "enhance", "promote", "interact", "up-regulate"]
negatives = ["down-regulate", "suppress", "prevent", "inhibit", "exploit"]

posis = []
negas = []
reDic = {}
for index, row in tt.iterrows():
    ## positive 또는 negative term을 가지고 있는 문장일 경우 저장 
    geneA = row[0]
    geneB = row[1]
    string = row[2]
    pmid = row[3]
    if any(ele for ele in positives if(ele in string)):
        re = [geneA, geneB]
        posis.append(re)
        reCur = reDic.get(geneA, {})
        reCC = reCur.get(geneB, [])
        if [string, pmid] not in reCC:
            reCC.append([string, pmid])
        reCur[geneB] = reCC
        reDic[geneA] = reCur
        reCur = reDic.get(geneB, {})
        reCC = reCur.get(geneA, [])
        if [string, pmid] not in reCC:
            reCC.append([string, pmid])
        reCur[geneA] = reCC
        reDic[geneB] = reCur
    if any(ele for ele in negatives if(ele in string)):
        re = [geneA, geneB]
        negas.append(re)
        reCur = reDic.get(geneA, {})
        reCC = reCur.get(geneB, [])
        if [string, pmid] not in reCC:
            reCC.append([string, pmid])
        reCur[geneB] = reCC
        reDic[geneA] = reCur
        reCur = reDic.get(geneB, {})
        reCC = reCur.get(geneA, [])
        if [string, pmid] not in reCC:
            reCC.append([string, pmid])
        reCur[geneA] = reCC
        reDic[geneB] = reCur

with open('dict.pickle', 'wb') as handle:
    pickle.dump(reDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

## pickle에 결과 저장, 나중에 flask에서 결과를 보여주기 위해 사용 

G = nx.Graph()
## networkx 객체 생성 

for edge in posis:
    ## positive edge는 weight를 1씩 더해주고 
    if edge in G.edges():
        G[edge[0]][edge[1]]['weight'] += 1
    else:
        G.add_edge(edge[0], edge[1], weight = 1)

for edge in negas:
    ## negative edge는 weight를 1씩 감소시킴 
    if edge in G.edges():
        G[edge[0]][edge[1]]['weight'] -= 1
    else:
        G.add_edge(edge[0], edge[1], weight = -1)

nx.write_gpickle(G, "Whole.pickle")
## network 저장 

preds = nx.jaccard_coefficient(G) ## Jacard distance 계산 
count = 0
ps = []
epred = []
res = []
for u, v, p in preds:
    count += 1 
    if count % 1000000 == 0:
        print(count, "edges calculated")
    if p > 0:
        res.append([u, v, p])

ress = pd.DataFrame(res, index = None)
ress.to_csv("Jacard.csv", header = None, index = None)

J = nx.Graph()
for index, row in ress.iterrows():
    if row[2] > 0.05:  ## Jacard distance가 0.05 이상인 경우를 prediction 
        J.add_edge(row[0], row[1], dist = row[2])

nx.write_gpickle(J, "Jacard.pickle")