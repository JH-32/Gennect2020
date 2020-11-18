from flask import Flask, render_template, request, url_for, send_file
from flask_bootstrap import Bootstrap
import networkx as nx 
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import pickle
import os
from io import BytesIO, StringIO
from google.cloud import storage

app = Flask(__name__)

## 미리 데이터 로드해두기 
with open('static/data/dict.pickle', 'rb') as handle:
    resData = pickle.load(handle)

Bootstrap(app) 



@app.route('/') #Main URL
def main():
	return render_template('main.html') 

@app.route('/result', methods=['POST']) # 결과 URL
def result():
    tGene = request.form['gene']

    G = nx.read_gpickle("static/data/Whole.pickle")

    print("Set genes for visualization..")
    target = [tGene]
    target.extend(list(G.neighbors(target[0])))
    H = G.subgraph(target)
    del G 

    ## positive, negative, neutral link에 대해 각각 시각화 
    eposi = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] > 0]
    eposiweight = [d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] > 0]
    resTablePosi = []
    for [u, v] in eposi:
        try:
            for [sen, pmid] in resData[u][v]:
                resTablePosi.append([u, v, sen, pmid])
        except KeyError:
            continue

    enega = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] < 0]
    print(enega)
    enegaweight = np.negative([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] < 0])
    resTableNega = []
    for [u, v] in enega:
        try:
            for [sen, pmid] in resData[u][v]:
                resTableNega.append([u, v, sen, pmid])
        except KeyError:
            continue
    
    eneu = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] == 0]
    print(eneu)
    resTableNeu = []
    for [u, v] in eneu:
        try:
            for [sen, pmid] in resData[u][v]:
                resTableNeu.append([u, v, sen, pmid])
        except KeyError:
            continue

    timestr = time.strftime("%Y%m%d_%H%M%S")
    figName = "static/fig/" + tGene + timestr + ".png"
    plt.figure(figsize=(20, 20))
    pos = nx.circular_layout(H)
    nx.draw_networkx_nodes(H, pos, node_size=20)
    nx.draw_networkx_labels(H, pos, font_size=30, font_family="sans-serif")
    nx.draw_networkx_edges(H, pos, edgelist = eposi, width=eposiweight, edge_color="r")
    nx.draw_networkx_edges(H, pos, edgelist = enega, width=enegaweight, edge_color="b")
    nx.draw_networkx_edges(H, pos, edgelist = eneu, width=5, edge_color="y")

    plt.savefig(figName, format="PNG", dpi=300)

    return render_template('result.html', fName = figName, tposi = resTablePosi, tnega = resTableNega, tneu = resTableNeu)


@app.route('/result_ai', methods=['POST']) # 결과 URL
def resultai():
    tGene = request.form['gene']

    G = nx.read_gpickle("static/data/Whole.pickle")

    with open("static/data/pickle_model.pkl", 'rb') as file:
        pickle_model = pickle.load(file)

    print("Set genes for visualization..")
    target = [tGene]
    target.extend(list(G.neighbors(target[0])))
    target_edges = []
    GAS = []
    GBS = []
    for GeneA in target:
        for GeneB in target:
            if GeneA == GeneB:
                continue
            target_edges.append((GeneA, GeneB))
            GAS.append(GeneA)
            GBS.append(GeneB)
    H = G.subgraph(target)

    pref_train = list(nx.preferential_attachment(G, target_edges))
    jacc_train = list(nx.jaccard_coefficient(G, target_edges))

    df_train = []
    for i in range(len(target_edges)):
        cur = [pref_train[i][2], jacc_train[i][2]]
        df_train.append(cur)
    
    df_train = pd.DataFrame(df_train, columns = ["pref", "jacc"])
    del G 

    pred = pickle_model.predict(df_train)
    prob = pickle_model.predict_proba(df_train)[:, 1]

    df_train["GeneA"] = GAS
    df_train["GeneB"] = GBS
    df_train["Prob"] = prob
    df_train["pred"] = pred

    df_train = df_train[df_train["pred"] == 1]

    epredf = []
    for i in range(len(target_edges)):
        if pred[i] == 1:
            epredf.append(target_edges[i])

    ## positive, negative, neutral link에 대해 각각 시각화 
    eposi = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] > 0]
    eposiweight = [d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] > 0]
    resTablePosi = []
    for [u, v] in eposi:
        try:
            for [sen, pmid] in resData[u][v]:
                resTablePosi.append([u, v, sen, pmid])
        except KeyError:
            continue

    enega = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] < 0]
    enegaweight = np.negative([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] < 0])
    resTableNega = []
    for [u, v] in enega:
        try:
            for [sen, pmid] in resData[u][v]:
                resTableNega.append([u, v, sen, pmid])
        except KeyError:
            continue
    
    eneu = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] == 0]
    resTableNeu = []
    for [u, v] in eneu:
        try:
            for [sen, pmid] in resData[u][v]:
                resTableNeu.append([u, v, sen, pmid])
        except KeyError:
            continue

    timestr = time.strftime("%Y%m%d_%H%M%S")
    figName = "static/fig/" + tGene + timestr + ".png"
    plt.figure(figsize=(20, 20))
    pos = nx.circular_layout(H)
    nx.draw_networkx_nodes(H, pos, node_size=20)
    nx.draw_networkx_labels(H, pos, font_size=30, font_family="sans-serif")
    nx.draw_networkx_edges(H, pos, edgelist = eposi, width=eposiweight, edge_color="r")
    nx.draw_networkx_edges(H, pos, edgelist = enega, width=enegaweight, edge_color="b")
    nx.draw_networkx_edges(H, pos, edgelist = eneu, width=5, edge_color="y")
    nx.draw_networkx_edges(H, pos, edgelist = epredf, edge_color="g", style = "dashed")

    plt.savefig(figName, format="PNG", dpi=300)

    return render_template('resultai.html', fName = figName, tpred = df_train.values.tolist(), tposi = resTablePosi, tnega = resTableNega, tneu = resTableNeu)


@app.route('/result_aionly', methods=['POST']) # 결과 URL
def resultaionly():
    tGene = request.form['gene']

    G = nx.read_gpickle("static/data/Whole.pickle")

    with open("static/data/pickle_model.pkl", 'rb') as file:
        pickle_model = pickle.load(file)

    print("Set genes for visualization..")
    target = [tGene]
    target.extend(list(G.neighbors(target[0])))
    target_edges = []
    GAS = []
    GBS = []
    for GeneA in target:
        for GeneB in target:
            if GeneA == GeneB:
                continue
            target_edges.append((GeneA, GeneB))
            GAS.append(GeneA)
            GBS.append(GeneB)
    H = G.subgraph(target)

    pref_train = list(nx.preferential_attachment(G, target_edges))
    jacc_train = list(nx.jaccard_coefficient(G, target_edges))

    df_train = []
    for i in range(len(target_edges)):
        cur = [pref_train[i][2], jacc_train[i][2]]
        df_train.append(cur)
    
    df_train = pd.DataFrame(df_train, columns = ["pref", "jacc"])
    del G 

    pred = pickle_model.predict(df_train)
    prob = pickle_model.predict_proba(df_train)[:, 1]

    df_train["GeneA"] = GAS
    df_train["GeneB"] = GBS
    df_train["Prob"] = prob
    df_train["pred"] = pred

    df_train = df_train[df_train["pred"] == 1]
    print(df_train)

    print(target_edges)
    print(prob)
    print(pred)
    epredf = []
    for i in range(len(target_edges)):
        if pred[i] == 1:
            print(prob[i])
            epredf.append(target_edges[i])

    print(epredf)

    timestr = time.strftime("%Y%m%d_%H%M%S")
    figName = "static/fig/" + tGene + timestr + ".png"
    plt.figure(figsize=(20, 20))
    pos = nx.circular_layout(H)
    nx.draw_networkx_nodes(H, pos, node_size=20)
    nx.draw_networkx_labels(H, pos, font_size=30, font_family="sans-serif")
    nx.draw_networkx_edges(H, pos, edgelist = epredf, edge_color="g", style = "dashed")

    plt.savefig(figName, format="PNG", dpi=300)

    return render_template('resultaionly.html', fName = figName, tpred = df_train.values.tolist())


if __name__ == '__main__':
	app.run(host='0.0.0.0', debug=True)
