from flask import Flask, render_template, request, url_for, send_file
from flask_bootstrap import Bootstrap
import networkx as nx 
import time
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
    eposiweight = np.log([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] > 0])
    resTablePosi = []
    for [u, v] in eposi:
        try:
            for [sen, pmid] in resData[u][v]:
                resTablePosi.append([u, v, sen, pmid])
        except KeyError:
            continue

    enega = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] < 0]
    enegaweight = np.log(np.negative([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] < 0]))
    resTableNega = []
    for [u, v] in enega:
        for [sen, pmid] in resData[u][v]:
            resTableNega.append([u, v, sen, pmid])
    
    eneu = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] == 0]
    resTableNeu = []
    for [u, v] in eneu:
        for [sen, pmid] in resData[u][v]:
            resTableNeu.append([u, v, sen, pmid])

    timestr = time.strftime("%Y%m%d_%H%M%S")
    figName = "static/fig/" + tGene + timestr + ".png"
    plt.figure(figsize=(20, 20))
    pos = nx.circular_layout(H)
    nx.draw_networkx_nodes(H, pos, node_size=20)
    nx.draw_networkx_labels(H, pos, font_size=30, font_family="sans-serif")
    nx.draw_networkx_edges(H, pos, edgelist = eposi, width=eposiweight, edge_color="r")
    nx.draw_networkx_edges(H, pos, edgelist = enega, width=enegaweight, edge_color="b")
    nx.draw_networkx_edges(H, pos, edgelist = eneu, width=enegaweight, edge_color="y")

    plt.savefig(figName, format="PNG", dpi=300)

    return render_template('result.html', fName = figName, tposi = resTablePosi, tnega = resTableNega, tneu = resTableNeu)


@app.route('/result_ai', methods=['POST']) # 결과 URL
def resultai():
    tGene = request.form['gene']

    G = nx.read_gpickle("static/data/Whole.pickle")

    with open(pkl_filename, 'rb') as file:
        pickle_model = pickle.load(file)

    print("Set genes for visualization..")
    target = [tGene]
    target.extend(list(G.neighbors(target[0])))
    H = G.subgraph(target)

    pref_train = list(nx.preferential_attachment(G, target))
    jacc_train = list(nx.jaccard_coefficient(G, target))

    df_train = []
    for i in range(len(train)):
        cur = [pref_train[i][2], jacc_train[i][2]]
        df_train.append(cur)
    
    df_train = pd.DataFrame(df_train, columns = ["pref", "jacc"])
    del G 
    thre = 0.08



    ## 특정 jacard distance 이상의 값을 사용 
    epredf = [(u, v) for (u, v, d) in Jsub.edges(data=True) if d["dist"] > thre]
    resTablePred = []
    for [u, v] in epredf:
        resTablePred.append([u, v])

    print("Start prediction..")
    print("prediction completed")

    ## positive, negative, neutral link에 대해 각각 시각화 
    eposi = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] > 0]
    eposiweight = np.log([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] > 0])
    resTablePosi = []
    for [u, v] in eposi:
        try:
            for [sen, pmid] in resData[u][v]:
                resTablePosi.append([u, v, sen, pmid])
        except KeyError:
            continue

    enega = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] < 0]
    enegaweight = np.log(np.negative([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] < 0]))
    resTableNega = []
    for [u, v] in enega:
        for [sen, pmid] in resData[u][v]:
            resTableNega.append([u, v, sen, pmid])
    
    eneu = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] == 0]
    resTableNeu = []
    for [u, v] in eneu:
        for [sen, pmid] in resData[u][v]:
            resTableNeu.append([u, v, sen, pmid])

    timestr = time.strftime("%Y%m%d_%H%M%S")
    figName = "static/fig/" + tGene + timestr + ".png"
    plt.figure(figsize=(20, 20))
    pos = nx.circular_layout(H)
    nx.draw_networkx_nodes(H, pos, node_size=20)
    nx.draw_networkx_labels(H, pos, font_size=30, font_family="sans-serif")
    nx.draw_networkx_edges(H, pos, edgelist = eposi, width=eposiweight, edge_color="r")
    nx.draw_networkx_edges(H, pos, edgelist = enega, width=enegaweight, edge_color="b")
    nx.draw_networkx_edges(H, pos, edgelist = eneu, width=enegaweight, edge_color="y")
    nx.draw_networkx_edges(H, pos, edgelist = epredf, edge_color="g", style = "dashed")

    plt.savefig(figName, format="PNG", dpi=300)

    return render_template('resultai.html', fName = figName, tpred = resTablePred, tposi = resTablePosi, tnega = resTableNega, tneu = resTableNeu)


if __name__ == '__main__':
	app.run(host='0.0.0.0', debug=True)
