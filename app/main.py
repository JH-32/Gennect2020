from flask import Flask, render_template, request, url_for, send_file
from flask_bootstrap import Bootstrap
import networkx as nx 
import time
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
from io import BytesIO, StringIO

app = Flask(__name__)
Bootstrap(app) 

@app.route('/') #Main URL
def main():
	return render_template('main.html') 

@app.route('/result', methods=['POST']) # 결과 URL
def result():
    tGene = request.form['gene']

    G = nx.read_gpickle("static/data/Whole.pickle")
    J = nx.read_gpickle("static/data/Jacard.pickle")

    print("Set genes for visualization..")
    target = [tGene]
    target.extend(list(G.neighbors(target[0])))
    H = G.subgraph(target)
    Jsub = J.subgraph(target)
    del G 
    del J
    thre = 0.08
    epredf = [(u, v) for (u, v, d) in Jsub.edges(data=True) if d["dist"] > thre]
    print(epredf)

    print("Start prediction..")
    print("prediction completed")

    eposi = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] > 0]
    eposiweight = np.log([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] > 0])
    enega = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] < 0]
    enegaweight = np.log(np.negative([d['weight'] for (u, v, d) in H.edges(data=True) if d["weight"] < 0]))
    eneu = [(u, v) for (u, v, d) in H.edges(data=True) if d["weight"] == 0]

    # timestr = time.strftime("%Y%m%d_%H%M%S")
    # figName = "static/fig/" + tGene + timestr + ".png"
    # print(figName, "saved")
    plt.figure(figsize=(20, 20))
    pos = nx.circular_layout(H)
    nx.draw_networkx_nodes(H, pos, node_size=20)
    nx.draw_networkx_labels(H, pos, font_size=30, font_family="sans-serif")
    nx.draw_networkx_edges(H, pos, edgelist = eposi, width=eposiweight, edge_color="r")
    nx.draw_networkx_edges(H, pos, edgelist = enega, width=enegaweight, edge_color="b")
    nx.draw_networkx_edges(H, pos, edgelist = eneu, width=enegaweight, edge_color="y")
    nx.draw_networkx_edges(H, pos, edgelist = epredf, edge_color="g", style = "dashed")

    img = BytesIO()
    plt.savefig(img, format="PNG", dpi=300)
    img.seek(0)

    return send_file(img, mimetype='image/png')
    #return render_template('result.html', fName = figName)

if __name__ == '__main__':
	app.run(host='0.0.0.0', debug=True)
