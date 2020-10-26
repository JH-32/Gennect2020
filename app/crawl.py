## pip install biopython
## pip install pandas

from Bio import Entrez
import pandas as pd
from timeit import default_timer as timer
from urllib.error import HTTPError
import time

gene_data = pd.read_csv('static/data/gene_with_protein_product.txt', sep = '\t')
genes = list(gene_data['symbol'])

def getRelations(ids, src):
    res = []
    try:
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                              id=ids)
        results = Entrez.read(handle)
    except:
        return res
    
    for paper in results['PubmedArticle']:
        try:
            curTexts = paper['MedlineCitation']['Article']['Abstract']['AbstractText']
            PMID = str(paper['MedlineCitation']['PMID'])
            sentences = []
            abs = ""
            for para in curTexts:
                sentences.extend(para.strip(' .').split('.'))
                abs += para.strip(' .').split('.')[0] + ". "
            for cur in sentences:
                words = cur.strip().split(' ')
                if src in words:
                    for tar in genes:
                        if src == tar:
                            continue
                        if tar in words:
                            tmp = [src, tar, cur, PMID, abs]
                            res.append(tmp)
        except:
            continue
    
    return res

if __name__ == "__main__":
    setAll = set()
    resSen = []
    meta = []
    totalQ = 0
    idx = 0

    start = timer()


    for src in genes:
        idx += 1
        try:
            handle = Entrez.esearch(db='pubmed', retmode='xml', restatrt = idx, retmax = 100000, term=src)
            results = Entrez.read(handle)
        except HTTPError:
            print('ha....')
            time.sleep(5)
            handle = Entrez.esearch(db='pubmed', retmode='xml', restatrt = idx, retmax = 100000, term=src)
            results = Entrez.read(handle)

        IdSet = set(results['IdList'])
        
        curSearch = IdSet - setAll
        
        tmp = getRelations(curSearch, src)
        resSen.extend(tmp)
        
        setAll.update(IdSet)
        totalQ += len(IdSet)
        end = timer()
        tmpMeta = [src, len(IdSet), totalQ, len(setAll), len(tmp), len(resSen), end - start]
        meta.append(tmpMeta)
        print(tmpMeta)
        
        if idx % 500 == 0:
            tt = pd.DataFrame(meta)
            tt.to_excel('drive/My Drive/Data/meta_' + str(idx) + '.xlsx', header = False, index = False)
            ttt = pd.DataFrame(resSen)
            ttt.to_excel('drive/My Drive/Data/res_' + str(idx) + '.xlsx', header = False, index = False)    