## pip install biopython
## pip install pandas

from Bio import Entrez
## 전체 코드의 핵심 library, Entrez 데이터베이스를 crawl 할 수 있는 api 
## https://biopython.org/docs/1.75/api/Bio.Entrez.html 참고 
import pandas as pd
from timeit import default_timer as timer
from urllib.error import HTTPError
import time

## Gene list 불러오기 
gene_data = pd.read_csv('static/data/gene_with_protein_product.txt', sep = '\t')
genes = list(gene_data['symbol'])

def getRelations(ids, src):
    ## 1. id 에 해당하는 article들을 불러와서 
    ## 2. 두개의 유전자(src와 다른 하나)가 포함된 문장들을 저장하기 위한 함수 
    res = []
    try:
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                              id=ids)
        results = Entrez.read(handle)
        ## id에 해당하는 article들을 불러옴, Dictionary 형태로 저장됨 
    except:
        return res
    
    for paper in results['PubmedArticle']:
        try:
            curTexts = paper['MedlineCitation']['Article']['Abstract']['AbstractText']
            ## Abstracttext만 긁어오기 
            PMID = str(paper['MedlineCitation']['PMID'])
            ## PMID도 함께 저장 ( 추가됨 )
            sentences = []
            abs = ""
            for para in curTexts:
                ## Abstract가 각 문장별로 저장되어있는경우와 아닌 경우가 있음( JAMA style, Nature style )
                ## 문장단위로 분할하여 저장 
                sentences.extend(para.strip(' .').split('.'))
                abs += para.strip(' .').split('.')[0] + ". "
            for cur in sentences:
                words = cur.strip().split(' ')
                ## 문장을 다시 단어단위로 나누어서 
                if src in words:
                    for tar in genes:
                        ## 두개의 Gene이 한 문장에 들어있는 경우만 저장하기 
                        if src == tar:
                            continue
                        if tar in words:
                            tmp = [src, tar, cur, PMID, abs]
                            res.append(tmp)
        except:
            continue
    
    return res

if __name__ == "__main__":
    """
        1. 특정 유전자로 검색되는 article의 id를 받아서 
        2. getRelations 함수를 통해 해당 gene과 다른 gene이 동시에 등장하는 문장을 검색하고 
        3. 이를 액셀에 저장 ( entrez api는 속도가 빠른 편이 아니기 때문에, 중간중간 다른 이름으로 저장하였음 )
    """
    setAll = set()
    ### crawl 시간을 줄이기 위한 방법. 이전에 탐색하지 않은 article만 분석하기 위해 set 사용 
    resSen = []
    ### 결과 저장을 위한 변수, 진짜 결과 
    meta = []
    ### 결과 저장을 위한 변수 222, 검색한 article 수, 유전자 수 등 meta 변수들을 저장함 
    totalQ = 0
    ### 결과 저장을 위한 변수 333, 현재까지 검색한 article 수  
    idx = 0
    ### 중간 저장을 위해 idx를 하나씩 더해가고, 500 단위로 중간 저장 

    start = timer()


    for src in genes:
        idx += 1
        try:
            handle = Entrez.esearch(db='pubmed', retmode='xml', restatrt = idx, retmax = 100000, term=src)
            results = Entrez.read(handle)
            ### src gene을 가진 article들의 id를 받아옴 
        except HTTPError:
            ### 간혹 너무 많은 query를 날려서 entrez 서버에서 505에러를 뱉어내는 경우가 있음, 이 경우 5초 뒤에 다시 쿼리하면 됨 
            print('ha....')
            time.sleep(5)
            handle = Entrez.esearch(db='pubmed', retmode='xml', restatrt = idx, retmax = 100000, term=src)
            results = Entrez.read(handle)

        IdSet = set(results['IdList'])
        
        curSearch = IdSet - setAll
        ### 기존에 쿼리하지 않은 article만 query하기 위해 차집합 사용 
        
        tmp = getRelations(curSearch, src)
        ### getRelations 함수를 통해 결과 받아옴 
        resSen.extend(tmp)
        
        ### 아래는 결과 저장을 위한 코드임 
        setAll.update(IdSet)
        totalQ += len(IdSet)
        end = timer()
        tmpMeta = [src, len(IdSet), totalQ, len(setAll), len(tmp), len(resSen), end - start]
        meta.append(tmpMeta)
        print(tmpMeta)
        
        # if idx % 500 == 0:
        #     tt = pd.DataFrame(meta)
        #     tt.to_excel('drive/My Drive/Data/meta_' + str(idx) + '.xlsx', header = False, index = False)
        #     ttt = pd.DataFrame(resSen)
        #     ttt.to_excel('drive/My Drive/Data/res_' + str(idx) + '.xlsx', header = False, index = False)    
    tt = pd.DataFrame(meta)
    tt.to_excel('drive/My Drive/Data/meta_' + str(idx) + '.xlsx', header = False, index = False)
    ttt = pd.DataFrame(resSen)
    ttt.to_excel('drive/My Drive/Data/res_' + str(idx) + '.xlsx', header = False, index = False)   