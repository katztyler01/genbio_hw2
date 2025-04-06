import mygene
import numpy as np

def map_genes(genes: list):
    mg = mygene.MyGeneInfo()
    results = mg.querymany(genes, scopes='symbol', fields=['ensembl.gene', 'entrezgene'], species='human')
    symbol_to_ensembl = {}
    symbol_to_entrez = {}
    for result in results:
        symbol = result.get("query")
        if isinstance(result.get("ensembl"), (list, np.ndarray)):
            ensembls = []
            for pair in result.get("ensembl"):
                if pair.get('gene'):
                    ensembls.append(pair.get('gene'))
            symbol_to_ensembl[symbol] = ensembls
        else:
            if result.get("ensembl") is not None:
                symbol_to_ensembl[symbol] = [result.get("ensembl").get("gene")]
        if result.get("entrezgene"):
            symbol_to_entrez[symbol] = result.get("entrezgene")
    return symbol_to_ensembl, symbol_to_entrez