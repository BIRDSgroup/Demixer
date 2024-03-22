from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer;
from sklearn.decomposition import NMF;
from sklearn.preprocessing import normalize;
import pandas as pd;
import pickle

def run_NMF(doc_term,K,strain,run=0):
    transformer = TfidfTransformer(smooth_idf=False);
    doc_tfidf = transformer.fit_transform(doc_term);

    xtfidf_norm = normalize(doc_tfidf, norm='l1', axis=1)
    
    model = NMF(n_components=K, init='nndsvd',max_iter=500)
    W=model.fit_transform(xtfidf_norm)
    #W=model.fit_transform(doc_term,setargs.n_m_z,setargs.n_z_t)
    #W=model.fit_transform(xtfidf_norm,W=custom_W,H=custom_H)
    H = model.components_
    
    est_theta=W/W.sum(axis=1)[:,None]

    
    pd.DataFrame(est_theta).to_csv(strain+str(run)+'_NMF'+'.csv',index=False)
  
    return W,H;
