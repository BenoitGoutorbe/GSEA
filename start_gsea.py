from gsea import *
import time


my_gsea = gsea("leukemia.txt", "pathways.txt", "output.txt", MIN_SIZE = 2)
"""tstart =time.time()
ES1 = my_gsea.get_ES([],False, 1.0)
t1 =time.time()-tstart

tstart =time.time()
ES2 = my_gsea.get_ES_fast([], False, 1.0)
t2 =time.time()-tstart

print("comparaison des 2 methodes :" )

print("scores      (",round(t1,3),"s ) :",ES1[:10])
print("scores fast (",round(t2,3),"s ) :",ES2[:10])
"""
ES  = my_gsea.get_ES_fast([], False, 1.0)
pval, norm_scores = my_gsea.get_pvalue(ES,size_sample=1000, p = 1.0, normalize=True, show=True)
my_gsea.write_output(pval, norm_scores, alpha = 0.05)


