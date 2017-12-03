from gsea import *

my_gsea = gsea("leukemia.txt", "pathways.txt", "output")
scores = my_gsea.get_ES([], False, 1.0)
H0 = my_gsea.get_random_distrib(50, True, 1.0)
pval = my_gsea.get_pvalue(scores, H0)
my_gsea.write_output(pval,0.02)


