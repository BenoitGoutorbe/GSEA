from gsea import *

my_gsea = gsea("leukemia.txt", "pathways.txt", "output.txt")
scores = my_gsea.get_ES([], False, 1.0)
(pval, norm_scores) = my_gsea.get_pvalue(scores,100, True)
my_gsea.write_output(pval, norm_scores, 0.05)


