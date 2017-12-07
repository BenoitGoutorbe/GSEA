# GSEA
Here is an implementation of the Gene Sets Enrichment Analysis (Subramanian et al.)  applied to Golub data set (transcriptomic profile of patiens with 2 kinds of leukemia : AML and ALL)

[IMPLEMENTED WITH PYTHON 3.6.3 and NUMPY 1.13.3]

Using from terminal : 
> python gsea.py leukemia.txt pathways.txt output_GSEA

Using as a library : 
(shown in start_gsea.py)
 
 Most important parameters to be set :
* MIN_SIZE (in the __init__) : minimum size of considered pathway (15 in the reference paper) 
* size_sample (in get_pvalue) : number of simulation de generate the null distribution
* show (in get_pvalue) : Create (or not) an figure saved as 'output_GSEA.png'

OUTPUT :
.txt : [p-value	NES	OverExpressingCategory	Pathway] for each set of the analysis, sorted by p-value
.png : figure of the observed NES and the null distribution
