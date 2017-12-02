import numpy as np
import matplotlib.pyplot as plt

# reading the expression file
file = open("leukemia.txt", "r")
content = file.readlines()
file.close()

patients = [] # patients category (ALL or AML)
genes = [] # name of the genes (headers)
expr = [] # table with expression levels
patients = content[0].split()[1:]
for line in content[1:] :
    genes.append(line.split()[0])
    expr.append([int(value) for value in line.split()[1:]])
NB_genes = len(genes)
NB_patients = len(patients)
print("Expression levels of", NB_genes, "genes for",NB_patients, "patients")

data = np.matrix(expr)

# reading the patwhays file
file = open("pathways.txt", "r")
content = file.readlines()
file.close()

MIN_SIZE = 15
pathways = [] # name of the pathways
genes_sets = [] # name of the genes implicated in each pathway
# NOTE : we keep only the genes sets with at least MIN_SIZE (15) genes for which we have expression data
for line in content :
    genes_implicated = [g for g in line.split()[1:] if g in genes]
    if len(genes_implicated) >= MIN_SIZE :
        pathways.append(line.split()[0])
        genes_sets.append(genes_implicated)
NB_sets = len(pathways)
print(NB_sets,"sets of genes used of the analysis")

# Calculate the Enrichment score
size = [len(p)for p in genes_sets] # number of genes implicated in each pathway
ALL_index = [i for i in range(NB_patients) if patients[i] == 'ALL']
AML_index = [i for i in range(NB_patients) if patients[i] == 'AML']
dif_expr = []

for i in range(NB_genes):
    expr_ALL = np.mean([data[i,j] for j in ALL_index])
    expr_AML = np.mean([data[i,j] for j in AML_index])
    dif_expr.append(expr_ALL-expr_AML) # negative values indicate a gene expressed (in average) more by AML patients

# sorting the genes by correlation with phenotye
index_sort = np.argsort(dif_expr)
genes = [genes[i] for i in index_sort]
data = np.matrix([data[i,].tolist()[0] for i in index_sort])
dif_expr = [dif_expr[i] for i in index_sort]
plt.plot(dif_expr)
plt.show()

ES = []
p = 1
for set in range(NB_sets) :
    N_H = size[set]
    N_R = np.sum([ pow(abs(dif_expr[i]),p) for i in range(NB_genes) if genes[i] in genes_sets[set] ])
    eq_0=[0]
    for i in range(len(genes)):
        if genes[i] in genes_sets[set]:
            eq_0.append(eq_0[-1] + pow(abs(dif_expr[i]), p) / N_R)
        else:
            eq_0.append(eq_0[-1]-1/(NB_genes-N_H))
    ES.append(np.max(np.abs(eq_0)))
plt.hist(ES, 20)
plt.show()
