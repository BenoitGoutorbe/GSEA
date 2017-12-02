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

print("Expression levels of", len(genes), "genes for",len(patients), "patients")
data = np.matrix(expr)
print(data.shape)