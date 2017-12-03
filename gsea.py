import numpy as np
import matplotlib.pyplot as plt

class gsea:
    def __init__(self, expression_file, pathway_file, output_file):
        self.output_file = output_file
        # reading the expression file
        file = open(expression_file, "r")
        content = file.readlines()
        file.close()

        self.patients = [] # patients category (ALL or AML)
        self.genes = [] # name of the genes (headers)
        expr = [] # table with expression levels
        self.patients = content[0].split()[1:]
        for line in content[1:] :
            self.genes.append(line.split()[0])
            expr.append([int(value) for value in line.split()[1:]])
        self.NB_genes = len(self.genes)
        self.NB_patients = len(self.patients)
        print("Expression levels of", self.NB_genes, "genes for",self.NB_patients, "patients")

        self.expr_matrix = np.matrix(expr)

        # reading the patwhays file
        file = open(pathway_file, "r")
        content = file.readlines()
        file.close()

        MIN_SIZE = 15
        self.pathways = [] # name of the pathways
        self.genes_sets = [] # name of the genes implicated in each pathway
        # NOTE : we keep only the genes sets with at least MIN_SIZE (15) genes for which we have expression data
        for line in content :
            genes_implicated = [g for g in line.split()[1:] if g in self.genes]
            if len(genes_implicated) >= MIN_SIZE :
                self.pathways.append(line.split()[0])
                self.genes_sets.append(genes_implicated)
        self.NB_sets = len(self.pathways)
        print(self.NB_sets,"sets of genes used of the analysis")

    def get_ES(self, patients = [], show = False,p = 1):
        if len(patients) == 0 :
            patients = self.patients

        ALL_index = [i for i in range(self.NB_patients) if patients[i] == 'ALL']
        AML_index = [i for i in range(self.NB_patients) if patients[i] == 'AML']
        dif_expr = []

        for i in range(self.NB_genes):
            expr_ALL = np.mean([self.expr_matrix[i,j] for j in ALL_index])
            expr_AML = np.mean([self.expr_matrix[i,j] for j in AML_index])
            dif_expr.append(expr_ALL-expr_AML) # negative values indicate a gene expressed (in average) more by AML patients

        self.correlation = list(dif_expr)

        # sorting the genes by correlation with phenotye
        index_sort = np.argsort(dif_expr)
        genes_sorted = [self.genes[i] for i in index_sort]
        dif_expr = [dif_expr[i] for i in index_sort]

        ES = []
        for set in self.genes_sets :
            N_H = len(set)
            N_R = np.sum([ pow(abs(dif_expr[i]),p) for i in range(self.NB_genes) if genes_sorted[i] in set ])
            dist2_0 =[0] #distance to 0
            for i in range(self.NB_genes):
                if genes_sorted[i] in set :
                    dist2_0.append(dist2_0[-1] + pow(abs(dif_expr[i]), p) / N_R)
                else:
                    dist2_0.append(dist2_0[-1]-1/(self.NB_genes-N_H))
            ES.append(np.max(np.abs(dist2_0)))
        if show :
            plt.subplot(121)
            plt.plot(dif_expr)
            plt.title("Genes expression distribution")
            plt.subplot(122)
            plt.hist(ES, 20)
            plt.title("Enrichment scores")
            plt.show()
        return ES


    def get_random_distrib (self, size, show = True,p = 1, normalize = 1):
        random_patients = np.array(self.patients)
        ES0 = []
        print('Generating the random distribution...')
        print('0%[' + 50 * ' ' + ']100%')
        print('   ', end='')
        cur = 0
        for i in range(size):
            np.random.shuffle(random_patients)
            random_scores = self.get_ES(random_patients, False, p)
            ES0.append(random_scores)
            if int(50.0 * i / size) != cur:
                print((int(50.0 * i / size) - cur) * '=', end='', flush = True)
                cur = int(50.0 * i / size)
        print((50 - cur) * '=')
        ES_means = np.mean(ES0, axis=0)
        distrib = []
        for line in ES0 :
            if normalize :
                distrib += np.divide(line, ES_means).tolist()
            else :
                distrib += line
        if show :
            plt.hist(distrib, 20)
            plt.title("Random distribution scores for " + str(size*self.NB_sets) + " sets of genes")
            plt.show()
        return (distrib, ES_means)

    def get_pvalue (self, scores, size_sample, normalize = True) :
        (random_distrib, ES_means) = self.get_random_distrib(size_sample, True, 1, normalize)
        if normalize :
            Norm_ES = np.divide(scores, ES_means).tolist()
        else :
            Norm_ES = scores
        pval = []
        for score in Norm_ES :
            pval.append(np.sum(np.greater(random_distrib, score))/size_sample)
        return (pval, Norm_ES)

    def write_output(self, p_values,  norm_scores, alpha):
        nb_output_sets = np.sum(np.less(p_values, alpha))
        hits = np.argsort(p_values)[0:nb_output_sets]
        if nb_output_sets == 0 : print('No pathways such as p_value<alpha ')
        for set in hits :
            print(self.pathways[set]+'      '+str(norm_scores[set])+'     '
                  +str(p_values[set]))



