import numpy as np
import matplotlib.pyplot as plt
import decimal

class gsea:
    def __init__(self, expression_file, pathway_file, output_file,MIN_SIZE = 2):
        self.output_file = output_file
        # reading the expression file
        file = open(expression_file, "r")
        content = file.readlines()
        file.close()

        self.patients = []  # patients category (ALL or AML)
        self.genes = []  # name of the genes (headers)
        expr = []  # table with expression levels
        self.patients = content[0].split()[1:]
        for line in content[1:]:
            self.genes.append(line.split()[0])
            expr.append([int(value) for value in line.split()[1:]])

        # reading the patwhays file
        file = open(pathway_file, "r")
        content = file.readlines()
        file.close()
        print("Input data :" , len(self.patients) , "patients," ,
              len(self.genes) , "genes, " , len(content) , "pathways")

        self.pathways = []  # name of the pathways
        self.genes_sets = []  # name of the genes implicated in each pathway
        genes_presence = [0] * len(self.genes)
        # NOTE : we keep only the genes sets with at least MIN_SIZE genes for which we have expression data
        for line in content:
            genes_implicated = [g for g in line.split()[1:] if g in self.genes]
            if len(genes_implicated) >= MIN_SIZE:
                self.pathways.append(line.split()[0])
                self.genes_sets.append(genes_implicated)
                for g in genes_implicated :
                    genes_presence[self.genes.index(g)] += 1
        # clearing the unused genes
        a = 0
        for i in range(len(self.genes)) :
            if genes_presence[i] == 0 :
                del self.genes[i-a]
                del expr[i-a]
                a = a+1

        self.expr_matrix = np.matrix(expr)

        self.index_genes_implicated = []
        for set in self.genes_sets :
            self.index_genes_implicated.append([self.genes.index(g) for g in set])

        self.NB_patients = len(self.patients)
        self.NB_genes = len(self.genes)
        self.NB_sets = len(self.pathways)
        print("Reduced data (only pathways with at least",MIN_SIZE,"genes and only genes belonging to a selected pathway)")
        print(self.NB_patients , "patients," , self.NB_genes , "genes, " , self.NB_sets , "pathways")

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
            N_R = np.sum([ pow(abs(dif_expr[i]),p) for i in range(self.NB_genes) if genes_sorted[i] in set])
            dist2_0 =[0] #distance to 0
            for i in range(self.NB_genes):
                if genes_sorted[i] in set :
                    dist2_0.append(dist2_0[-1] + pow(abs(dif_expr[i]), p) / N_R)
                else:
                    dist2_0.append(dist2_0[-1] - 1.0/(self.NB_genes-N_H) )
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

    def get_ES_fast(self, patients = [], show = False,p = 1):
        if len(patients) == 0 :
            remind = True
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
        index_sort = np.argsort(dif_expr).tolist()
        positions = [index_sort.index(i) for i in range(self.NB_genes)]
        dif_expr_sorted = [dif_expr[i] for i in index_sort]
        ES = []

        for set in range(self.NB_sets) :
            N_H = len(self.genes_sets[set])
            N_R = np.sum(pow(abs(dif_expr[i]),p) for i in self.index_genes_implicated[set])
            penalty = -1.0/(self.NB_genes-N_H)
            positions_genes_implicated = np.sort([positions[i] for i in self.index_genes_implicated[set]]).tolist()
            path = [0]
            if positions_genes_implicated[0] != 0 :
                path.append((positions_genes_implicated[0])* penalty)
                path.append(path[-1] + pow(abs(dif_expr_sorted[positions_genes_implicated[0]]), p) / N_R)
            else :
                path = [pow(abs(dif_expr[positions_genes_implicated[0]]), p) /N_R]
            for i in range(1,len(positions_genes_implicated)) :
                path.append(path[-1]+(positions_genes_implicated[i]-positions_genes_implicated[i-1]-1)*penalty)
                path.append(path[-1] + pow(abs(dif_expr_sorted[positions_genes_implicated[i]]), p) /N_R)
            path.append(path[-1]+float((self.NB_genes-positions_genes_implicated[-1]-1)*penalty))
            max_deviation_position = np.argmax(np.abs(path))
            ES.append(path[max_deviation_position])


        if show :
            plt.subplot(121)
            plt.plot(np.sort(dif_expr))
            plt.title("Genes expression distribution")
            plt.subplot(122)
            plt.hist(ES, 20)
            plt.title("Enrichment scores")
            plt.show()
        return ES


    def get_random_distrib (self, size,p = 1, normalize = True):
        random_patients = np.array(self.patients)
        ES0 = []
        print('Generating the random distribution...')
        print('0%[' + 50 * ' ' + ']100%')
        print('   ', end='')
        cur = 0
        for i in range(size):
            np.random.shuffle(random_patients)
            random_scores = self.get_ES_fast(random_patients, False, p)
            ES0.append(random_scores)
            if int(50.0 * i / size) != cur:
                print((int(50.0 * i / size) - cur) * '=', end='', flush = True)
                cur = int(50.0 * i / size)
        print((50 - cur) * '=')
        if normalize : print('Normalization...', end = '', flush = True)
        ES_pos_means = []
        ES_neg_means = []
        for j in range(self.NB_sets):
            pos = []
            neg = []
            for i in range(size) :
                if  ES0[i][j] > 0 :
                    pos.append(ES0[i][j])
                else :
                    neg.append(ES0[i][j])
            if len(pos) == 0 :
                pos.append(1.0)
            if len(neg) == 0 :
                neg.append(-1.0)
            ES_pos_means.append(np.mean(pos))
            ES_neg_means.append(-np.mean(neg))

        distrib = []
        for ES in ES0 :
            if normalize :
                for i in range(self.NB_sets) :
                    if ES[i] >0 :
                        distrib.append(ES[i]/ES_pos_means[i])
                    else :
                        distrib.append(ES[i] / ES_neg_means[i])
            else :
                distrib += ES
        return (distrib, ES_pos_means, ES_neg_means)

    def get_pvalue (self, ES, size_sample,p=1, normalize = True, show = False) :
        random_distrib, ES_pos_means,ES_neg_means = self.get_random_distrib(size_sample, p, normalize)
        if normalize :
            NES= []
            for i in range(self.NB_sets) :
                if ES[i] > 0 :
                    NES.append( ES[i] / ES_pos_means[i])
                else :
                    NES.append(ES[i] / ES_neg_means[i])
        else :
            NES = ES
        pval = []
        for i in range(self.NB_sets) :
            if NES[i] > 0 :
                pval.append(np.sum(np.greater(random_distrib,NES[i])) / len(random_distrib))
            else :
                pval.append(np.sum(np.less(random_distrib,NES[i])) / len(random_distrib))
        if normalize: print('  DONE')
        if show :
            xmin, xmax = np.min(random_distrib)-0.2, np.max(random_distrib)+0.2
            x = np.arange(xmin,xmax,0.1).tolist()
            y = [ sum([val >= xi and val < xi+0.1 for val in random_distrib]) for xi in x]
            y = (np.array(y) / size_sample).tolist()
            plt.plot(np.array(x)+0.05,y, 'r--')
            plt.title("Random distribution scores for " + str(size_sample*self.NB_sets) + " sets of genes")
            plt.hist(NES, bins=x)
            plt.savefig('distribution.png')
        return (pval, NES)

    def write_output(self, p_values,  norm_scores, alpha = 0.05):
        file = open(self.output_file, 'w')
        file.write("Pathway p-value NES Leukemia Correlation \n")
        nb_hits = max(10, len(np.argwhere(np.array(p_values) < alpha)))
        hits = np.argsort(p_values).tolist()[:nb_hits]
        for set in hits :
            mean_correlation = np.mean([self.correlation[g] for g in self.index_genes_implicated[set]])
            if mean_correlation < 0 :
                leu = "\tAML\t"
            else :
                leu = "\tALL\t"
            print(self.pathways[set]+"\t"+str(p_values[set])+"\t"+str(norm_scores[set])+leu+ str(mean_correlation))
            file.write(self.pathways[set]+"\t"+str(p_values[set])+"\t"+str(norm_scores[set])+leu+ str(mean_correlation)+" \n")
        file.close()
        if p_values[hits[0]] >= alpha: print('NO SETS WITH p-value < ', alpha, '  ->  WE SAVED THE 10 BEST HITS')