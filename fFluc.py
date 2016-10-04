#!/usr/bin/env python
'''
fFluc.py will exclude genes with low fluctuating expression with variance cutoff [1], i.e. less variation among expression levels.
fFluc.py will seperately investiate and identify genes that expressed fluctuatedly across all cancer types or in certain cancer types.
fFluc.py will then build the correlation matrix based on spearman's rank correlation, with na cutoff [5].

Finally, fFluc.py will return a set of genes that with correlation > [critera] upon at least [n] drugs.
'''

import os
import os.path
import numpy
import math
import argparse
from scipy.stats.stats import spearmanr as spr
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdrc

__Author__ = 'Yue Wang'
__Date__ = 'Sept 7, 2016'


def b_matrix(expath, dpath):
    '''
    Traverse all expression path and drug information path, and build correlation matrix for all.
    '''
    # get cancer type in dpath
    canType = []
    for parent, dirnames, filenames in os.walk(dpath):
        for name in dirnames:
            canType.append(name)

    for ty in canType:
        # apply cor() function
        cor(expath + ty + '_nc_fluc_expression.txt', dpath + ty + '/')

    print 'Matrix built. Cheers!'

    return


def cor(expression, drugPerfm):
    '''
    cor() will read in two file track, and do correlation test while traversing.
    Input: two file track (fFluc_sigma_/f_divided_expression/[cancertype]_nc_fluc_expression.txt, divided_drug/[cancertype]/)
    Output: a correlation matrix file (col = drug, row = gene)
    '''
    # create output directory
    # re = 'cor_pharm_geno_fdr_bool/' + drugPerfm.split('/')[1] + '/'
    re = 'cor_pharm_geno_fdr_corrected/' + drugPerfm.split('/')[1] + '/'
    isExist = os.path.exists(re)

    if not isExist:
        os.mkdir(re)
    else:
        print 'The directory ' + re + ' already exists. Continue to next step.'

    # factors that will be separatedly tested across cell lines
    factors = ['EC50', 'IC50', 'Amax', 'ActArea']

    def output(factor, col_list, row_dict, fdr_p):
        '''
        Output metrix into output
        '''
        with open(re + 'cor_' + factor + '.txt', 'w') as o:
            # write the headline
            o.write('\t'.join(col_list) + '\n')

            # write correlation coefficients along with gene IDs
            for gene in row_dict.keys():
                out = [gene]
                for i in range(len(row_dict[gene])):
                    if fdr_p[gene][i] != 'false' and row_dict[gene][i] != 'na':
                        out.append(str(row_dict[gene][i]))

                    elif row_dict[gene][i] == 'na':
                        out.append('na')

                    else:
                        # rejected by FDR control
                        out.append('fnan')

                o.write('\t'.join(out) + '\n')

    def fdrCorr(pvector):
        '''
        FDR correction for p-value, only output genes that pass FDR correction, with a cutoff by 0.10 (10%)
        '''
        # pvector without nan
        pNonan = []
        # index for pNonan
        inpN = []

        # pop 'nan' out of pvector for fdr correction
        for i in range(len(pvector)):
            if pvector[i] != 'nan':
                pNonan.append(pvector[i])
                inpN.append(i)
            elif pvector[i] == 'nan':
                pvector[i] = 'false'
            else:
                continue

        # bool vector for results of fdr correction, method = Benjamini/Hochberg for tests are independent??
        bl = list(fdrc(pNonan, alpha=0.1, method='indep', is_sorted=False)[0])

        # merge the bl into pvector when true
        for i in range(len(inpN)):
            # if true
            if bl[i]:
                continue
            # loosen the criteria with p-value < 0.05
            elif pvector[inpN[i]] <= 0.05:
                continue
            else:
                # if false
                pvector[inpN[i]] = 'false'

        return pvector

    # calculate correlation matrix
    for f in factors:
        d = open(drugPerfm + f + '.txt', 'r')

        # initial col factors: drug list
        col = []

        # initial row factors and row panel: gene list + correlation value
        row = {}

        # initial p-value vectors: {gene:[p1, p2, ...]} for FDR correction
        pvalue = {}

        # skip heads
        d.readline()

        for drug in d:
            # every time need re-open the expression file
            e = open(expression, 'r')
            # skip headlines
            e.readline()

            drug = drug.rstrip().split('\t')
            col.append(drug[0])
            nd = drug[1:]
            cd = []

            # pop cell lines with 'NAN' in drug sensitivity metrix
            topop = []
            for i in range(len(nd)):
                if nd[i] == 'NAN' or nd[i] == 'NA':
                    topop.append(i)
                else:
                    continue
            # topop = [i for i in range(len(nd)) if nd[i] == 'NAN' or 'NA']

            for i in range(len(nd)):
                if i not in topop:
                    cd.append(nd[i])
                else:
                    continue

            for gene in e:
                # get current gene expression vector
                gene = gene.rstrip().split('\t')
                ng = [float(ge) for ge in gene[1:]]
                cg = []

                # pop cell lines with 'NAN' in drug sensitivity metrix
                for i in range(len(ng)):
                    if i not in topop:
                        cg.append(ng[i])
                    else:
                        continue

                # Spearman Rank Correlation Test: return(spearman correlation coefficient, 2 tailed p-value)
                # record the score into row dictionary
                # remove lines with too many NAN or NA
                if len(cd) >= 5:
                    # here use cutoff 5; could be further regulated
                    if gene[0] in row.keys():
                        row[gene[0]].append(list(spr(cd, cg))[0])
                        pvalue[gene[0]].append(list(spr(cd, cg))[1])
                    else:
                        row[gene[0]] = [list(spr(cd, cg))[0]]
                        pvalue[gene[0]] = [list(spr(cd, cg))[1]]
                else:
                    if gene[0] in row.keys():
                        row[gene[0]].append('na')
                        pvalue[gene[0]].append('nan')
                    else:
                        row[gene[0]] = ['na']
                        pvalue[gene[0]] = ['nan']

        for gcor in pvalue.keys():
            # apply fdr correction
            cr = fdrCorr(pvalue[gcor])
            pvalue[gcor] = cr

        output(f, col, row, pvalue)

        e.close()
        d.close()

    return


def var_across_type(path, sigma):
    '''
    Gene expression variation across all cancer types, globally calculating the standard deviation.
    Writing to './merge_nc_fluc_expression.txt'
    '''
    o = open('fFluc_sigma_' + str(sigma) + '/merge_nc_fluc_expression.txt', 'w')

    # record genes that pass cut-off
    cgene = []
    with open(path, 'r') as p:
        lname = p.readline().rstrip()
        o.write(lname + '\n')

        for expression in p:
            expression = expression.rstrip().split('\t')

            # t for std calculation
            t = []

            for i in range(len(expression[1:])):
                # logarithm and pseudo count normalization
                expression[i + 1] = str(math.log((float(expression[i + 1]) + 1), 2))
                t.append(float(expression[i + 1]))

            if numpy.std(t, ddof=1) >= sigma:
                # write this gene into output
                o.write('\t'.join(expression) + '\n')
                cgene.append(expression[0])
            else:
                # exclude this gene
                continue

    o.close()

    print 'var_across_type finished.'

    return cgene


def var_in_type(path, sigma):
    '''
    Gene expression varied in one certain type, only calculating the standard deviation in one cancer type.
    Gene expressions have been already normalized by log2(counts + 1), hence the standard deviation here is actually in fold change form
    Return: write to 'f_divided_expression/cancertype_nc_fluc_expression.txt'
    Read expression files by traversing path '/divided_expression/' and filter ones have low standard deviation.
    '''
    track = []
    for parent, dirnames, filenames in os.walk(path):
        for name in filenames:
            track.append(parent + name)

    # create output directory
    d = 'fFluc_sigma_' + str(sigma) + '/f_divided_expression/'
    isExist = os.path.exists(d)
    if not isExist:
        os.mkdir(d)
    else:
        print 'directory ' + d + ' already existed, continue to next step OvO!'

    # record genes that pass cut-off, keys = cancer type
    cgene = {}
    # record genes that pass cut-off with keys = geneID, value = cancertype or cancertypes
    igene = {}

    # traverse all subpath in tracks
    for sp in track:
        # create output file
        cname = sp.split('/')[-1].split('_')[0]

        # initialize the cgene sets
        cgene[cname] = []

        # create output file
        o = open(d + cname + '_nc_fluc_expression.txt', 'w')

        with open(sp, 'r') as ex:
            # get cell line names
            lname = ex.readline().rstrip()
            o.write(lname + '\n')

            # calculate top sigma percent genes cutoff
            # sigma: 0.25
            stdl = []
            for lines in ex:
                lines = lines.rstrip().split('\t')
                # str() to float()
                t = [float(l) for l in lines[1:]]

                # corrected sample standard deviation
                stdl.append(numpy.std(t, ddof=1))

            cutoff = numpy.percentile(stdl, sigma)

            for lines in ex:
                lines = lines.rstrip().split('\t')
                # str() to float()
                t = [float(l) for l in lines[1:]]

                if numpy.std(t, ddof=1) >= cutoff:
                    # write this gene into output
                    o.write('\t'.join(lines) + '\n')
                    # record in cgene sets
                    cgene[cname].append(lines[0])

                    # record in igene sets
                    if lines[0] in igene.keys():
                        igene[lines[0]].append(cname)
                    else:
                        igene[lines[0]] = [cname]

                else:
                    # exclude this gene
                    continue

        o.close()

    print 'var_in_type finished.'

    return cgene, igene


def main():
    '''
    Input data and parameters, output and reform the results
    '''
    # Acquire parameter from outside
    parser = argparse.ArgumentParser()
    parser.add_argument("sigma", help='Cutoff to identify genes expressed differentially across all cancer types or in certain types', type=float)
    args = parser.parse_args()
    sigma = args.sigma

    # apply filtering

    rInType = var_in_type('divided_expression/', sigma)
    rAcrType = var_across_type('noncoding_CCLE_EXP.txt', sigma)

    # get overlap genes among types and in/across types
    # genes of var_in_type
    # cgeneIn = rInType[0]
    # cancer type index for var_in_type genes
    # igeneIn = rInType[1]
    # genes of var_across_type
    # cgeneAcr = rAcrType
    # get overlap genes among cancer types
    # print len([gene for gene in igeneIn.keys() if len(igeneIn[gene]) == 1])
    # print len([gene for gene in igeneIn.keys() if len(igeneIn[gene]) >= 10])
    # print len([gene for gene in igeneIn.keys() if len(igeneIn[gene]) >= 15])
    # print len([gene for gene in cgeneAcr if gene in igeneIn.keys()])

    # build correlation matrix
    b_matrix('fFluc_sigma_' + str(sigma) + '/f_divided_expression/', 'divided_drug/')


if __name__ == '__main__':
    main()
