#!/usr/bin/env python
'''
ex2type.py will divide noncoding_CCLE_EXP.txt into cancer-type-dominated expression files.
'''
import os
import os.path
import math

__Author__ = 'Yue Wang'
__Date__ = 'Aug 31, 2016'


def ex2file(allType, gorder, exIn):
    '''
    ex2file() is the output of ex2type.py, writing expressions of noncoding genes in cell lines in their corresponding cancer type files.
    '''
    for ct in allType.keys():
        o = open('divided_expression/' + ct + '_nc_expression.txt', 'w')
        title = allType[ct]
        o.write('\t'.join(title) + '\n')

        for i in range(len(gorder)):
            seq = [gorder[i]]

            for lines in title:
                seq.append(str(exIn[lines][i]))

            o.write('\t'.join(seq) + '\n')

        o.close()

    print 'ex2type.py process finished.'

    return


def exIndex(expression):
    '''
    exIndex() will build index of expression file, meanwhile, maintaining the order of gene IDs.
        Keys: cell lines, Value: a vector(list) contain log2(read counts + 1), in IDs' order
    '''
    with open(expression, 'r') as e:
        # cell line names
        title = e.readline().rstrip().split('\t')
        exIn = {}
        for k in title:
            exIn[k] = []

        # gene ID series
        gorder = []

        # build the index
        for line in e:
            line = line.rstrip().split('\t')
            gorder.append(line[0])

            for i in range(len(title)):
                # log2 normalization
                exIn[title[i]].append(str(math.log((int(line[i + 1]) + 1), 2)))

    print 'Build expression index finished.'

    return gorder, exIn


def catg(op, index):
    '''
    categorize cell lines in op according to their cancer type in index
    '''
    allCat = {}
    for cell in op:
        if index[cell] in allCat.keys():
            allCat[index[cell]].append(cell)
        else:
            allCat[index[cell]] = [cell]

    print 'Build cancer type index finished.'

    return allCat


def bl(path):
    '''
    bl() will build a dictionary for cancer cell lines according to their cancer type
    Dict structure: keys = cell line ID, values = cancer type
    Return: dict(cell line index)
    '''
    # build track
    track = []
    for parent, dirnames, filenames in os.walk(path):
        for name in filenames:
            if name[-4:] == '.txt':
                track.append(parent + name)

    # build index
    cellIndex = {}
    for t in track:
        b = open(t, 'r')
        for lines in b:
            cellIndex[lines.rstrip()] = t.split('/')[1].split('.')[0]

        b.close()

    print 'Build cell line index finished.'
    return cellIndex


def opCL(file):
    '''
    opCL() will read overlap cell lines recorded in ppp.txt or modified_ppp.txt
    '''
    with open(file, 'r') as f:
        op = []
        for lines in f:
            op.append(lines.rstrip())

    print 'Read cell lines finished.'

    return op


def main():
    # for matching cell lines in drug sensitivity
    opDrug = opCL('ppp.txt')
    # for matching cell lines in Atlas
    opAtlas = opCL('modified_ppp.txt')

    # build cell-type index
    cellIndex = bl('blocks_8_29/')

    # categorize
    allType = catg(opAtlas, cellIndex)

    # build index of expression file
    r = exIndex('noncoding_CCLE_EXP.txt')
    gorder = r[0]
    exIn = r[1]

    # process the categorization
    ex2file(allType, gorder, exIn)


if __name__ == '__main__':
    main()
