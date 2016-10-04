#!/usr/bin/env python
'''
getOpCL.py will get overlap cellines from cell line files, together with their cancer type, and compounds that have been treated on these lines.
'''

__Author__ = 'Yue Wang'
__Date__ = 'Aug 30, 2016'

import os
import os.path


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

    return cellIndex


def opCL(file):
    '''
    opCL() will read overlap cell lines recorded in ppp.txt or modified_ppp.txt
    '''
    with open(file, 'r') as f:
        op = []
        for lines in f:
            op.append(lines.rstrip())

    return op


def main():
    # read overlap cell lines in ppp.txt
    opDrug = opCL('ppp.txt')
    opAtlas = opCL('modified_ppp.txt')

    # build cell-type index
    cellIndex = bl('blocks_8_29/')

    # categorize
    allType = catg(opAtlas, cellIndex)

    with open('CLstat.csv', 'w') as o:
        for ct in allType.keys():
            o.write(ct + '\t' + str(len(allType[ct])) + '\n')


if __name__ == '__main__':
    main()
