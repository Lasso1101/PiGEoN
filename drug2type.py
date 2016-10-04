#!/usr/bin/env python
'''
drug2type.py will divide drug_response_detail.csv into cancer-type-dominated drug sensitivity files.
'''
import os
import os.path

__Author__ = 'Yue Wang'
__Date__ = 'Sept 1, 2016'


def drug2file(cancerIndex, indexDrug, drugList):
    '''
    drug2file() is the output function of drug2type.py.
    For each cancer types, drug2file will create a directory, and write EC50, IC50, Amax, ActArea in seperated files, along with cell line names.
    '''
    for cancer in cancerIndex:
        # create directory
        d = 'divided_drug/' + cancer + '/'
        isExist = os.path.exists(d)
        if not isExist:
            os.mkdir(d)
        else:
            print 'The directory ' + d + ' already exists. Continue to next step.'

        print 'Writing output to ' + d

        # define cell line order
        lines = cancerIndex[cancer]

        # write drug responses
        factors = ['EC50', 'IC50', 'Amax', 'ActArea']
        for f in factors:
            # create file
            with open(d + f + '.txt', 'w') as o:
                # write title
                o.write('\t'.join(lines) + '\n')

                # get drug names drug response values
                for drug in drugList:
                    # get drug name
                    l = [drug]

                    # get drug response values
                    for cell in lines:
                        # if this drug has been tested on this cell line
                        if drug in indexDrug[cell].keys():
                            l.append(str(indexDrug[cell][drug][f]))

                        else:
                            # if this drug hasn't been tested on this cell line, return NAN
                            l.append('NAN')

                    # if all NAN, abandon this drug and continue
                    if l[1:] != ['NAN' for i in range(len(lines) - 1)]:
                        o.write('\t'.join(l) + '\n')
                    else:
                        continue

    return


def drugIndex(file, indexID):
    '''
    drugIndex will build index for drug responses and cell lines.
    Return a dictionary: keys = Cell line names (trans to Atlas format), values = [drugs{EC50, IC50, Amax, AArea}].
    '''
    indexDrug = {}
    with open(file, 'r') as f:
        # EC50, IC50, Amax, ActArea
        factors = f.readline().rstrip().split('\t')[2:]
        factors[0] = 'EC50'
        factors[1] = 'IC50'

        # read each line and build the dictionary
        for lines in f:
            lines = lines.rstrip().split('\t')

            # if cell line in CCLE_Atlas
            if lines[0] in indexID.keys():
                # trans cell line name into Atlas format
                aName = transID(lines[0], indexID)

                # initial drug dict
                drug = {}
                # trans data type from str to float; if NA, remain str
                for i in range(len(lines[2:])):
                    if lines[i + 2] != 'NA':
                        drug[factors[i]] = float(lines[i + 2])

                    else:
                        drug[factors[i]] = lines[i + 2]

                # Pop drug dict into indexDrug[aName]
                if aName in indexDrug.keys():
                    indexDrug[aName][lines[1]] = drug

                else:
                    indexDrug[aName] = {}
                    indexDrug[aName][lines[1]] = drug

            else:
                continue

    return indexDrug


def transID(opDrug, indexID):
    '''
    Transform id for cell lines into id for drug sensitivities.
    Return: opDrug.
    '''
    return indexID[opDrug]


def binID(opAtlas, opDrug):
    '''
    Build index for id transformation.
    '''
    indexID = {}
    for i in range(len(opDrug)):
        indexID[opDrug[i]] = opAtlas[i]

    return indexID


def readCancer(direct):
    '''
    Read expression files in separated cancer types.
    Return a dictionary: Keys = cancer type, values = cell line names
    '''
    # build track
    track = []
    for parent, dirnames, filenames in os.walk(direct):
        for name in filenames:
            if name[-4:] == '.txt':
                track.append(parent + name)

    # build index
    cancerIndex = {}
    for t in track:
        b = open(t, 'r')
        cancerIndex[t.split('/')[-1].split('.')[0].split('_')[0]] = b.readline().rstrip().split('\t')
        b.close()

    print 'Build cancer index finished.'
    return cancerIndex


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

    # ID trans index
    indexID = binID(opAtlas, opDrug)

    # divided_expression directory
    d = 'divided_expression/'
    cancerIndex = readCancer(d)

    # build drug index
    indexDrug = drugIndex('drug_response_detail.csv', indexID)
    # print indexDrug['KYSE-70']['Topotecan']['EC50']

    # get all drug names
    with open('drug.txt', 'r') as g:
        drugList = []
        for lines in g:
            lines = lines.rstrip()
            drugList.append(lines)

    # write output
    drug2file(cancerIndex, indexDrug, drugList)


if __name__ == '__main__':
    main()
