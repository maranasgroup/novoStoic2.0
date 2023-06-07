import os
from optstoicpy.core.reaction import *
#import csv

# def parse_gams_input_remove_comment_and_quotes(data_filepath, data_filename, begin=1, end=-1):
#     f = open(os.path.join(data_filepath, data_filename), 'rU')
#     #remove quotes if begin and end is specified
#     # don't read in lines with gams comment "*"
#     data = []
#     # for line in f.readlines():
#     #     if (not line.startswith("*")) and (len(line.strip()) > 0):
#     #         try:
#     #             data.append(line.strip().split()[0][begin:end])
#     #         except:
#     #             print line
#     #             continue

#     data = [line.strip().split()[0][begin:end] for line in f.readlines()
#         if (not line.startswith("*") and len(line.strip())> 0)]
#     return data


def convert_set_to_list(filename):
    """
    replacing function parse_gams_input_remove_comment_and_quotes

    generate a list from gams set data file
    """
    f = open(filename, 'rU')
    # remove quotes if begin and end is specified
    # don't read in lines with gams comment "*"
    data = []
    for line in f.readlines():
        if (not line.startswith("*")) and (len(line.strip()) > 0):
            if line.startswith("/"):
                continue
            data.append(line.replace("\'", "").split()[0])

    return data


def convert_parameter_table_to_dict(filename, Sdict=None):
    """
    old function name : parse_gams_input_Smatrix
    generate a 2D dictionary based on Gams parameter table

    output is a python dictionary:
    e.g. 'i'.'j' value  =>  S[i][j] = value


    """
    f = open(filename, 'rU')
    # remove quotes and don't read in lines with gams comment "*"
    if Sdict is None:
        Sdict = {}
    for line in f.readlines():
        if (not line.startswith("*") and len(line.strip()) > 0):
            if line.startswith("/"):
                continue
            entries = line.strip().split()
            met, rxn = entries[0][1:-1].split("'.'")
            if met not in Sdict:
                Sdict[met] = {}
            Sdict[met][rxn] = float(entries[1])

    return Sdict


def convert_parameter_list_to_dict(filename, datadict=None):
    """
    old function name: parse_gams_input_rxntype
    generate a dictionary from gams parameter input
    args:
        filename: filepath of data in gams parameter format
                    'R01425'    0\n,
                    'R01426'    2\n,
        datadict: Existing dictionary (optional) Data is append to existing dictionary if provided.
    """
    f = open(filename, 'rU')

    if datadict is None:
        datadict = {}

    # remove quotes and don't read in lines with gams comment "*"
    for line in f.readlines():
        if (not line.startswith("*") and len(line.strip()) > 0):
            if line.startswith("/"):
                continue
            # remove single quotes
            entries = line.replace("'", "").split()
            datadict[entries[0]] = float(entries[1])
    return datadict


def write_list_to_file(reaction_list, outputfilename, quotes=False):
    f = open(outputfilename, 'w')
    if quotes:
        for rxn in reaction_list:
            f.write("'%s'\n" % rxn)
    else:
        for rxn in reaction_list:
            f.write("%s\n" % rxn)
    f.close()
    return 1


def write_dict_to_file(paramdict, outputfilename, quotes=False):
    f = open(outputfilename, 'w')
    if quotes:
        for k, v in paramdict.items():
            f.write("'%s' %s\n" % (k, v))
    else:
        for k, v in paramdict.items():
            f.write("%s %s\n" % (k, v))
    f.close()
    return 1


def write_nested_dict_to_file(data, outputfilename, orient='second'):
    """
    Args:
    data = nested dictionary (2D)
    {
        "FACOAE120": {
            "C00001": -1.0,
            "C00010": 1.0,
            "C01832": -1.0,
            "C02679": 1.0
        }
    }
    orient =
             1. first: write the first key
             2. second: write the second key first
    """
    f = open(outputfilename, 'w')

    if orient == 'first':
        for k1, g in data.items():
            for k2, v in g.items():
                f.write("'%s'.'%s'  %f\n" % (k1, k2, v))

    elif orient == 'second':
        for k1, g in data.items():
            for k2, v in g.items():
                f.write("'%s'.'%s'  %f\n" % (k2, k1, v))
    else:
        raise ValueError("orient must be either 'first' or 'second'!")

    f.close()

    return 1


if __name__ == '__main__':

    data_filepath = "data/"

    # rewrite reaction file

    all_rxn = convert_set_to_list(
        os.path.join(
            data_filepath,
            'reactions_modified.txt'))
    #convert_set_to_list(all_rxn, 'reactions_modified_noquotes.txt')
