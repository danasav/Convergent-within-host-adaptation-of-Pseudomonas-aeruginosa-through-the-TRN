"""
Convert dict of TF:targets from 4-letter annotations/locus annotations to RefSeq Protein annotations
"""

import gzip as gz
from Bio import SeqIO
import warnings


def parse_cds(cds_file, extra_file):
    """
    Parse CDS file to dict with gene_name:refseq and locus:refseq
    """
    # initialize
    cds_dict = {}
    # parse extra file
    with open(extra_file, 'rt') as fl:
        for i in fl.readlines():
            if len(i.split()) > 1: cds_dict[i.strip().split()[0].lower()] = i.strip().split()[1]
    # open file
    with gz.open(cds_file, 'rt') as fl:
        seqs = list(SeqIO.parse(fl, "fasta"))
    # iterate
    for seq in seqs:
        # add locus
        locus_tag = [i for i in seq.description.split() if "locus_tag" in i][0].split("=")[1].replace("]", "").lower()
        # add refseq name
        # print locus_tag
        protein_id = [i for i in seq.description.split() if "protein_id" in i][0].split("=")[1].replace("]", "")
        # add to dict
        cds_dict[locus_tag] = protein_id
        # add gene if available
        pregene = [i for i in seq.description.split() if "gene=" in i]
        if len(pregene) > 0:
            # gene name
            gene = pregene[0].split("=")[1].replace("]", "").lower()
            # add to dict
            cds_dict[gene] = protein_id
        # add if protein X] in description
        for spliter in ["protein", "transporter"]:
            preprot = seq.description.split(" %s " % spliter)
            if len(preprot) == 2:
                prot = preprot[1].split("]")[0].lower()
                if prot not in cds_dict: cds_dict[prot] = protein_id
    # return
    return cds_dict


def parse_locus_tag(cds_file, extra_file):
    """
    Parse CDS file to dict with gene_name:refseq and locus:refseq
    """
    # initialize
    cds_dict = {}
    # parse extra file
    with open(extra_file, 'rt') as fl:
        for i in fl.readlines():
            if len(i.split()) > 1: cds_dict[i.strip().split()[0].lower()] = i.strip().split()[1]
    # open file
    with gz.open(cds_file, 'rt') as fl:
        seqs = list(SeqIO.parse(fl, "fasta"))
    # iterate
    for seq in seqs:
        # add locus
        locus_tag = [i for i in seq.description.split() if "locus_tag" in i][0].split("=")[1].replace("]", "").lower()
        # add refseq name
        # print locus_tag
        protein_id = [i for i in seq.description.split() if "protein_id" in i][0].split("=")[1].replace("]", "")
        # add to dict
        cds_dict[locus_tag] = protein_id
    return cds_dict


def convert_to_refseq(tf_dict, cds_dict):
    """
    Convert all names in tf_dict to their refseq protein counterparts
    """
    # initialize
    converted_tf_dict = {}
    # iterate through keys
    for TF in tf_dict:
        # convert tf
        try:
            converted_tf = [cds_dict[TF.lower()]]
        except  KeyError:
            if TF == "ihf":
                converted_tf = [cds_dict[i] for i in ["hima", "himd"]]
            else:
                # print TF
                a = [i for i in cds_dict if TF.lower() in i]
                if len(a) > 1: print(a)
                continue
        if len(converted_tf) == 1 and converted_tf[0] in ["None", "ncRNA"]: continue
        # convert targets
        converted_targets = []
        for target in tf_dict[TF]:
            # convert target
            try:
                converted_target = [(cds_dict[target[0].lower()], target[1])]
            except KeyError:
                if target[0] == "cupB":
                    converted_target = [(cds_dict[i], target[1]) for i in
                                        ["cupb1", "cupb2", "cupb3", "cupb4", "cupb5", "cupb6"]]
                elif target[0] == "cupC":
                    converted_target = [(cds_dict[i], target[1]) for i in ["cupc1", "cupc2", "cupc3"]]
                elif target[0] == "cupA":
                    converted_target = [(cds_dict[i], target[1]) for i in ["cupa1", "cupa2", "cupa3", "cupa4", "cupa5"]]
                elif target[0] == "phzG":
                    converted_target = [(cds_dict[i], target[1]) for i in ["phzg1", "phzg2"]]
                elif target[0] == "fumC":
                    converted_target = [(cds_dict[i], target[1]) for i in ["fumc1", "fumc2"]]
                elif target[0] == "astA":
                    converted_target = [(cds_dict[i], target[1]) for i in ["aruf", "arug"]]
                else:
                    warnings.warn("%s not properly in cds dict" % target[0])
                    a = [i for i in cds_dict if target[0].lower() in i]
                    if len(a) > 1: warnings.warn(
                        "%s found as substring of multiple entries in the cds dict though: %s" % (target[0], a))
                    continue
            # add to targets
            if len(converted_target) == 1 and converted_target[0][0] in ["None", "ncRNA"]: continue
            converted_targets += converted_target
        # add to results dict
        for tf in converted_tf:
            converted_tf_dict[tf] = converted_targets
    # return
    return converted_tf_dict


def main(tf_dict, cds_file, extra_file):
    # create cds_dict
    cds_dict = parse_cds(cds_file, extra_file)
    # convert dict to refseq names
    converted_tf_dict = convert_to_refseq(tf_dict, cds_dict)
    # return
    return converted_tf_dict


if __name__ == "__main__":
    raise Exception("This script is not meant to be called from the command line.")
