import pandas as pd
import re
from config import *


def get_proper_samples():
    """
    :return:  dictionary of {(trail,patient): [(TimepointX, TimepointY),...], ...}
    """
    comparisons = {}
    trials = {}
    result_proper = open(RESULT_PROPER_FILE, "r")
    for row in result_proper.readlines():
        row = row.rstrip(os.linesep)
        row = row.replace('//', '/')
        if PSEUDOMONAS in row:
            if "patient" in row:
                file_path = os.path.join(row, 'trees', 'formatted_comparisons_converted.txt')
                parse_formatted_comparisons_converted(file_path, comparisons, trials)
            else:
                trial_num = row.split('_')[-1]
                for directory in sorted(os.listdir(HIGH_LOST_PATH_P)):
                    if "trial_" + trial_num + "_" in directory:
                        file_path = os.path.join(HIGH_LOST_PATH_P, directory, 'trees',
                                                 'formatted_comparisons_converted.txt')
                        parse_formatted_comparisons_converted(file_path, comparisons, trials)
    return comparisons, trials


def parse_formatted_comparisons_converted(file_path, comparisons, trials):
    trial = int(file_path[:-42].split('_')[-3])
    patient = int(file_path[:-42].split('_')[-1])
    # get the trial numbers and patients
    if trial not in trials.keys():
        trials[trial] = [patient]
    else:
        trials[trial].append(patient)
    # get the proper timepoints
    comps = []
    try:
        with open(file_path, 'r') as f:
            for row in f.readlines():
                comps.append((row.split()[0], row.split()[1]))
        comparisons[(trial, patient)] = comps
    except IOError:
        pass


def clean_s1():
    # cleaning the relevant patients and trials (from pseudomonas only):
    _, trials = get_proper_samples()
    df = pd.read_csv(TABLE_S1, sep='\t')
    df = df[df.Bacteria == "Pseudomonas aeruginosa"].reset_index()
    cleaned_df = df
    patients = []
    for idx_row, row in cleaned_df.iterrows():
        experiment = row["Experiment"]
        if experiment not in trials.keys():
            cleaned_df = cleaned_df.drop(idx_row, axis=0)
        else:
            original_exp_index = df[df.Experiment == experiment].index[0]
            patient_index = idx_row - original_exp_index + 1
            if (patient_index not in trials[experiment]):
                cleaned_df = cleaned_df.drop(idx_row, axis=0)
            else:
                patients.append(patient_index)

    cleaned_df = cleaned_df.reset_index()
    cleaned_df.insert(loc=3, column="patient", value=patients)
    cleaned_df = cleaned_df.drop(["level_0", "index"], axis=1)
    columns = list(cleaned_df.columns)

    # changing the timepoints to be uniformed
    for idx_row, row in cleaned_df.iterrows():
        for idx_col, col in enumerate(row):
            if type(col) is str and '/' in col:
                if re.findall(re.compile(r"\-(SRR.*)\/"), col):
                    cleaned_df.at[idx_row, columns[idx_col]] = re.findall(re.compile(r"\-(SRR.*)\/"), col)[0]
                elif '/SRR' in col or '/ERR' in col:
                    cleaned_df.at[idx_row, columns[idx_col]] = col.split("/")[-1]
            elif re.findall(re.compile(r"(Sommer_setup).*"), str(col)):
                new_val = re.findall(re.compile(r"(Sommer_setup).(.*)"), str(col))[0]
                cleaned_df.at[idx_row, columns[idx_col]] = new_val[0] + "_" + new_val[1]

    return cleaned_df


def calc_time_intervals(df, comparisons):
    """
    :param df: the cleaned data frame
    :param comparisons: dictionary of {(trail,patient): [(TimepointX, TimepointY),...], ...}
    :return: updated comparisons such that each tuple pair became a triple containing the timeinterval in days between 
       timepoint X to timepoint Y: {(trail,patient): [(TimepointX, TimepointY, intervalInDays),...], ...}
    """
    time_df = df.filter(regex=("Timepoint.*"))
    for key, pairs_list in comparisons.iteritems():
        for idx, pair in enumerate(pairs_list):
            prev_row = time_df.loc[time_df.where(time_df == pair[0]).dropna(how='all').dropna(axis=1).index]
            prev_col = time_df.where(time_df == pair[0]).dropna(how='all').dropna(axis=1).columns
            prev_time = prev_row[prev_col[0][:-9] + "time"].iloc[0]
            curr_col = time_df.where(time_df == pair[1]).dropna(how='all').dropna(axis=1).columns
            curr_time = prev_row[curr_col[0][:-9] + "time"].iloc[0]
            # print pair[0], prev_time, pair[1], curr_time
            pairs_list[idx] += (time_converter(curr_time) - time_converter(prev_time),)

    return comparisons


def time_converter(time_string):
    """
    convert time (string) to days (int)
    """
    if time_string.endswith('y'):
        return int(time_string[:-1]) * 365
    elif time_string.endswith('m'):
        return int(time_string[:-1]) * 30
    elif time_string.endswith('d'):
        return int(time_string[:-1])
    else:
        return -1


def get_network_df(comparisons_time_intervals):
    """
    :return: pandas dataframe where each row is a comparison between two timepoints 
    """
    df = pd.DataFrame(
        columns=["trial", "patient", "t1", "t2", "time_interval", "t1_exist_genes", "lost_genes_t1_to_t2"])
    for key, pairs_list in comparisons_time_intervals.iteritems():
        trial = key[0]
        patient = key[1]
        for triple in pairs_list:
            t1, t2, time_interval = triple
            patient_path = "Data/Pseudomonas_aeruginosa/trial_{0}_patient_{1}/breseq".format(
                trial, patient)
            t1_exist_genes_path = os.path.join(patient_path, t1, "output", "{0}_high_non_lost_genes.txt".format(t1))
            with open(t1_exist_genes_path, "r") as f:
                t1_exist_genes = [line[:-1] for line in f.readlines()]

            t2_exist_genes_path = os.path.join(patient_path, t2, "output", "{0}_high_non_lost_genes.txt".format(t2))
            with open(t2_exist_genes_path, "r") as f:
                t2_exist_genes = [line[:-1] for line in f.readlines()]

            lost_genes = set(t1_exist_genes) - set(t2_exist_genes)
            df = df.append({"trial": key[0], "patient": key[1], "t1": t1, "t2": t2, "time_interval": time_interval,
                            "t1_exist_genes": t1_exist_genes, "lost_genes_t1_to_t2": lost_genes}, ignore_index=True)

    df = df.sort_values(["trial", "patient"])
    df.to_csv(RAW_DATA_P, index=False)
    return df
