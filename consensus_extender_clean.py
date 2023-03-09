import os
import re
import io
import pandas as pd
import numpy as np
import sys
import getopt
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import DBSCAN
from Bio.Phylo.TreeConstruction import DistanceCalculator
import subprocess
import math
from itertools import chain
from Bio.Align import substitution_matrices
from Bio.SeqUtils import molecular_weight
from math import log, sqrt
from Bio.Align import substitution_matrices as sm

# Get path of the script, to use as base path also for the other in bash.
script_name = sys.argv[0]
script_basename = os.path.dirname(script_name)

spec = [['input', 'i', 1, "character"],
    ['output', 'o', 1, "character"],
    ['genome', 'g', 1, "character"],
    ['max_sequences', 'x', 1, "integer"],
    ['help', 'h', 0, "logical"],
    ['min_plurality', 'p', 1, "integer"],
    ['min_saturation', 's', 1, "integer"],
    ['cluster_factor', 'c', 1, "integer"],
    ['min_cluster', 'm', 1, "integer"],
    ['group_outliers', 'l', 0, "logical"],
    ['interactive', 'n', 0, "logical"],
    ['top_rounds', 'r', 1, "integer"],
    ['min_seqs_per_cluster', 'q', 1, "integer"],
    ['extend', 'e', 1, "integer"],
    ['identity', 'd', 1, "double"]]

optlist, args = getopt.getopt(sys.argv[1:], 'i:o:g:x:hp:s:c:m:lnr:q:e:d:'), [s[0] for s in spec]

opt = {}
for option in optlist[0]:
    long_name_i = [x[1] for x in spec].index(option[0][1:])
    long_name = [x for x in spec[long_name_i]][0]
    opt[long_name] = option[1]

if 'help' in opt:
    print(getopt.getopt(spec, usage=True))
    sys.exit(1)

if 'input' not in opt:
    print("-input missing")
    sys.exit(1)

if 'genome' not in opt:
    print("-genome missing")
    sys.exit(1)

if 'output' not in opt:
    print("-output missing")
    sys.exit(1)

if 'max_sequences' not in opt:
    opt['max_sequences'] = 200
else:
    opt['max_sequences'] = int(opt['max_sequences'])

if 'min_plurality' not in opt:
    opt['min_plurality'] = 40
else:
    opt['min_plurality'] = int(opt['min_plurality'])

if 'min_saturation' not in opt:
    opt['min_saturation'] = 80
else:
    opt['min_saturation'] = int(opt['max_sequences'])

if 'cluster_factor' not in opt:
    opt['cluster_factor'] = 10
else:
    opt['cluster_factor'] = int(opt['cluster_factor'])

if 'min_cluster' not in opt:
    opt['min_cluster'] = 8
else:
    opt['min_cluster'] = int(opt['min_cluster'])

if 'end_threshold' not in opt:
    opt['end_threshold'] = 100
else:
    opt['end_threshold'] = int(opt['end_threshold'])

if 'group_outliers' not in opt:
    opt['group_outliers'] = True
else:
    opt['group_outliers'] = int(opt['group_outliers'])

if 'interactive' not in opt:
    opt['interactive'] = False
else:
    opt['interactive'] = bool(opt['interactive'])

if 'top_rounds' not in opt:
    opt['top_rounds'] = 20
else:
    opt['top_rounds'] = int(opt['top_rounds'])

if 'extend' not in opt:
    opt['extend'] = 500
else:
    opt['extend'] = int(opt['extend'])

if 'identity' not in opt:
    opt['identity'] = 0.9
else:
    opt['identity'] = float(opt['identity'])


def read_maf(maf_file, output_path):
    name = re.sub(".*/", "", maf_file)
    name = re.sub("\\..*$", "", name)
    print("[+++] Processing:", name)
    # Avoid processing a maf already processed.
    if os.path.getsize(maf_file) != 0 and not (os.path.exists(os.path.join(output_path, name + "_alt_1.con.fa")) or os.path.exists(os.path.join(output_path, name + "_alt_0.con.fa")) or os.path.exists(os.path.join(output_path, name + ".con.fa"))): 
         # This function takes a fasta file path as input and returns a dictionary
         # with the sequence id as the key and the sequence as the value.
        sequences = {}
        with open(maf_file) as f:
            sequence = ''
            sequence_id = ''
            for line in f:
                if line.startswith('>'):
                    if sequence_id:
                        sequences[sequence_id] = sequence
                        sequence = ''
                    sequence_id = line.strip().lstrip('>')
                else:
                    sequence += line.strip()
            sequences[sequence_id] = sequence
        # Read sequences from MAF file.
        data = []
        for sequence_id, sequence in sequences.items():
         row = list(sequence)
         row.insert(0, name)
         data.append(row)
        headers = ['id'] + list(range(len(data[0]) - 1))
        df = pd.DataFrame(data, columns=headers)
        # Return data table with the name of the maf as name of the first column to pass it between functions.
        return df
    else:
        return None


def K2Pdistance(maf):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    distance_matrix = np.zeros((maf.shape[0], maf.shape[0]), dtype=float)
    ts = ['AG', 'GA', 'CT', 'TC']  # transition pairs
    tv = ['AC', 'CA', 'AT', 'TA', 'CG', 'GC']  # transversion pairs
    for i in range(maf.shape[0]):
        seq1 = "".join(maf.iloc[i, 1:]).upper()
        for j in range(maf.shape[0]):
            seq2 = "".join(maf.iloc[j, 1:]).upper()
            pairs = []
            # collect ungapped pairs
            for x in zip(seq1, seq2):
                if '-' not in x: pairs.append(x)
            ts_count = 0
            tv_count = 0
            length = len(pairs)
            transitions = ["AG", "GA", "CT", "TC"]
            transversions = ["AC", "CA", "AT", "TA",
                             "GC", "CG", "GT", "TG"]
            for (x, y) in pairs:
                if x + y in transitions:
                    ts_count += 1
                elif x + y in transversions:
                    tv_count += 1
            p = float(ts_count) / length
            q = float(tv_count) / length
            try:
                d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
            except:
                d = np.nan
            distance_matrix[j][i] = d
    return distance_matrix

def seq_clus(maf, opt):
    if maf is None:
        return None
    else:
        name = maf.iloc[0, 0]
        """d = kimura_distance(maf)"""
        d = K2Pdistance(maf)
        kd = np.nanmean(d)
        print("[+++] Kimura Distance of full alignment:", kd)
        if maf.shape[0] > opt['min_cluster']:  # Don't try to cluster if less than min_cluster sequences
            print("[+++] Proceeding to cluster", name)
            if maf.shape[0] > opt['max_sequences']:  # If more than max_sequences get a sample of them
                maf = maf.sample(n=opt['max_sequences'], axis=0)
                d = K2Pdistance(maf)
                print("[+++] Sampling", opt['max_sequences'], "sequences from the alignment.")
            num_seqs = maf.shape[0]
            # Calculating mp value for DBSCAN based on the clustering factor and the minimum size of cluster
            mp = num_seqs // opt['cluster_factor']
            mp = max(mp, opt['min_cluster'] // 2)
            print("[+++] Mpoints =", mp)
            # Calculating the eps for DBSCAN as the one that produces more "stable" clusters.
            # We define "stable as clusters that remain for a distance of at least 0.04 (3 eps "steps")
            # The idea is to ignore clusters that are happen just at specific values as outliers.
            points = []
            # Convert NAs to maximum distance.
            d[d==-0.]=0
            # min eps granularity of 0.2. Could be made a parameter. 
            for eps in np.arange(0.02, 1.00, 0.02):
                candidate = len(np.unique(DBSCAN(eps=eps, min_samples=mp, metric='precomputed').fit_predict(d)))
                points.append(candidate)
            # 3 means the cluster must extend at least over 3 eps candidate steps ("stable cluster") 
            # Could be made a parameter.
            best_points = []
            for point in set(points):
                if points.count(point) > 3:
                    best_points.append(point)
            nclusters = max(best_points) 
            if not opt['group_outliers']:
                print("NOT GROUP")
                eps = (points.index(nclusters)) * 0.02 - 0.02
            else: # If this version is used the 0 clusters are probably just junk and should be
                eps = (50 - (points[::-1].index(nclusters))) * 0.02 -0.02
            dbscan = DBSCAN(eps=eps, min_samples=mp, metric='precomputed')
            labels = dbscan.fit_predict(d)
            clusters = labels + 1 # add 1 to start clusters at 1 instead of -1 for noise points
            print("[+++] Eps:", eps)
            if len(set(points)) > 0:
                print("[+++] Detected clusters:")
                print(set(clusters))
            else:
                print("[+++] No clusters detected.")
            mafs = []
            for i in set(clusters):
                if i != 0:
                    clus = maf[clusters == i].reset_index(drop=True)
                    clus.iloc[:, 0] = name + '_alt_' + str(i)
                    mafs.append(clus)
            if (0 in set(clusters)) and (sum(clusters == 0) > mp-1) and (max(clusters) > 1) or (eps == 0.02):
                clus = maf[clusters == 0].reset_index(drop=True)
                clus.iloc[:, 0] = name + '_alt_0'
                mafs.append(clus)
        else:
            print("[+++] No need to cluster", name)
            maf.iloc[:, 0] = name + '_alt_1'
            mafs = [maf]
        return mafs, np.round(kd, 3)


def process_maf(list_r, opt):
    if list_r is None:
        return None
    else:
        stats = []
        list_r, kd = list_r
        for i in range(len(list_r)):
            r = list_r[i].T
            seqs = r.shape[1]
            cols = r.shape[0]
            name = r.iloc[0, 0]
            # Calculate number of empty places and create an index of rows (bases)
            r["res"] = r.isin(["-"]).sum(axis=1)
            r["ID"] = np.arange(len(r))
            # remove fully empty positions (can happen after clustering)
            r = r[r["res"] != seqs]
            l = r.shape[0]
            r["resA"] = (r == "a").sum(axis=1)
            r["resC"] = (r == "c").sum(axis=1)
            r["resG"] = (r == "g").sum(axis=1)
            r["resT"] = (r == "t").sum(axis=1)
            r["mb"] = r[["resA", "resC", "resG", "resT"]].max(axis=1)  # Majority base
            # "Saturation" is the ratio of bases forming the consensus compared to the total
            # number of sequences.
            r["satur"] = r["mb"] / seqs
            # Remove empty positions and positions with just one base
            r = r[r["res"] + 1 < seqs]
            # Get majority base.
            r["base"] = ""
            r.loc[r["resA"] == r["mb"], "base"] = "a"
            r.loc[r["resC"] == r["mb"], "base"] = "c"
            r.loc[r["resG"] == r["mb"], "base"] = "g"
            r.loc[r["resT"] == r["mb"], "base"] = "t"
            r = r.reset_index(drop=True)
            # Define edges
            # lm and rm are equivalent to saturation, but requiring that the previous (lm) or
            # next (rm) position forms the consensus too. 
            r['lm'] = 0
            for i in range(1, r.shape[0]):
                r.at[i, 'lm'] = np.nansum((r.iloc[i, :seqs] == r.at[i, 'base']) & (r.iloc[i - 1, :seqs] == r.at[i - 1, 'base'])) / seqs
            r['rm'] = r['lm'].shift(-1)
            m1 = r.loc[r['lm'] > 0.5]['lm'].mean()
            sd1 = r.loc[r['lm'] > 0.5]['lm'].std()
            cutpoint1 = m1 - sd1 * 2
            # Select edges (lm and rm could be unified as the max of both)
            min_r = r[(r['lm'] >= cutpoint1) | (r['rm'] >= cutpoint1)]
            # Correct edges for low number of seqs (need more confidence -> more consecutive bases)
            # Take the edges that match too this rule
            confidence_needed = 20//seqs + 1
            v = min_r.iloc[0:confidence_needed]['rm'] >= cutpoint1
            v = v.tolist()
            while sum(v) < confidence_needed:
                min_r = min_r[(confidence_needed - (len(v) - v[::-1].index(False)) + 1):]
                v = min_r[0:confidence_needed]['rm'] >= cutpoint1
                v = v.tolist()
            v = min_r[(len(min_r) - confidence_needed):(len(min_r))]['lm'] >= cutpoint1
            v = v.tolist()
            while sum(v) < confidence_needed:
                min_r = min_r[0:(len(min_r) - confidence_needed + v.index(False))]
                v = min_r[(len(min_r) - confidence_needed):(len(min_r))]['lm'] >= cutpoint1
                v = v.tolist()
            # Get the final set of bases that we will use to build the consensus.
            minmin_r = r[(r['ID'] >= min(min_r['ID'])) & (r['ID'] <= max(min_r['ID']))]
            # kimura distance of the final alignment BUT STILL WITH INTERNAL LOW COUNT BASES THAT WILL NOT BE PART OF THE CONSENSUS.
            #d = kimura_distance(minmin_r.iloc[:, 0:seqs].T)
            d = K2Pdistance(minmin_r.iloc[:, 0:seqs].T)
            kdc = np.nanmean(d)
            print("[+++] Kimura distance for {}: {}".format(name, kdc))
            # Build the consensus with a simple plurality model. It can be made as sophisticated as needed.
            minmin_r = minmin_r[minmin_r['satur'] > (opt['min_plurality'] / 100)].reset_index(drop=True)
            s = ''.join(minmin_r['base'])
            # Draw plots and open cluster alignments in AliView when in interactive mode.
            if opt['interactive']:
                plt.plot(r['ID'], r['lm'])
                plt.title(name)
                plt.xlabel("left = {} - right = {}".format(min(minmin_r['ID']), max(minmin_r['ID'])))
                plt.axhline(y=cutpoint1, color='blue')
                plt.axvline(x=min(minmin_r['ID']), color='red')
                plt.axvline(x=max(minmin_r['ID']), color='green')
                plt.axhline(y=(opt['min_plurality'] / 100), color='black')
                plt.savefig("plot_{}.png".format(name))
                plt.show()
                #aliview(as.DNAbin(t(minmin_r.iloc[:, 0:(ncol(minmin_r) - 9)])))
                # aliview(as.DNAbin(t(r.iloc[:, 0:(ncol(minmin_r) - 9)])))
                input("Press a key to continue")
            min_l = minmin_r.shape[0]
            if kd is None:
                kd = float('inf')
            stats.append(pd.DataFrame({
                'seq': [name],
                'seqs_for_consensus': [seqs],
                'alig_size': [l],
                'cons_size': [len(s)],
                'start_con': [minmin_r.loc[0, 'ID']],
                'end_con': [minmin_r.loc[min_l-1, 'ID']],
                'surplus_left': [r[r['ID'] < minmin_r.loc[0, 'ID']].shape[0]],
                'surplus_right': [r[r['ID'] > minmin_r.loc[min_l-1, 'ID']].shape[0]],
                'dust_left': [1-(np.sum(r[r['ID'] < minmin_r.loc[0, 'ID']]['res'])/(seqs*r[r['ID'] < minmin_r.loc[0, 'ID']].shape[0]))],
                'dust_right': [1-(np.sum(r[r['ID'] > minmin_r.loc[min_l-1, 'ID']]['res'])/(seqs*r[r['ID'] > minmin_r.loc[min_l-1, 'ID']].shape[0]))],
                'maf_kd': [kd],
                'cluster_kd': [kdc],
                'end_l': [r[r['ID'] < minmin_r.loc[0, 'ID']].shape[0] - np.sum((r[r['ID'] < minmin_r.loc[0, 'ID']]['res']+r[r['ID'] < minmin_r.loc[0, 'ID']]['mb'])==seqs) - (np.sum(r[r['ID'] < minmin_r.loc[0, 'ID']]['res'])/seqs) > opt['end_threshold']],
                'end_r': [r[r['ID'] > minmin_r.loc[min_l-1, 'ID']].shape[0] - np.sum((r[r['ID'] > minmin_r.loc[min_l-1, 'ID']]['res']+r[r['ID'] > minmin_r.loc[min_l-1, 'ID']]['mb'])==seqs) - (np.sum(r[r['ID'] > minmin_r.loc[min_l-1, 'ID']]['res'])/seqs) > opt['end_threshold']]
            }))
            # Write fasta for final result
            if stats[-1].loc[0, 'end_l'] and stats[-1].loc[0, 'end_r']:
                with open(opt['output'] + "/final/" + name + ".con.fa", 'w') as f:
                    f.write(">" + name + "\n" + s + "\n")
            else:
                end = False
                with open(opt['output'] + "/to_extend_" + str(round) + "/" + name + ".con.fa", 'w') as f:
                    f.write(">" + name + "\n" + s + "\n")
    return stats

# Main body

if __name__ == "__main__":
    round = 0 # Counter of extension iterations
    end = False # TRUE if no more sequences to extend.
    os.makedirs(os.path.join(".", opt["output"]), exist_ok=True)
    os.makedirs(os.path.join(".", opt["output"], f"final"), exist_ok=True)

    while (round < opt["top_rounds"]) and not end:
        round += 1
        os.makedirs(os.path.join(".", opt["output"], f"mafs_{round}"), exist_ok=True)
        os.makedirs(os.path.join(".", opt["output"], f"to_extend_{round}"), exist_ok=True)

        print(f"[+++] STARTING ROUND {round}")
        print(f"[+++] Sequences folder: {opt['input']}")
        print(f"[+++] Mafs folder: {os.path.join(opt['output'], f'mafs_{round}')}")

        # Fasta table for the special first case (all should try to be extended)
        if round == 1:
            cmd = f"cat {opt['input']}/*.fa | awk '$0 ~ \">\" {{if (NR > 1) {{print c \"\\tFALSE\\tFALSE\";}} c=0;printf substr($0,2,100) \"\\t\"; }} $0 !~ \">\" {{c+=length($0);}} END {{ printf c \"\\tFALSE\\tFALSE\";}}'"
            buffer = io.StringIO(subprocess.check_output(cmd, shell=True, text=True))
            fasta_table = pd.read_csv(buffer, sep="\t", names=["seq", "cons_size", "end_l", "end_r"])
        # Read fa files and generate their maf. Uses a modified version of make_align_from_blast from TE_AID
        filenames = [f for f in os.listdir(opt['input']) if re.match(r'.*\.fa', f)]
        for i in range(len(filenames)):
            cmd = f"../../{script_basename}/make_align_from_blast_alt.sh {opt['genome']} ../../{opt['input']}/{filenames[i]} {fasta_table.iloc[[i]]['cons_size'].item()*opt['identity']} {opt['extend']} {fasta_table.iloc[[i]]['end_l'].item()} {fasta_table.iloc[[i]]['end_r'].item()}"
            answer = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=f"{opt['output']}/mafs_{round}", check=True)
            print(answer.stdout)
        # Delete potential empty files on the maf folder.
        subprocess.run("find . -size 0 -print -delete", shell=True, cwd=f"{opt['output']}/mafs_{round}", check=True)

        filenames = os.listdir(opt["output"]+"/mafs_"+str(round))
        filenames = [opt["output"]+"/mafs_"+str(round)+"/"+x for x in filenames if x.endswith(".maf.fa")]

        # Divide in sets of alignments to avoid overflow processing them. Should be a variable
        lf = np.array_split(filenames, math.ceil(len(filenames)/20))
        stats_all_list = []

        # Main loop for round: reads maf, generates a list of sub mafs, and finally
        # writes the consensus for each one and returns a list of stats about them.
        for i in lf:
            stats = [process_maf(seq_clus(read_maf(x, opt["output"]), opt), opt) for x in i]
            for k in range(len(stats)):
                for l in range(len(stats[k])):
                    stats_all_list.append(stats[k][l])
        # Write the final stats
        stats_all = stats_all_list[0].iloc[:, :]
        for i in range(1, len(stats_all_list)):
            stats_all = pd.concat([stats_all, stats_all_list[i].iloc[:, :]])

        with open(os.path.join(opt['output'], f'stats_round_{round}.tsv'), 'a') as f:
            stats_all.to_csv(f, sep='\t', index=False, header=True)
        print(f"[+++] Round {round} completed. {stats_all.shape[0]} consensus sequences generated.")
        fasta_table = stats_all.loc[:, ["seq", "cons_size", "end_l", "end_r"]]
        fasta_table = fasta_table[~((stats_all['end_l'] == True) & (stats_all['end_r'] == True))].reset_index(drop=True)
        opt['input'] = f"{opt['output']}/to_extend_{round}"
        if fasta_table.shape[0] == 0:
            end = True