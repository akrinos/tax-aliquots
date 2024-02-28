# Python script to read in cluster data from a certain DeepClust cluster,
# read in sequences, and compute intersections between all given sequences.
# 
#
#from ...dependencies.faSomeRecords.faSomeRecords import main as fasomerecords
import faSomeRecords as fasomerecords
import os
import pandas as pd
import numpy as np
import json
import sys
import time
import subprocess
from Bio import SeqIO
import itertools
import dask
from dask import delayed
from dask.distributed import Client
import difflib
import scipy
import random
import pysais
from difflib import SequenceMatcher
from sklearn.cluster import AgglomerativeClustering
from sklearn.datasets import make_blobs
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import ahocorasick
from ahocorasick import Automaton
import argparse, textwrap
import parsl
from parsl import python_app, bash_app

@python_app
def read_in_sequences(behaved_fasta,list_of_seqs):
    sequence_list = fasomerecords.main(behaved_fasta,list_of_seqs,True)
    recs = SeqIO.parse(sequence_list,"fasta")
    return recs

@python_app
def compute_intersections(strings,kmer_len):
    interecting_strings=dict()
    len_intersec_strings=[]
    n=int(kmer_len)
    keywords = [string[i:i+n] for string in strings for i in range(len(string) - n + 1)]

    automaton = Automaton()
    for i, keyword in enumerate(keywords):
        automaton.add_word(keyword, (i, keyword))
    automaton.make_automaton()

    for i, string in enumerate(strings):
        interecting_strings[string] = set([curr[1][1] for curr in list(automaton.iter(string))])
        len_intersec_strings.append(len(interecting_strings[string]))
    intersection_info_curr = np.zeros((len(strings),len(strings)))
    
    
    for i,string1 in enumerate(strings):
        for j,string2 in enumerate(strings):
            min_length=min(len(interecting_strings[string1]),
                           len(interecting_strings[string2]))
            if min_length==0:
                min_length=1
            intersection_info_curr[i,j] = (min_length-len(interecting_strings[string1].\
                                          intersection(interecting_strings[string2])))/min_length


    intersection_info_curr  = np.triu(intersection_info_curr) + np.triu(intersection_info_curr,1).T
    intersection_info_curr = pd.DataFrame(intersection_info_curr)

    intersection_info_curr.index=cluster_seq_ids
    intersection_info_curr.columns=cluster_seq_ids
    return intersection_info_curr
    
@python_app
def get_list_of_seqs(behaved_fasta,list_of_seqs):
    return [curr.seq for curr in read_in_sequences(behaved_fasta,list_of_seqs)]
    
@python_app
def hier_clust(clust_rep,subfile_clust,kmer_len,
               intermed_file_dir,hom_clust_id,behaved_fasta,
               tax_table):
    strings=get_list_of_seqs(behaved_fasta,subfile_clust.loc[subfile_clust.Rep==clust_rep,"Mem"])
    hier_clust_input=compute_intersections(strings,kmer_len)
    hier_clust_input.index=[curr.replace(".","N") for curr in hier_clust_input.index]
    sequence_names = hier_clust_input.index
    hier_clust_input.index.rename('sequence', inplace=True)

    dataframe_intersections_unmelt= hier_clust_input.merge(subfile_clust,
                                                           left_on="sequence",
                                                           right_on="Mem")
    
    dataframe_intersections_unmelt = dataframe_intersections_unmelt.merge(tax_table.drop_duplicates(),
                                                                          left_on="Source_ID",
                                        right_on="Source_ID",how="left")

    if len(dataframe_intersections_unmelt.index) > 1:
        to_drop = list(tax_table.columns)
        a=np.array(dataframe_intersections_unmelt.\
                   drop(to_drop,axis="columns"))
        np.fill_diagonal(a,0)
        Z = linkage(scipy.spatial.distance.squareform(a),
                method='weighted')
        max_d=float(distance_threshold)
        clusters = fcluster(Z, max_d, criterion='distance')
        dataframe_intersections_unmelt["cluster_blob"] = clusters
    else:
        dataframe_intersections_unmelt["cluster_blob"] = 1

    dataframe_intersections_unmelt["HomClust"]=hom_clust_id.split("/")[-1].split(".")[0]
    dataframe_intersections_unmelt["SeqNames"] = sequence_names
    dataframe_intersections_unmelt.loc[:,["Source_ID","HomClust","cluster_blob","Division","Domain"]].\
        to_csv(os.path.join(intermed_file_dir,str(hom_clust_id)+"_k_"+str(kmer_len)+".intermed"))
    return hom_clust_id


# main function
def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="process_clusters.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="Runs the hierarchical clustering steps of tax-aliquots following DeepClust clustering.",
        epilog=textwrap.dedent('''\
    Check out the documentation for each of the flags.
    '''))
    parser.add_argument(
    '--combined_fasta', '-c', metavar='combined_fasta', type=str, required=True,
    help='FASTA file where the combination of all the sequences is.')
    parser.add_argument(
    '--output_dir', metavar='output_dir', type=str,
    help='Directory in which to write output.')
    parser.add_argument(
    '--output_file', metavar='output_file', type=str,
    help='Name of the output files.')
    parser.add_argument(
    '--cluster_file', '-i', metavar='cluster_file', type=str,
    help='Name of the file where DeepClust clusters are stored.')
    parser.add_argument(
    '--tax_table', '-t', metavar='tax_table', type=str,
    help='Location of the database taxonomy table.')
    parser.add_argument(
    '--prot_map', '-p', metavar='prot_map', type=str,
    help='Location of the json protein map for the database.')
    parser.add_argument(
    '--level_of_interest', '-l', metavar='level_of_interest', default="domain", type=str,
    help='Level of taxonomy to operate on.')
    parser.add_argument(
    '--list_of_interest', metavar='list_of_interest', default="domain", type=str,
    help='Underscore-separated list of taxa of interest (e.g. Skeletonematacae_Rhizosoleniaceae).')
    parser.add_argument(
    '--sample_id', '-s', metavar='sample_id', default="Samp", type=str,
    help='Sample delimiter in sequence file (e.g. NarBay).')
    parser.add_argument(
    '--kmer_len', '-k', metavar='kmer_len', default=3, type=int,
    help='kAAmer length to use for algorithm.')
    args = parser.parse_args()

    # read in taxonomy table 
    tax_table = pd.read_csv(args.tax_table,sep="\t")
    tax_table=tax_table.drop([curr for curr in tax_table.columns if "Unnamed" in curr],
                             axis="columns")
    tax_table.columns = [curr.lower() for curr in tax_table.columns]
    tax_table.source_id = tax_table.source_id.astype(str)
    
    # read in cluster file
    clusters_read=pd.read_csv(args.cluster_file,sep="\t",header=None,names=["Rep","Mem"])
    
    # reformat FASTA file
    os.system(" ".join(["mkdir -p",args.output_dir]))
    behaved_fasta = os.path.join(args.output_dir,
                                 "behaved_fasta_"+str(args.sample_id)+".fasta")
    os.system(" ".join(["bash","remove_newlines.sh",args.combined_fasta,">",behaved_fasta]))
    
    # make protein file linebreak-separated so we can just rg from it
    prot_map_lb = os.path.join(args.output_dir,
                               "linebreak_"+str(args.sample_id)+".fasta")
    os.system(" ".join(["cat",args.prot_map,"|","tr",",",repr(os.linesep),"|",
                        "sed","'s/\"//g'","|",
                        "sed","'s/{//g'","|","sed","'s/}//g'","|","sed","'s/ //g'","|",
                        "tr",":",",",">",
                        prot_map_lb]))
    
    print(" ".join(["cat","args.protmap","|","tr",",",repr("\n"),"|",
                        "sed","'s/\"//g'","|",
                        "sed","'s/{//g'","|","sed","'s/}//g'","|","sed","'s/ //g'","|",
                        "tr",":",",",">",
                        "protmaplb"]),flush=True)
   
    ## filter out organism of interest
    source_ids=list(set(tax_table.loc[tax_table[args.level_of_interest].\
                                      isin(str(args.list_of_interest).split("_")),"source_id"]))
    
    ## idiosyncracies of the problems I have with my formatted labels
    if args.sample_id == "CAMPEP":
        clusters_read["Stem"] = clusters_read["Mem"]
    elif (args.sample_id == "P_")|(args.sample_id=="k"):
        clusters_read[["Stem","Leaf"]] = clusters_read["Mem"].str.split(pat="\.p",expand=True,n=2)
        clusters_read["Stem"] = clusters_read.Stem.str.replace(".","N") 
    else:
        clusters_read[["Stem","Leaf"]] = clusters_read.Mem.str.split(pat="_CONTIGID",expand=True,n=2)
        
    ## let's extract our situation of interest
    temp_file=os.path.join(args.output_dir,"temp_file.txt")
    source_id_fr=pd.read_csv(prot_map_lb,sep=",", header=None, names=["Stem","Source_ID"])
    source_id_fr=source_id_fr.merge(clusters_read,how="inner",left_on="Stem",
                                    right_on="Stem")
    
    clusters_read = clusters_read.merge(source_id_fr,how="left",
                                        left_on="Stem",right_on="Stem")
    os.system("rm "+temp_file)
    
    ## get just the clusters that align with our interest
    interest_mem=clusters_read.loc[clusters_read.Source_ID.isin(source_ids)]
    interest_clust=clusters_read.loc[clusters_read.Rep.isin(list(interest_mem.Rep))]
    reps_with_samp = interest_clust.loc[[sample_id in curr for curr in interest_clust.Mem]]

    clust_rep_set=list(set(reps_with_samp.Rep))
    cluster_ids = ["HomClust_"+str(curr) for curr in range(len(clust_rep_set))]
    cluster_id_dict = dict(zip(clust_rep_set,cluster_ids))
    
    intermed_file_dir=os.path.join(args.output_dir,
                                   args.output_file+"_intermed")
    os.system("mkdir -p "+intermed_file_dir)
    completion_nums = []
    for clust_rep in range(clust_rep_set):
        subfile_clust=clusters_read.loc[clusters_read.Rep==clust_rep,:]
        hom_clust_id=cluster_id_dict[clust_rep]
        completion_nums.append(hier_clust(clust_rep,subfile_clust,args.kmer_len,
                                          intermed_file_dir,hom_clust_id,
                                          behaved_fasta,tax_table))

    # Main function: wait for all apps to finish and collect the results
    outputs = [i.result() for i in completion_nums]

## EXAMPLE COMMAND:
## python process_clusters.py -i /vortexfs1/omics/alexander/akrinos/2022-euk-diversity/code/snakemake-workflows/diamond-deepclust/deepclust-workflow/deepclust_contigname.mad.50.out -c /vortexfs1/omics/alexander/akrinos/2022-euk-diversity/code/snakemake-workflows/hierarchical-clustering-samples/combined_seqs.fasta --output_file="test" --output_dir="test_dir" -t /vortexfs1/omics/alexander/data/databases/marmmetsp-5Dec2022/tax-table.txt --prot_map=/vortexfs1/omics/alexander/data/databases/marmmetsp-5Dec2022/prot-map.json --kmer_len=3 -s test --level_of_interest family --list_of_interest Hemiaulaceae_Rhizosoleniaceae_Thalassiosiraceae_Skeletonemataceae
if __name__ == "__main__":
    main()