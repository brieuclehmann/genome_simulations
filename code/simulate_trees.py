import os
import argparse
import sys
import pandas as pd
import numpy as np
import msprime
import demes
import tskit
from scipy.linalg import eigh

sys.path.insert(0, "code")
import pedigree_tools
import grm_tools

def main(args):
    ###################
    #### Load data ####
    ###################

    # input pedigree file name
    pedigree_file_name = args.pedigree_name
    # load input pedigree
    txt_ped = pedigree_tools.load_and_verify_pedigree(pedigree_file_name)

    # Two population out-of-Africa
    # demography model from Tennessen et al., 2012
    yaml_file = '{}/code/Tennessen_ooa_2T12_rescale.yaml'.format(args.dir)
    graph = demes.load(yaml_file)
    ooa_2T12 = msprime.Demography.from_demes(graph)
    # define recombination map file name from stdpopsim
    map_file_name = '{}/code/maps/genetic_map_GRCh37_chr{}.txt'.format(args.dir, args.chromosome)
    map = msprime.intervals.RateMap.read_hapmap(map_file_name)
    # identify centromere
    centromeres=pd.read_csv('{}/code/maps/centromeres.csv'.format(args.dir))
    centromere_intervals=centromeres.loc[centromeres['chrom'] == str(args.chromosome),["start","end"]].values
    # get chromosome length
    # from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
    chromosome_length=pd.read_csv('{}/code/maps/GRCh37_chromosome_length.csv'.format(args.dir))
    assembly_len=chromosome_length.loc[chromosome_length['chrom'] == str(args.chromosome),["length_bp"]].values
    # sequence_length from map file
    len =  map.right.max()

    print('pedigree: {}, chromosome: {}, length: {}, mutation_rate: {}, censor: {}'.format(
          args.pedigree_name, args.chromosome, len, args.mut_rate, args.censor))

    ################################
    #### Run genome simulations ####
    ################################
    ts, ts_recap = pedigree_tools.simulate_genomes_with_known_pedigree(
         text_pedigree = txt_ped,
         demography = ooa_2T12,
         model = 'hudson',
         mutation_rate = args.mut_rate,
         rate_map = map,
         sequence_length = len,
         sequence_length_from_assembly = assembly_len,
         centromere_intervals = centromere_intervals,
         censor = args.censor,
         seed = args.repetition + args.seed_offset
         )

    # run some basic sanity checks
    pedigree_tools.simulation_sanity_checks(ts, txt_ped)
    pedigree_tools.simulation_sanity_checks(ts_recap, txt_ped)

    # Save output tree sequence
    ts_dir_out = 'tree_sequences/balsac/chr{}/sim{}'.format(
        args.chromosome, args.repetition)
    if not os.path.exists(ts_dir_out):
        os.makedirs(ts_dir_out)

    ts.dump(ts_dir_out + "/trees.ts")
    ts_recap.dump(ts_dir_out + "/trees_recap.ts")

    print("Finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument("-d", "--dir",
        help="directory of input text pedigree"
        )
    parser.add_argument("-o", "--out",
        help="directory of output files"
        )
    parser.add_argument("-ped", "--pedigree_name",
        help="name of input text pedigree (no file extention)"
        )
    # optional arguments
    parser.add_argument("-sfx", "--suffix",
        default="sim",
        help="output file suffix"
        )
    parser.add_argument("-censor", "--censor",
        action="store_true",
        help="removes pedigree information from tree sequence"
        )
    parser.add_argument("-chr", "--chromosome",
        type=int,
        help="specify chromosome number to be simulated"
        )
    parser.add_argument("-m", "--mut_rate",
        default=1.66e-8,
        type=float,
        help="specify mutation rate"
        )
    parser.add_argument("-seed", "--seed_offset",
        default=0,
        type=int,
        help="specify random seed"
        )
    parser.add_argument("-rep", "--repetition",
        default=1,
        type=int,
        help="specify repetition"
        )
    args = parser.parse_args()

    main(args)
