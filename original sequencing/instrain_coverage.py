#!/home/ederrick/virtual_envs/inStrain/bin/python

import pandas as pd
import pickle
import inStrain
import inStrain.SNVprofile
import inStrain.profile.profile_utilities

profile_list = ["subsamp_I4_pulse1_instrain", "subsamp_I8_pulse1_instrain", "subsamp_K1_pulse1_instrain", "subsamp_L2_pulse1_instrain", "subsamp_L3_pulse1_instrain",
                "subsamp_L4_pulse1_instrain", "subsamp_L6_pulse1_instrain", "subsamp_L7_pulse1_instrain", "subsamp_L8_pulse1_instrain", "subsamp_K1_pulse0_instrain",
                "subsamp_L3_pulse0_instrain", "subsamp_L4_pulse0_instrain", "subsamp_L2_pulse0_instrain"]

for profile in profile_list:

        IS = inStrain.SNVprofile.SNVprofile(profile)
        covT = IS.get('covT')
        all_scaf = {}
        scaffolds = open("all_contigs.txt", "r")
        for scaffold in scaffolds:
                scaf = scaffold.strip()
                try:
                    	cov_scaf = inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(covT[scaf])
                        all_scaf[scaf] = [cov_scaf]
                except KeyError:
                        pass
        output_file = profile + '_depth.pickle'
        with open(output_file, 'wb') as handle:
                pickle.dump(all_scaf, handle, protocol = pickle.HIGHEST_PROTOCOL)
