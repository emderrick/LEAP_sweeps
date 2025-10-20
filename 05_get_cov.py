#!/mfs/ederrick/miniconda3/envs/instrain/bin/python

import pandas as pd
import pickle
import inStrain
import inStrain.SNVprofile
import inStrain.profile.profile_utilities

profile_list = ["LEAP_META_01_T1_refined_inStrain","LEAP_META_02_T1_refined_inStrain","LEAP_META_03_T1_refined_inStrain","LEAP_META_04_T1_refined_inStrain","LEAP_META_05_T1_refined_inStrain",
		"LEAP_META_06_T1_refined_inStrain","LEAP_META_07_T1_refined_inStrain","LEAP_META_08_T1_refined_inStrain","LEAP_META_09_T1_refined_inStrain","LEAP_META_10_T1_refined_inStrain",
		"LEAP_META_11_T1_refined_inStrain","LEAP_META_12_T1_refined_inStrain","LEAP_META_13_T1_refined_inStrain","LEAP_META_14_T1_refined_inStrain","LEAP_META_15_T1_refined_inStrain",
		"LEAP_META_16_T1_refined_inStrain","LEAP_META_17_T1_refined_inStrain","LEAP_META_18_T1_refined_inStrain",]

for profile in profile_list:
	IS = inStrain.SNVprofile.SNVprofile(profile)
	covT = IS.get('covT')
	all_scaf = {}
	scaffolds = open("good_contig_list.txt", "r")
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
