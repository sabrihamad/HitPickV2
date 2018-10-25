'''
This script is a modification from the original one writen by Dr. Ingo Vogt.
run using the following command:
python hitpick_parallel.py <input_file> <output_destination>


The database compounds is final_database.pickle saved in pickle format for
faster loading.

Query compounds should be in SMILES and tab seperated.

The script creates a folder names 'Results' in the output destination and creates
the output file in that folder. The output file is named exactly as the input file
but with '_output' appended to it.

Because of parallel computing, this script will take 1.2s per compound when submitting
1000 compounds.

Sabri Hamad
sabrihamad@outlook.com
'''
# Import packages
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit import RDLogger   #to suppress the warnings
from multiprocessing import Pool
import gzip
import time,sys
import os
import random
import csv
import math
import pickle
import numpy as np
import os.path
import json
import csv

# Definitions
mffp_radius = 3
mffp_nBits=9192
max_predictions = 5000
parallel=True
sessionIDLength=25
#minimumSimilarity=0.1  # old version
minimumSimilarity=0.05
maxQueryCompounds=500000
workers=7
max_predictions=10
draw_input=True #draw input molecules
draw_sim_comp=True #draw most similar compounds
#parameters:
cwd = os.getcwd()
db=cwd+'/databases/final_database.pickle'
modelPath=cwd+'/databases/bayesian_models'
sessionID=''
#path='results'
ensg2symbol=''

#def functions is to define the functions before running main to solve the problem of multiprocessing on the VM
def functions():

	def truncate(f, n):
		'''Truncates/pads a float f to n decimal places without rounding'''
		s = '{}'.format(f)
		if 'e' in s or 'E' in s:
			return '{0:.{1}f}'.format(f, n)
		i, p, d = s.partition('.')
		return '.'.join([i, (d+'0'*n)[:n]])

	def occur_cat(x):
		if 2<=x<4:
			y='occur: [2-4)'
			return y
		elif 4<=x<6:
			y='occur: [4-6)'
			return y
		elif 6<=x<8:
			y='occur: [6-8)'
			return y
		elif 8<=x<10:
			y='occur: [8-10)'
			return y
		elif x>=10:
			y='occur: >10'
			return y
		elif 1<=x<2:
			y='occur: 1'
			return y

	def read_input_smiles(filename):
		##Input format: <molecule id><tab><smiles> . Identical molecule names will not create separate entries, rather only contain the last read molecule.##

		#open and read the files:
		fIN = file(filename,'r')
		data = fIN.readlines()
		fIN.close()

		mols = {}
		c = 0		#counter
		for line in data:
			c += 1
			#vals = line.strip().split("\t")

			vals = line.strip().split() #Monica
			# if more than one column are provided, first is taken as name and second as smiles string, if only one column, this should be the smiles string, name will then be the molecule count
			if len(vals)>1:
				mol = Chem.MolFromSmiles(vals[1].strip(),sanitize=True)
			else:
				mol = Chem.MolFromSmiles(vals[0].strip(),sanitize=True)

			if mol is not None:
				if len(vals)>1:
					name = ''.join(vals[0])
				else:
					name = "query_"+str(c)

				#creat the mol dictionary: mol[name]:{mol,fp,onBits}
				mols[name] = {}
				mols[name]['smiles'] = vals[1]
				mols[name]['mol'] = mol	# molecule object
				mols[name]['fp'] = AllChem.GetMorganFingerprintAsBitVect(mol,mffp_radius,nBits=mffp_nBits,useFeatures=True)	# create fingerprint
				mols[name]['onBits']=list(mols[name]['fp'].GetOnBits())

		print ("read input: Done!")
		return mols

	def read_db_file():
		##DB file format
		with open(db,'rb') as pIN:
			database_compounds=pickle.load(pIN)
		return database_compounds
			#q.put([dbCompound[0],AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(dbCompound[1].strip(),sanitize=True),mffp_radius,nBits=mffp_nBits,useFeatures=True),dbCompound[2]])

	def search_10_target(query_all):
		querycounter=0
		query_molecules={}

		for query_comp in query_all:
			querycounter+=1
			query_molecules[query_comp]={}

			all_sim_compounds={}
			all_sim_compounds[query_comp]=[]
			flag=0
			dbcounter=0
			#print ("searching for 10 targets, query comp# ",querycounter, ' of ', len(query_all) #,"dbcompound# ", dbcounter)
			for db_comp in dbComps:
				dbcounter+=1

				tc_sim=DataStructs.TanimotoSimilarity(query_all[query_comp]['fp'],dbComps[db_comp][2])
				#this produces a list of tuples of (   database compound  ,  Tanimoto similarity   ,  targets )
				all_sim_compounds[query_comp].append((db_comp,tc_sim,dbComps[db_comp][3].split(";")))
				if tc_sim==float(1):

					if dbComps[db_comp][1]== Chem.rdinchi.MolToInchi(query_all[query_comp]['mol'])[0]:
						flag=1
						break

			#sort the tuples by tanimoto similarity

			if flag==1:
				query_molecules[query_comp]['mol']=query_all[query_comp]['mol']
				query_molecules[query_comp]["sim_comp"]=db_comp
				query_molecules[query_comp]['targets_sorted']=dbComps[db_comp][3].split(";")
				query_molecules[query_comp]['in_db']=True
			else:
				all_sim_compounds_sorted = sorted(all_sim_compounds[query_comp], key=lambda s:s[1],reverse=True)
				targets={}
				for sim_comp in all_sim_compounds_sorted:
					for target in sim_comp[2]:
						if target in targets:
							targets[target]["#occur"]+=1
						else:
							targets[target]={}
							targets[target]["#occur"]=1
							targets[target]["tc_sim"]=round(sim_comp[1],3)
							targets[target]["sim_comp"]=sim_comp[0]
					if len(targets)>=max_predictions:
						break
				query_molecules[query_comp]['mol']=query_all[query_comp]['mol']
				query_molecules[query_comp]['targets']=targets
				query_molecules[query_comp]['onBits']=query_all[query_comp]['onBits']
				query_molecules[query_comp]['scores']=[]
				query_molecules[query_comp]['targets_sorted']=[]
				query_molecules[query_comp]['in_db']=False

		return query_molecules

	def outputjson(query_file,output_path):
		  txtfile = ''.join([output_path.rstrip("/"),"/Results/",query_file.split("/")[-1].split(".")[0] , "_targets.txt"])
		  jsonfile = ''.join([output_path.rstrip("/"),"/Results/",query_file.split("/")[-1].split(".")[0], "_targets.json"])

		  # save the json output as emp.json

		  jsfile = file(jsonfile, 'w')
		  jsfile.write('{\r\n')
		  jsfile.write('    "data\": [\r\n')
		  with open(txtfile,'r') as f:
			# Monica : do not skip the heading
		#next(f) # skip headings
			  reader=csv.reader(f,delimiter='\t')

			  # get the total number of rows excluded the heading
			  row_count = len(list(reader))
			  ite = 0

			  # back to first position
			  f.seek(0)
			 # next(f) # skip headings # Monica skip headings

			  for a1,b1,c1,d1,e1,f1 in reader:

	  #            a2=''.join(['<img src="',output_path,'/Results/images/input/',a1,'.svg">'])

				  ite+= 1

				  jsfile.write('        [\r\n')

				  n = '            "' + a1 + '\",\r\n'
				# n = '            "' + a2 + '\",\r\n'
				  i = '            "' + b1 + '\",\r\n'
				  d = '            "' + c1 + '\",\r\n'
				  n2 = '            "' + d1 + '\",\r\n'
				  i2 = '            "' + e1 + '\",\r\n'
				  d2 = '            "' + f1 + '\"\r\n'

				  jsfile.write(n)
				  jsfile.write(i)
				  jsfile.write(d)
				  jsfile.write(n2)
				  jsfile.write(i2)
				  jsfile.write(d2)

				  jsfile.write('        ]')
				  # omit comma for last row item
				  if ite < row_count:
					   jsfile.write(',')

				  jsfile.write('\r\n')

		  jsfile.write('    ]\r\n}')
		  jsfile.close()

	def calculate_scores(query_mols):
		#print ('calc')
		#print (query_mols)
		for mol in query_mols:

			#print ("mol#: ",counter)
			if not query_mols[mol]['in_db']:
				original_scores = []
				for target in query_mols[mol]['targets']:
					# open and read the coresponding target file
					f = gzip.open(modelPath.rstrip("/")+'/bitScores_'+target+'_MFFP3.csv.gz', 'rb')
					model = f.read().split("\n")[1:]
					f.close()

					tmp = 0.0
					for idx in query_mols[mol]['onBits']:
						tmp += float(model[idx].split("\t")[3])
					original_scores.append((tmp,query_mols[mol]['targets'][target]['tc_sim'],target))  #score=[(score1,target1)  ,  (score2,target2)  ,  ]

				#the scores needs to be shifted so that the minimum is zero
				inc=min(original_scores)[0]*-1
				x=[]
				# x=[  (original_score  + increment  + 10* number of occurence  ,   target   )    ]
				for score in original_scores:
					x.append((score[0]+inc+10*query_mols[mol]['targets'][score[2]]['#occur'],score[1],score[2]))

				scores=x
				scores.sort(reverse=True)
				scores=scores[0:max_predictions]
				scores=sorted(scores, key=lambda scores:scores[1],reverse=True)

				top_scores=(mol,scores)		#top_scores=(    mol,    [ (score1,target1) , ...., (score10,target10)  ]    )
				query_mols[mol]['targets_sorted']=[x[2] for x in top_scores[1]]
				query_mols[mol]['scores']=[x[0] for x in top_scores[1]]
		return query_mols

	def molImage_to_directory(mol,name,directory='',width=300,height=300,filetype="svg"):
		'''Write molecule image to specified directory. Return name of created file'''
		Draw.MolToFile(mol,''.join([directory.rstrip("/"),"/",name,".",filetype]),size=(width, height), kekulize=True, wedgeBonds=True, title=name, imageType=filetype)
		return ''.join([directory.rstrip("/"),"/",name,".",filetype])

	def output(list_mols,query_file,output_path):

		if not os.path.exists(output_path+'/Results'):
			os.system("mkdir "+output_path+"/Results")
		if draw_input:
			if not os.path.exists(output_path+'/Results/images'):
				os.system("mkdir "+output_path+"/Results/images")
				os.system("mkdir "+output_path+"/Results/images/input")
				os.system("mkdir "+output_path+"/Results/images/db")

		with open(output_path.rstrip("/")+"/Results/"+query_file.split("/")[-1].split(".")[0] + "_targets.txt",'wb') as fOUT:
			for mols in list_mols:
				for mol in mols:
					if draw_input:
						# write input molecule image
						file1 = molImage_to_directory(mols[mol]['mol'],mol.replace(";","_"),directory=output_path.rstrip()+"/Results/images/input")


					#print (mols[mol]['targets_sorted'])
					for indx,target in enumerate(mols[mol]['targets_sorted']):
						if not mols[mol]['in_db']:
							if truncate(mols[mol]['targets'][target]['tc_sim'],1) not in precision_table[">%d actual targets" %(indx+1)][occur_cat(mols[mol]['targets'][target]['#occur'])]:
								precision_table[">%d actual targets" %(indx+1)][occur_cat(mols[mol]['targets'][target]['#occur'])][truncate(mols[mol]['targets'][target]['tc_sim'],1)]=np.nan

							fOUT.write(  "\t".join(   [mol,  comp_name[mols[mol]['targets'][target]['sim_comp']], target,  str(mols[mol]['targets'][target]['tc_sim']) ,str(precision_table[">%d actual targets" %(indx+1)][occur_cat(mols[mol]['targets'][target]['#occur'])][truncate(mols[mol]['targets'][target]['tc_sim'],1)]) ,str(mols[mol]['targets'][target]['#occur']) ]     )   + "\n"   )
						else:
							fOUT.write(  "\t".join(   [mol,   comp_name[mols[mol]['sim_comp']] , target,  str(1.0) , str(100),str(np.nan) ]     )   + "\n"   )
						if draw_sim_comp:
							# write most similar db molecule image

							try:
								sim_comp_name=(comp_name[mols[mol]['targets'][target]['sim_comp']][:30]) if len(mols[mol]['targets'][target]['sim_comp']) > 30 else mols[mol]['targets'][target]['sim_comp']
								file2 = molImage_to_directory(Chem.MolFromSmiles(dbComps[mols[mol]['targets'][target]['sim_comp']][0].strip(),sanitize=True),comp_name[mols[mol]['targets'][target]['sim_comp']],directory=output_path.rstrip()+"/Results/images/db")
							except:
								sim_comp_name=(comp_name[mols[mol]['sim_comp']][:30]) if len(mols[mol]['sim_comp']) > 30 else comp_name[mols[mol]['sim_comp']]
								file2 = molImage_to_directory(Chem.MolFromSmiles(dbComps[mols[mol]['sim_comp']][0].strip(),sanitize=True),comp_name[mols[mol]['sim_comp']],directory=output_path.rstrip()+"/Results/images/db")
		fOUT.close()



		#with open(path+"targetPredictionall.txt",'wb') as fOUT2:
			#for mol in mols:
				#fOUT2.write(  "\t".join(   [mol, ";".join([mols[mol]['targets'][target]['sim_comp'] for target in mols[mol]['targets']]), ";".join([target for target in mols[mol]['targets']]),  ";".join([str(mols[mol]['targets'][target]['tc_sim']) for target in mols[mol]['targets']])  ,";".join([str(mols[mol]['targets'][target]['#occur']) for target in mols[mol]['targets']]) ]     )   + "\n"   )

	def dictionary():
		'''this function matches the compound id (SBSMID) to its name in the dictionary comp_name.pickle'''
		global gene_name
		gene_name={}
		with open('/home/sabri/Dropbox/hitpick/databases/Genes.csv','rb') as fIN:
			data=csv.reader(fIN,delimiter="\t")
			for row in data:
				intrez=row[0]
				gene_name[row[0]]=row[4]
		fIN.close()

		#comp_names.pickle matches the SBSMID to the compound's name or ID
		global comp_name
		comp_name={}
		with open('/home/sabri/Dropbox/hitpick/databases/comp_names.pickle','rb') as pIN:
			comp_name=pickle.load(pIN)

	def precision():
		'''this function loads the pricision matrix which stored as a dictionary'''
		global precision_table
		with open('/home/sabri/Dropbox/hitpick/databases/precision_table.pickle','rb') as pIN1:
			precision_table=pickle.load(pIN1)

def truncate(f, n):
	'''Truncates/pads a float f to n decimal places without rounding'''
	s = '{}'.format(f)
	if 'e' in s or 'E' in s:
		return '{0:.{1}f}'.format(f, n)
	i, p, d = s.partition('.')
	return '.'.join([i, (d+'0'*n)[:n]])

def occur_cat(x):
	if 2<=x<4:
		y='occur: [2-4)'
		return y
	elif 4<=x<6:
		y='occur: [4-6)'
		return y
	elif 6<=x<8:
		y='occur: [6-8)'
		return y
	elif 8<=x<10:
		y='occur: [8-10)'
		return y
	elif x>=10:
		y='occur: >10'
		return y
	elif 1<=x<2:
		y='occur: 1'
		return y

def read_input_smiles(filename):
	##Input format: <molecule id><tab><smiles> . Identical molecule names will not create separate entries, rather only contain the last read molecule.##

	#open and read the files:
	fIN = file(filename,'r')
	data = fIN.readlines()
	fIN.close()

	mols = {}
	c = 0		#counter
	for line in data:
		c += 1
		#vals = line.strip().split("\t")

		vals = line.strip().split() #Monica
		# if more than one column are provided, first is taken as name and second as smiles string, if only one column, this should be the smiles string, name will then be the molecule count
		if len(vals)>1:
			mol = Chem.MolFromSmiles(vals[1].strip(),sanitize=True)
		else:
			mol = Chem.MolFromSmiles(vals[0].strip(),sanitize=True)

		if mol is not None:
			if len(vals)>1:
				name = ''.join(vals[0])
			else:
				name = "query_"+str(c)

			#creat the mol dictionary: mol[name]:{mol,fp,onBits}
			mols[name] = {}
			mols[name]['smiles'] = vals[1]
			mols[name]['mol'] = mol	# molecule object
			mols[name]['fp'] = AllChem.GetMorganFingerprintAsBitVect(mol,mffp_radius,nBits=mffp_nBits,useFeatures=True)	# create fingerprint
			mols[name]['onBits']=list(mols[name]['fp'].GetOnBits())

	print ("read input: Done!")
	return mols

def read_db_file():
	##DB file format
	with open(db,'rb') as pIN:
		database_compounds=pickle.load(pIN)
	return database_compounds
		#q.put([dbCompound[0],AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(dbCompound[1].strip(),sanitize=True),mffp_radius,nBits=mffp_nBits,useFeatures=True),dbCompound[2]])

def search_10_target(query_all):
	querycounter=0
	query_molecules={}

	for query_comp in query_all:
		querycounter+=1
		query_molecules[query_comp]={}
		query_molecules[query_comp]['in_db']=(False,False)
		query_molecules[query_comp]['indb_targets_sorted']=[]
		all_sim_compounds={}
		all_sim_compounds[query_comp]=[]
		flag=0
		dbcounter=0
		#print ("searching for 10 targets, query comp# ",querycounter, ' of ', len(query_all) #,"dbcompound# ", dbcounter)
		for db_comp in dbComps:
			dbcounter+=1

			tc_sim=DataStructs.TanimotoSimilarity(query_all[query_comp]['fp'],dbComps[db_comp][2])
			if tc_sim==float(1):

				if (dbComps[db_comp][1]== Chem.rdinchi.MolToInchi(query_all[query_comp]['mol'])[0]) | (Chem.MolToSmiles(Chem.MolFromSmiles(dbComps[db_comp][0]))== Chem.MolToSmiles(query_all[query_comp]['mol'])) :
					if len(dbComps[db_comp][3].split(";"))>=10:
						query_molecules[query_comp]['mol']=query_all[query_comp]['mol']
						query_molecules[query_comp]["indb_sim_comp"]=db_comp
						query_molecules[query_comp]['indb_targets_sorted']=dbComps[db_comp][3].split(";")
						query_molecules[query_comp]['in_db']=(True,True)
						break
					else:
						#query_molecules[query_comp]['mol']=query_all[query_comp]['mol']
						query_molecules[query_comp]["indb_sim_comp"]=db_comp
						query_molecules[query_comp]['indb_targets_sorted']=dbComps[db_comp][3].split(";")
						query_molecules[query_comp]['in_db']=(True,False)
						continue
			#this produces a list of tuples of (   database compound  ,  Tanimoto similarity   ,  targets )
			all_sim_compounds[query_comp].append((db_comp,tc_sim,dbComps[db_comp][3].split(";")))

		#sort the tuples by tanimoto similarity
		if not query_molecules[query_comp]['in_db'][1]:
			all_sim_compounds_sorted = sorted(all_sim_compounds[query_comp], key=lambda s:s[1],reverse=True)
			targets={}
			for sim_comp in all_sim_compounds_sorted:
				for target in sim_comp[2]:
					if target not in query_molecules[query_comp]['indb_targets_sorted']:
						if target in targets:
							targets[target]["#occur"]+=1
						else:
							targets[target]={}
							targets[target]["#occur"]=1
							targets[target]["tc_sim"]=round(sim_comp[1],3)
							targets[target]["sim_comp"]=sim_comp[0]
				if query_molecules[query_comp]['in_db'][0]:
					if (len(query_molecules[query_comp]['indb_targets_sorted'] ) + len(targets))  >= max_predictions:
						break
				if (len(targets)>=max_predictions)  :
					break
			query_molecules[query_comp]['mol']=query_all[query_comp]['mol']
			query_molecules[query_comp]['targets']=targets
			query_molecules[query_comp]['onBits']=query_all[query_comp]['onBits']
			query_molecules[query_comp]['scores']=[]
			query_molecules[query_comp]['targets_sorted']=[]


	return query_molecules

def outputjson(query_file,output_path):
	  txtfile = ''.join([output_path.rstrip("/"),"/Results/",query_file.split("/")[-1].split(".")[0] , "_targets.txt"])
	  jsonfile = ''.join([output_path.rstrip("/"),"/Results/",query_file.split("/")[-1].split(".")[0], "_targets.json"])

	  # save the json output as emp.json

	  jsfile = file(jsonfile, 'w')
	  jsfile.write('{\r\n')
	  jsfile.write('    "data\": [\r\n')
	  with open(txtfile,'r') as f:
		# Monica : do not skip the heading
	#next(f) # skip headings
		  reader=csv.reader(f,delimiter='\t')

		  # get the total number of rows excluded the heading
		  row_count = len(list(reader))
		  ite = 0

		  # back to first position
		  f.seek(0)
		 # next(f) # skip headings # Monica skip headings

		  for a1,b1,c1,d1,e1,f1 in reader:

  #            a2=''.join(['<img src="',output_path,'/Results/images/input/',a1,'.svg">'])

			  ite+= 1

			  jsfile.write('        [\r\n')

			  n = '            "' + a1 + '\",\r\n'
			# n = '            "' + a2 + '\",\r\n'
			  i = '            "' + b1 + '\",\r\n'
			  d = '            "' + c1 + '\",\r\n'
			  n2 = '            "' + d1 + '\",\r\n'
			  i2 = '            "' + e1 + '\",\r\n'
			  d2 = '            "' + f1 + '\"\r\n'

			  jsfile.write(n)
			  jsfile.write(i)
			  jsfile.write(d)
			  jsfile.write(n2)
			  jsfile.write(i2)
			  jsfile.write(d2)

			  jsfile.write('        ]')
			  # omit comma for last row item
			  if ite < row_count:
				   jsfile.write(',')

			  jsfile.write('\r\n')

	  jsfile.write('    ]\r\n}')
	  jsfile.close()

def calculate_scores(query_mols):
	#print ('calc')
	#print (query_mols)
	for mol in query_mols:

		#print ("mol#: ",counter)
		if not query_mols[mol]['in_db'][1]:
			original_scores = []
			for target in query_mols[mol]['targets']:
				# open and read the coresponding target file
				f = gzip.open(modelPath.rstrip("/")+'/bitScores_'+target+'_MFFP3.csv.gz', 'rb')
				model = f.read().split("\n")[1:]
				f.close()

				tmp = 0.0
				for idx in query_mols[mol]['onBits']:
					tmp += float(model[idx].split("\t")[3])
				original_scores.append((tmp,query_mols[mol]['targets'][target]['tc_sim'],target))  #score=[(score1,target1)  ,  (score2,target2)  ,  ]

			#the scores needs to be shifted so that the minimum is zero
			inc=min(original_scores)[0]*-1
			x=[]
			# x=[  (original_score  + increment  + 10* number of occurence  ,   target   )    ]
			for score in original_scores:
				x.append((score[0]+inc+10*query_mols[mol]['targets'][score[2]]['#occur'],score[1],score[2]))

			scores=x
			scores.sort(reverse=True)
			scores=scores[0:max_predictions]
			scores=sorted(scores, key=lambda scores:scores[1],reverse=True)

			top_scores=(mol,scores)		#top_scores=(    mol,    [ (score1,target1) , ...., (score10,target10)  ]    )
			query_mols[mol]['targets_sorted']=[x[2] for x in top_scores[1]]
			query_mols[mol]['scores']=[x[0] for x in top_scores[1]]
	return query_mols

def molImage_to_directory(mol,name,directory='',width=300,height=300,filetype="svg"):
	'''Write molecule image to specified directory. Return name of created file'''
	Draw.MolToFile(mol,''.join([directory.rstrip("/"),"/",name,".",filetype]),size=(width, height), kekulize=True, wedgeBonds=True, title=name, imageType=filetype)
	return ''.join([directory.rstrip("/"),"/",name,".",filetype])

def output(list_mols,query_file,output_path):

	if not os.path.exists(output_path+'/Results'):
		os.system("mkdir "+output_path+"/Results")
	if draw_input:
		if not os.path.exists(output_path+'/Results/images'):
			os.system("mkdir "+output_path+"/Results/images")
			os.system("mkdir "+output_path+"/Results/images/input")
			os.system("mkdir "+output_path+"/Results/images/db")

	with open(output_path.rstrip("/")+"/Results/"+query_file.split("/")[-1].split(".")[0] + "_targets.txt",'wb') as fOUT:
		for mols in list_mols:
			for mol in mols:
				targets_perfect_match=[]
				if draw_input:
					# write input molecule image
					file1 = molImage_to_directory(mols[mol]['mol'],mol.replace(";","_"),directory=output_path.rstrip()+"/Results/images/input")
				#print (mols[mol]['targets_sorted'])
				if mols[mol]['in_db'][0]:
					targets_perfect_match=mols[mol]['indb_targets_sorted']
					for indx,target in enumerate(targets_perfect_match):
						fOUT.write(  "\t".join(   [mol,   comp_name[mols[mol]['indb_sim_comp']] , target,  str(1.0) , str(100),str(np.nan) ]     )   + "\n"   )
						if draw_sim_comp:
							# write most similar db molecule image
							sim_comp_name=(comp_name[mols[mol]['indb_sim_comp']][:30] if len(mols[mol]['indb_sim_comp']) > 30 else comp_name[mols[mol]['indb_sim_comp']])
							file2 = molImage_to_directory(Chem.MolFromSmiles(dbComps[mols[mol]['indb_sim_comp']][0].strip(),sanitize=True),comp_name[mols[mol]['indb_sim_comp']].replace("/","-"),directory=output_path.rstrip()+"/Results/images/db")
				if not mols[mol]['in_db'][1]:
					for indx,target in enumerate(mols[mol]['targets_sorted']):

						if indx+len(targets_perfect_match)>9:
							break
						if truncate(mols[mol]['targets'][target]['tc_sim'],1) not in precision_table[">%d actual targets" %(indx+1)][occur_cat(mols[mol]['targets'][target]['#occur'])]:
							precision_table[">%d actual targets" %(indx+1)][occur_cat(mols[mol]['targets'][target]['#occur'])][truncate(mols[mol]['targets'][target]['tc_sim'],1)]=np.nan
						fOUT.write(  "\t".join(   [mol,  comp_name[mols[mol]['targets'][target]['sim_comp']], target,  str(mols[mol]['targets'][target]['tc_sim']) ,str(precision_table[">%d actual targets" %(indx+1)][occur_cat(mols[mol]['targets'][target]['#occur'])][truncate(mols[mol]['targets'][target]['tc_sim'],1)]) ,str(mols[mol]['targets'][target]['#occur']) ]     )   + "\n"   )
						if draw_sim_comp:
							try:
								sim_comp_name=(comp_name[mols[mol]['targets'][target]['sim_comp']][:30]) if len(mols[mol]['targets'][target]['sim_comp']) > 30 else mols[mol]['targets'][target]['sim_comp']
								file2 = molImage_to_directory(Chem.MolFromSmiles(dbComps[mols[mol]['targets'][target]['sim_comp']][0].strip(),sanitize=True),comp_name[mols[mol]['targets'][target]['sim_comp']].replace("/","-"),directory=output_path.rstrip()+"/Results/images/db")
							except:
								sim_comp_name=(comp_name[mols[mol]['sim_comp']][:30].replace("/","-")) if len(mols[mol]['sim_comp']) > 30 else comp_name[mols[mol]['sim_comp']]
								file2 = molImage_to_directory(Chem.MolFromSmiles(dbComps[mols[mol]['sim_comp']][0].strip(),sanitize=True),comp_name[mols[mol]['sim_comp']].replace("/","-"),directory=output_path.rstrip()+"/Results/images/db")


		fOUT.close()



	#with open(path+"targetPredictionall.txt",'wb') as fOUT2:
		#for mol in mols:
			#fOUT2.write(  "\t".join(   [mol, ";".join([mols[mol]['targets'][target]['sim_comp'] for target in mols[mol]['targets']]), ";".join([target for target in mols[mol]['targets']]),  ";".join([str(mols[mol]['targets'][target]['tc_sim']) for target in mols[mol]['targets']])  ,";".join([str(mols[mol]['targets'][target]['#occur']) for target in mols[mol]['targets']]) ]     )   + "\n"   )

def dictionary():
	'''this function matches the compound id (SBSMID) to its name in the dictionary comp_name.pickle'''
	global gene_name
	gene_name={}
	with open(cwd+'/databases/Genes.csv','rb') as fIN:
		data=csv.reader(fIN,delimiter="\t")
		for row in data:
			intrez=row[0]
			gene_name[row[0]]=row[4]
	fIN.close()

	#comp_names.pickle matches the SBSMID to the compound's name or ID
	global comp_name
	comp_name={}
	with open(cwd+'/databases/comp_names.pickle','rb') as pIN:
		comp_name=pickle.load(pIN)

def precision():
	'''this function loads the pricision matrix which stored as a dictionary'''
	global precision_table
	with open(cwd+'/databases/precision_table.pickle','rb') as pIN1:
		precision_table=pickle.load(pIN1)

def main(parameters):
	start=time.time()
	#suppress warnings
	lg=RDLogger.logger()
	lg.setLevel(RDLogger.ERROR)

	query_file=parameters[0] #location of input file
	output_path=parameters[1] #path for output file

	#load the drug name dictionary and the precision matrix
	dictionary()
	precision()

	# STEP #1: read the query compounds
	print ('reading input file...........')
	input_molecules=read_input_smiles(query_file.lstrip("-"))
	# STEP #2: load the database compounds
	print ('loading database file........')
	global dbComps
	dbComps=read_db_file()
	print ('Done!')

	# STEP #3: search for 10 targets
	print ('searching for 10 targets........')
	###divide the input molecules into a list of dictionaries to be processed in parallel###
	input_molecules_list=[{k:input_molecules[k]} for k in input_molecules.keys()]
	pool = Pool(processes=workers)
	target_dict=pool.map(search_10_target, input_molecules_list)
	pool.close()
	pool.join()
	print ("search: Done!")

	# STEP #4: calculate and sort the scores
	print ('calculating scores.........')

	pool = Pool(processes=workers)
	input_molecules=pool.map(calculate_scores, target_dict)
	pool.close()
	pool.join()
	print ("scores calculation: Done!")
	# STEP #5: output a txt file


	print ('creating output directory and file.......')
	output(input_molecules,query_file,output_path)
	##file in json format
	outputjson (query_file,output_path)
	print ("output: Done")

	usedTime=time.time() -start

	if len(input_molecules)>0:
		print ("\nTime elapsed: "+str(usedTime)+" seconds, "+str(usedTime/float(len(input_molecules)))+" seconds per input molecule.\n")
		#print ("\nTime elapsed to read the database: "+ str(x_time))
	else:
		print ("\nTime elapsed: "+str(usedTime)+" seconds, 0 seconds per input molecule.\n")


if __name__ == '__main__':

	functions()

	main(sys.argv[1:])
