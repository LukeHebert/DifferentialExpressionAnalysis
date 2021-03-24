'''
Author: Luke Hebert
Date begun: February 20th, 2021?
Description:
	take expression profile results (log-fold change data from experiments)
		and {transcripts:genenames} pickled dictionary
		and more raw .csv files that include padj values
	create variations of bar graphs to depict data
	output each graph version style in a separate folder
'''

#import necessary libraries
import seaborn as sns
import matplotlib.pyplot as plt
import sys, os, pickle

#assign a slash symbol based on the current machine's operating system
slash = '\\' if os.name == 'nt' else '/'
#store input file pathway and also save input folder
in_folder_path = sys.argv[1]
in_files_list = []
for root,dirs,files in os.walk(in_folder_path):
	in_files_list.extend(files)
	break

#load the pickled dictionary of {'transcripts':'genenames'}
for file_name in in_files_list:
	if '.pkl' in file_name:
		gene_dict = pickle.load(open(file_name,"rb"))


#load adjusted p values from the "padj" files, filtered by padj .05 & log2fold 2/-2
padj_dict = {}
for file_name in in_files_list:
	if 'padj' in file_name:
		with open(in_folder_path+slash+file_name) as in_file:
			for i,line in enumerate(in_file):
				line_list = line.replace('\r','').replace('\n','').split('\t')
				if i == 0:
					pass
				elif line_list[2]!='NA' and line_list[6]!='NA':
					transcript = line_list[0].replace('\"','').replace(',',';')
					log2fold = float(line_list[2])
					padj = float(line_list[6])
					#line_list2 is the logbase2foldchange & line_list6 is adjusted p value
					if (log2fold >= 2 or log2fold <= -2) and padj < 0.05:
						padj_dict[gene_dict[transcript]]=padj


#read in logfold2 data to a dictionary for making graphs
out_dict = {}
for file_name in in_files_list:
	if 'Fig' in file_name:
		out_dict[file_name[:-4]] = {'gene':[],'log2float':[],'label':[], 'padj':[]}
		with open(in_folder_path+slash+file_name) as in_file:
			for line in in_file:
				line_list = line.replace('\r','').replace('\n','').split(',')
				out_dict[file_name[:-4]]['gene'].append('$\it{'+line_list[0]+'}$')
				out_dict[file_name[:-4]]['log2float'].append(float(line_list[1]))
				out_dict[file_name[:-4]]['label'].append(line_list[2].replace('WNT','$\it{WNT}$'))
				out_dict[file_name[:-4]]['padj'].append("{:.2e}".format( float(padj_dict[line_list[0]]) ))
					#formatting tip from here: https://stackoverflow.com/questions/6913532/display-a-decimal-in-scientific-notation



#for each input file, cycle through the dictionary and call the graphing function for each item (and output)
for graph in out_dict:

	'''
	#sort the graphing dictionary by foldchange
	zipped_list = list( zip(out_dict[graph]['gene'], out_dict[graph]['log2float'], out_dict[graph]['label'], out_dict[graph]['padj']) )
	zipped_list.sort(key=lambda x:x[1]) #sort based on log2float
	out_dict[graph]['gene'] = [x for x, this, that, other in zipped_list]
	out_dict[graph]['label'] = [x for this, that, x, other in zipped_list]
	out_dict[graph]['padj'] = [x for this, that, other, x  in zipped_list]
	'''

	out_folder=in_folder_path+slash+graph+slash
	if not os.path.exists(out_folder):
		os.makedirs(out_folder)
		print('\nOutputting:\t' + str(out_folder))###
	else:
		print('This folder already exists and was therefore not overwritten: '+out_folder)


	#make a list of all combination tuples, which will be iterated over to create
		#graphs of varying properties; tuples should be of the format:
		#(gridcolor, palette name, orientation, dodge boolean)
	graph_versions = (
	('darkgrid','Accent','h',False),
	('darkgrid','Set1_r','h',False),
	('darkgrid','Set3','h',False),
	('darkgrid','colorblind','h',False),
	('darkgrid','pastel','h',False)
	)
	#loop through all graph version combinations
		#create, save, and clear the figure for each version
	for details_tup in graph_versions:
		gridcolor, palette_name, orientation, dodge_bool = details_tup
		if len(set(out_dict[graph]['gene']))>25:
			plt.figure(figsize=(8,10))
		sns.set_theme(style=gridcolor,palette=palette_name)
		#avoid repeat y axis value issue by manually setting y tick labels and
			#using the indexes of x values for y values (stored in range_list)
		range_list=[]
		for i,item in enumerate(out_dict[graph]['log2float']):
			range_list.append(i)
		ax=sns.barplot(x=out_dict[graph]['log2float'],
						y=range_list,
						dodge=dodge_bool,#avoid tiny bar widths
						hue=out_dict[graph]['label'],
						orient=orientation)
		ax.set_yticklabels(out_dict[graph]['gene'])
		#add padj values near the bars
		for i in range(len(out_dict[graph]['gene'])):
			if out_dict[graph]['log2float'][i] > 0:
				horiz_adj = 'left'
			else:
				horiz_adj = 'right'

			ax.text(out_dict[graph]['log2float'][i], i+0.25, str(out_dict[graph]['padj'][i]),
					fontdict=dict(color='red',fontsize=6),
					ha=horiz_adj, va='baseline')
		#add gene names inside the bars
		for i in range(len(out_dict[graph]['gene'])):
			ax.text(out_dict[graph]['log2float'][i]/2, i+0.25, str(out_dict[graph]['gene'][i]),
					fontdict=dict(color='black',fontsize=6),
					ha='center', va='baseline')
		ax.figure.savefig(''.join([out_folder,slash,graph,'_',orientation,'_',gridcolor[0:4],'_',palette_name,'.png']),dpi=800)
		plt.close() #clear the current figure to make room for next one


#keeping this code below in case I ever need to switch back to "vertical" barplot orientation
'''
	#horizontal, light =========================================================
	if len(set(out_dict[graph]['gene']))>25:
		plt.figure(figsize=(8,5))
	sns.set_theme(style="whitegrid",palette='colorblind')
	ax=sns.barplot(x=out_dict[graph]['gene'],
					y=out_dict[graph]['log2float'],
					dodge=boolboi,#avoid tiny bar widths
					hue=out_dict[graph]['label'],
					orient='v')
	if len(set(out_dict[graph]['gene']))>25:
		plt.setp(ax.get_xticklabels(), rotation=-45)
	else:
		plt.setp(ax.get_xticklabels(), rotation=-30)
	ax.figure.savefig(''.join([out_folder,slash,graph,'_v','_light.png']),dpi=800)
	plt.close() #clear the current figure to make room for next one
'''
