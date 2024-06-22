# -*- coding: utf-8 -*-
"""windowing-method.ipynb
"""

#How to run the pipeline
#python3 pipeline_for_frustration_clean.py patho_to_pdbs path_to_results #frames #protein_length #sim_start
#python3 pipeline_for_frustration_clean.py /home/titanx1/Documents/sims/frustra/rfah/traj1/ /home/titanx1/Documents/sims/frustra/rfah/traj1/Results/ 162 0

#sys.argv[1] -> path to pdbs files
#sys.argv[2] -> path where do you want to save the frustration results
#sys.argv[3] -> protein length
#sys.argv[4] -> start frame

#import libraries
import os
import numpy as np
import sys

#Preparing frustration calculations:
#The R scritp is generated to run the frustrations of the trajectory frames

out_r=open('r_frustration.R','w')

os.system('cd '+sys.argv[1]+';ls *.pdb | wc -l > aux')
aux=open('aux')
laux=aux.readline()
n_pdbs=int(laux.rstrip('\n'))
os.system('rm '+sys.argv[1]+'aux')

out_r.write('library(reticulate)\n')
out_r.write('use_python("/usr/bin/python3")\n')
out_r.write('Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")\n')
out_r.write('reticulate::py_config()\n')
out_r.write('library(frustratometeR)\n')
out_r.write('PdbsDir <- "'+sys.argv[1]+'"\n')
out_r.write('ResultsDir <- "'+sys.argv[2]+'"\n')
out_r.write('OrderList <-c()\n')
out_r.write('for(i in as.numeric('+sys.argv[4]+'):as.numeric('+str(n_pdbs)+')){OrderList <- c(OrderList, paste("pdb",i,".pdb",sep=""))}\n')
out_r.write('Dynamic_mutational <- dynamic_frustration(PdbsDir = PdbsDir, ResultsDir = ResultsDir,\n')
out_r.write('                                    GIFs = FALSE, Mode = "mutational")\n')
out_r.write('Dynamic_configurational <- dynamic_frustration(PdbsDir = PdbsDir, ResultsDir = ResultsDir,\n')
out_r.write('                                    GIFs = FALSE, Mode = "configurational")\n')
out_r.close()
print('Start frustration calculations')
os.system('Rscript r_frustration.R > Frustration')
print('End')

#Windowing method

l_protein=int(sys.argv[3])
sim_start=int(sys.argv[4])

out_ref=open(sys.argv[1]+'Reference.csv','w')
out_ref.write('Res Min Max Neu CMin CMax CNeu\n')


tam_vent=int(float(n_pdbs)*0.05) # Here we calculate the window size
ref_min=np.zeros(l_protein+1)
ref_neu=np.zeros(l_protein+1)
ref_max=np.zeros(l_protein+1)
mode=['configurational','mutational']
folder_name=sys.argv[1]+'/pngs-all/'

os.system('mkdir '+folder_name)
residues=[]
print('Start SD and mean calculations')
for x in range(0, len(mode)):
	for i in range(sim_start, int(laux),int(tam_vent)):
	  cmin=np.zeros(l_protein+1)
	  cmax=np.zeros(l_protein+1)
	  cneu=np.zeros(l_protein+1)
	  for j in range(i, int(tam_vent)+i):
	     if os.path.exists(sys.argv[2]+'pdb'+str(j)+'.done/FrustrationData/pdb'+str(j)+'.pdb_'+mode[x]+'_5adens'):
	       fst=open(sys.argv[2]+'pdb'+str(j)+'.done/FrustrationData/pdb'+str(j)+'.pdb_'+mode[x]+'_5adens')
	       for line in fst.readlines():
	        if not 'nMinimallyFrst' in line:
	          sp=line.rstrip('\n').split()
	         # nHighlyFrst nNeutrallyFrst nMinimallyFrst 3, 4 ,5
	          cmin[int(sp[0])]+= int(sp[5]) #here we add up the number of Highly frustrated contacts around a sphere of 5 ams for each residue within the window
	          cmax[int(sp[0])]+= int(sp[3]) #here we add up the number of Neutral frustrated contacts around a sphere of 5 ams for each residue within the window
	          cneu[int(sp[0])]+= int(sp[4]) #here we add up the number of Minimally frustrated contacts around a sphere of 5 ams for each residue within the window

	       fst.close()

	  if i == sim_start: #if the windows y the first one, we save this parameters because we need it for the future comparations (initia window)
	    for k in range(sim_start, l_protein+1):
	      ref_neu[k] = cneu[k]
	      ref_min[k] = cmin[k]
	      ref_max[k] = cmax[k]
	  else: # if te windows is not the first, we compare this with the initial window
	    for k in range(1, l_protein+1):
		if ref_min[k]!= 0:
		      div_ref_min=(cmin[k])/ref_min[k] #here we calculate th ratio between the Wn and the W0, for the minimally frustrated residues
		else:
			div_ref_min = 0     
		if ref_max[k]!= 0:
			div_ref_max=(cmax[k])/ref_max[k] # for highly frustrated
		else:
			div_ref_max=0
		if ref_neu[k]!= 0:
			div_ref_neu=(cneu[k])/ref_neu[k] # for neutral
		else:
			div_ref_neu=0
		out_ref.write(str(k)+' '+str(div_ref_min)+' '+str(div_ref_max)+' '+str(div_ref_neu)+' '+str(cmin[k]/tam_vent)+' '+str(cmax[k]/tam_vent)+' '+str(cneu[k]/tam_vent)+'\n')


	import pandas as pd

	df = pd.read_csv(sys.argv[1]+'Reference.csv',sep=' ')

	std_dev_grouped = df.groupby('Res')['Min'].std()
	mean_ref_grouped = df.groupby('Res')['Min'].mean()
	median_dev_grouped = df.groupby('Res')['CMin'].median()

        residues_filterd=[]
	for i in range(0, len(std_dev_grouped)):
		if mean_ref_grouped.iloc[i] !=0:
			if std_dev_grouped.iloc[i]/mean_ref_grouped.iloc[i] > 0.4 and float(median_dev_grouped.iloc[i]) > 4: #here we calculates the coefficient of variation if this number is above to 0.4 and the mean of the contacts is higher than 4 the residue change their frustration along the MD
				print(i+1,median_dev_grouped.iloc[i], std_dev_grouped.iloc[i])
				if not str((i+1)) in residues_filterd:
					 residues_filterd.append(str((i+1)))

#Frustration plots for selected residues using the filters:
#A scritp the R is generated to generate the output graphs
out_r=open('r_residues.R','w')

out_r.write('library(reticulate)\n')
out_r.write('use_python("/usr/bin/python3")\n')
out_r.write('Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")\n')
out_r.write('reticulate::py_config()\n')
out_r.write('library(frustratometeR)\n')
out_r.write('PdbsDir <- "'+sys.argv[1]+'"\n')
out_r.write('ResultsDir <- "'+sys.argv[2]+'"\n')
out_r.write('OrderList <-c()\n')
out_r.write('for(i in as.numeric('+sys.argv[4]+'):as.numeric('+str(n_pdbs)+')){OrderList <- c(OrderList, paste("pdb",i,".pdb",sep=""))}\n')
out_r.write('Dynamic_sing <- dynamic_frustration(PdbsDir = PdbsDir, ResultsDir = ResultsDir, OrderList = OrderList,\n')
out_r.write('                                    GIFs = FALSE, Mode = "singleresidue")\n')
for j in range(0,len(residues_filterd)):
   out_r.write('Dynamic_sing <- dynamic_res(Dynamic = Dynamic_sing, Resno = '+residues_filterd[j]+', Chain = "X", Graphics = TRUE)\n')
out_r.write('Dynamic_mutational <- dynamic_frustration(PdbsDir = PdbsDir, ResultsDir = ResultsDir, OrderList = OrderList,\n')
out_r.write('                                    GIFs = FALSE, Mode = "mutational")\n')
for j in range(0,len(residues_filterd)):
   out_r.write('Dynamic_mutational <- dynamic_res(Dynamic = Dynamic_mutational, Resno = '+residues_filterd[j]+', Chain = "X", Graphics = TRUE)\n')
out_r.write('Dynamic_configurational <- dynamic_frustration(PdbsDir = PdbsDir, ResultsDir = ResultsDir, OrderList = OrderList,\n')
out_r.write('                                    GIFs = FALSE, Mode = "configurational")\n')
for j in range(0,len(residues_filterd)):
   out_r.write('Dynamic_configurational <- dynamic_res(Dynamic = Dynamic_configurational, Resno = '+residues_filterd[j]+', Chain = "X", Graphics = TRUE)\n')

out_r.close()

os.system('Rscript r_residues.R')
for j in range(0,len(residues_filterd)):
   os.system('cp '+sys.argv[2]+'/Dynamic_plots_res_'+residues[j]+'_X/dynamic5adens_mutational_Res'+residues[j]+'.png '+sys.argv[1]+'/pngs-all/')
   os.system('cp '+sys.argv[2]+'/Dynamic_plots_res_'+residues[j]+'_X/dynamic5adens_configurational_Res'+residues[j]+'.png '+sys.argv[1]+'/pngs-all/')
   os.system('cp '+sys.argv[2]+'/Dynamic_plots_res_'+residues[j]+'_X/dynamic_IndexFrustration_singleresidue_Res'+residues[j]+'.png '+sys.argv[1]+'/pngs-all/')
