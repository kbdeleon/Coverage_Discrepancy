from __future__ import division
import sys
import subprocess
import os
from os.path import isfile
from os.path import join as jn
from os import listdir
import time
from numpy import mean
from numpy import median
import numpy
from scipy.stats import kurtosis,skew  #2018-10-24 scipy.stats throwing error with conda


def fmt_data(loci,cov):
	ct=0
	scaffolds_loci=[]
	scaffolds_cov=[]
	temp_loci=[]
	temp_cov=[]
	loci_trim=[]
	cov_trim=[]
	scaffold_avg_cov=[]
	scaffold_stdev=[]
	length=[]
	num_scaffolds=[]
	count_0_cov=[]
	for i in range(len(loci)):
		if loci[i]>ct:
			temp_loci.append(loci[i])
			temp_cov.append(cov[i])
			ct=loci[i]
		else:  #new scaffold
			scaffolds_loci.append(temp_loci)
			scaffolds_cov.append(temp_cov)
			temp_loci=[loci[i]]
			temp_cov=[cov[i]]
			ct=loci[i]
	scaffolds_loci.append(temp_loci)
	scaffolds_cov.append(temp_cov)
	num_scaffolds=len(scaffolds_loci)
	print "Number of Scaffolds:",  len(scaffolds_loci)
	for i in range(len(scaffolds_loci)):
		length.append(len(scaffolds_loci[i]))
		if len(scaffolds_loci[i])<401:
			print "Scaffold <400 bp"
		else:
			l_trim=scaffolds_loci[i][199:-200]  #remove the 1st and last 200bp because lower coverage at ends of assembly
			c_trim=scaffolds_cov[i][199:-200]
			count_0_cov.append(scaffolds_cov[i].count(0))
			loci_trim.append(l_trim)
			cov_trim.append(c_trim) 
#			loci_sub=loci_trim[0::1000]
#			cov_sub=cov_trim[0::1000]
			cov_avg=mean(c_trim)
			cov_stdev=numpy.std(c_trim)
			scaffold_avg_cov.append(cov_avg)
			scaffold_stdev.append(cov_stdev)
			print "Average coverage: ", round(cov_avg,2)
			print "Standard Deviation: ", round(cov_stdev,2)
	print count_0_cov
	return loci_trim, cov_trim, scaffold_avg_cov, scaffold_stdev,length, num_scaffolds, count_0_cov

def Window_Step(data,window,step):
	start=0
	out=[]
	while start<=len(data)-window:
		out.append(data[start:start+window])
		start+=step
#	print len(data[start:])
	out.append(data[start:])  #should pick up last fragment
	return out


def Find_Discontinuities2(cov_step,loci_step, cov_avg,cov_stdev):
	avg_cov_step=[]
	mean_based=[]
	for i in cov_step:
		avg_cov_step.append(mean(i))
		if mean(i)>(1.5*cov_avg):   #this is still the window size set below
			ind=cov_step.index(i)
			# loci_range=[min(loci_step[ind]),max(loci_step[ind])]
			mean_based.append([min(loci_step[ind]),max(loci_step[ind])]) #,mean(i)])  #removed mean for combine3 input
	return mean_based
"""
def combine(mean_based,step_size, length_req):
	to_combine=[]
	mean_based_out=["Min_Loci\tMax_Loci\tAverage Coverage"]
	count=0
	for a in range(len(mean_based)):  #can sum up range
		if a==len(mean_based)-1 and mean_based[a][0]-50==mean_based[a-1][0] : #have come to end of list and part of increased coverage
			to_combine.append(mean_based[a]) #have come to end of list and not part of increased coverage
			count+=1
			newmin=to_combine[0][0]
			newmax=to_combine[-1][1]
			if newmax-newmin+1>=length_req:
				avgs=[]
				for n in to_combine:
					avgs.append(str(n[2]))
				mean_based_out.append("\t".join([str(newmin),str(newmax),";".join(avgs)]))
				to_combine=[]
				count+=1
			else:
				print "Did not meet length req", to_combine	
				to_combine=[]  #didn't meet length requirement	
		elif a==len(mean_based)-1:
			mean_based_out.append("\t".join(map(str,mean_based[a])))
			count+=1
			break
		elif mean_based[a][0]+step_size==mean_based[a+1][0]:  #increased ranges are adjacent
			to_combine.append(mean_based[a])
		elif len(to_combine)>0:  #increased ranges are not adjacent but previous were, so combine
			newmin=to_combine[0][0]
			newmax=to_combine[-1][1]
			if newmax-newmin+1>=length_req:
				avgs=[]
				for n in to_combine:
					avgs.append(str(n[2]))
				mean_based_out.append("\t".join([str(newmin),str(newmax),";".join(avgs)]))
				to_combine=[]
			else:
				print "Did not meet length req", to_combine	
				to_combine=[]  #didn't meet length requirement		
		else:
			 continue #not nearby any adjacent increase in coverage, thus won't meet length req. (since step size is currently 50)
	return mean_based_out


def combine2(mean_based,step_size, length_req):  #the purpose of this now is to combine ranges so that we can then do Webster's to find Discontinuities for final approximate size of increased coverage
	to_combine=[]
	mean_based_out=[]
	count=0
	for a in range(len(mean_based)):  #can sum up range
		if a==len(mean_based)-1 and mean_based[a][0]-50==mean_based[a-1][0] : #have come to end of list and part of increased coverage
			to_combine.append(mean_based[a]) #have come to end of list and not part of increased coverage
			count+=1
			newmin=to_combine[0][0]
			newmax=to_combine[-1][1]
			if newmax-newmin+1>=length_req:
				avgs=[]
				for n in to_combine:
					avgs.append(str(n[2]))
				mean_based_out.append([newmin,newmax])
				to_combine=[]
				count+=1
			else:
				print "Did not meet length req", to_combine	
				to_combine=[]  #didn't meet length requirement	
		elif a==len(mean_based)-1:
			mean_based_out.append(mean_based[a][:2])
			count+=1
			break
		elif mean_based[a][0]+step_size==mean_based[a+1][0]:  #increased ranges are adjacent
			to_combine.append(mean_based[a])
		elif len(to_combine)>0:  #increased ranges are not adjacent but previous were, so combine
			newmin=to_combine[0][0]
			newmax=to_combine[-1][1]
			if newmax-newmin+1>=length_req:
				avgs=[]
				for n in to_combine:
					avgs.append(str(n[2]))
				mean_based_out.append([newmin,newmax])
				to_combine=[]
			else:
				print "Did not meet length req", to_combine	
				to_combine=[]  #didn't meet length requirement		
		else:
			 continue #not nearby any adjacent increase in coverage, thus won't meet length req. (since step size is currently 50)
	return mean_based_out  #output is ranges of regions of interest

def combine2(mean_based,step_size):  #the purpose of this now is to combine ranges so that we can then do Webster's to find Discontinuities for final approximate size of increased coverage
	mean_based_out=[]
	starts=[]
	stops=[]
	for i in mean_based:
		starts.append(i[0])
		stops.append(i[1])
	groups=numpy.split(starts,numpy.where(numpy.diff(starts) !=step_size)[0]+1)
#	print len(groups)
	for g in groups:
#		print g
		ls=g.tolist()
		minimum=ls[0]
		maximum=stops[starts.index(ls[-1])]
#		print minimum,maximum
		mean_based_out.append([minimum,maximum])
	return mean_based_out	
"""

def combine3(mean_based,step_size):  #combine any overlapping windows; mean_based is [start,end,avg_cov]
	merged=[mean_based[0]]
	for current in mean_based:
		previous=merged[-1]
		if current[0]<=previous[1]:
			previous[1]=max(previous[1],current[1])
		else:  #start new range
			merged.append(current)
	mean_based_out=merged
	print "Mean_based_out", mean_based_out
	return mean_based_out
			
			

def Websters_Discontinuities(mean_based_out,loci_trim,cov_trim,cutoff):  #NOT USING WEBSTER'S ANYMORE
	out=["Start_Discontinuity\tEnd_Discontinutity\tLength\tAvg_Cov"]
	for t in mean_based_out:
		index_start=loci_trim.index(t[0])
		index_end=loci_trim.index(t[1])
		region_cov=cov_trim[index_start:(index_end+1)]  #add 1 so includes last value
		region_loc=loci_trim[index_start:(index_end+1)]
		#now find discontinutities within region_cov
#		print t, len(region_cov),region_cov
		window=50
		step=50
		st=0
		cov_step=[]
		loci_step=[]
		avg_step=[]
		while st<=len(region_cov)-window:  
			cov_step.append(region_cov[st:st+window])
			loci_step.append(region_loc[st:st+window])
			avg_step.append(mean(region_cov[st:st+window]))		
			st+=step
		if len(region_cov[st:st+window])>0.5*window:  #will only be >0 if change window and step size
			cov_step.append(region_cov[st:])  #should pick up last fragment	
			loci_step.append(region_loc[st:])	
			avg_step.append(mean(region_cov[st:st+window]))	
#		print "avg_step", avg_step
		for a in avg_step:
			if a>cutoff:  #start of region of interest
				s_ind=min(loc for loc, val in enumerate(avg_step) if val==a)
				break
		for b in reversed(avg_step):
			if b>cutoff:
				e_ind=max(loc for loc, val in enumerate(avg_step) if val==b)
				break
		min_loci=min(loci_step[s_ind])
		max_loci=max(loci_step[e_ind])
		len_discrep=max_loci-min_loci
#		print "Range: %s - %s; Length: %s" % (min_loci, max_loci,len_discrep)
		if len_discrep<1500:  #CAN PLAY WITH THIS NUMBER.  THIS WOULD BE A 500BP GENE
			print "Too short: (%s - %s)" % (min_loci,max_loci)
			continue  #don't want it
		else:
			cov_discrep=mean([cov_trim[loci_trim.index(min_loci):loci_trim.index(max_loci)+1]])
#			print "Coverage Discrepancy:", cov_discrep
			out.append("\t".join(map(str,[min_loci,max_loci,len_discrep,cov_discrep])))
	return out


"""
#This used websters to find discrepancies for range.  Was too stringent	
		#go from 5' until hit max
		#then go from 3' until hit max
		#??? check for any discrepancies in middle (purely for flagging)????
		print len(cov_step)
		webster=[]
		for i in cov_step:
			A=i[:len(i)/2]
			B=i[len(i)/2:]
			discrepancy=mean(A)-mean(B) #this splits the segment into two windows (so 1/2 of window_size) and compares the two
			webster.append(discrepancy)
		print len(webster)
		print webster
		print max(webster), min(webster)
		#find most neg at beginning (start) then most positive at end (end)
		s=webster[0]
		s_ind=0
		for w in webster[1:]:
			if w<s:
				s=w
				s_ind=webster.index(w)
			else:
				break
		e=webster[-1]
		for w in reversed(webster[:-1]):
			if w>e:
				e=w
				e_ind=max(loc for loc, val in enumerate(webster) if val==w)  
			else:
				break
		min_loci=int(median(range(loci_step[s_ind][0],loci_step[e_ind][1]+1)))
		max_loci=int(median(range(loci_step[e_ind][0],loci_step[e_ind][1]+1)))
		discrep_loci=[min_loci,max_loci]  #records location of A/B junction
		print "Discrep_loci:", discrep_loci	
		len_discrep=max_loci-min_loci
		if len_discrep<1500:  #CAN PLAY WITH THIS NUMBER.  THIS WOULD BE A 500BP GENE
			print "Too short: (%s - %s)" % (min_loci,max_loci)
			continue  #don't want it
		else:
			print len([cov_trim[loci_trim.index(min_loci):loci_trim.index(max_loci)+1]])
			cov_discrep=mean([cov_trim[loci_trim.index(min_loci):loci_trim.index(max_loci)+1]])
			print cov_discrep
			out.append("\t".join(map(str,[min_loci,max_loci,len_discrep,cov_discrep])))
	return out
"""		

window_size=3000
step_size=100
#length_req=3000  #the length of a discontinuity to be considered  (this set to window size means no length req?)
o="Mean-based_stats:W%s_S%s" % (window_size,step_size)
print o
sra_csv=open(sys.argv[1],"r")  #csv of SRA IDs of specific alignments to look at
cov_dir=sys.argv[2]  #directory where Depth.txt files are housed
out_dir=jn(cov_dir,o)
if os.path.isdir(out_dir):
	print "The output directory for Stats already exists"
else:
	subprocess.call(["mkdir",out_dir])

sras=sra_csv.readline().strip("\r\n").split(",")
#sras=["SRR3493382_GCF_001877035.1_blasr_sorted_unique_depth.txt"]

onlyfiles=[f for f in listdir(cov_dir) if isfile(jn(cov_dir,f))]  #1st one always seems to be ".DS_Store" so will ignore 1st file
if onlyfiles[0]==".DS_Store" or onlyfiles[0]=="._.DS_Store":
	files=onlyfiles[1:]
else:
	files=onlyfiles

print len(sras)


summary_out=jn(out_dir, sys.argv[3])
if os.path.exists(summary_out)==False:
	handle=open(summary_out,"w")
	titleline=["SRA\tNum_Scaffolds\tScaffold\tLength\tAverage_Coverage\tStDev\tCount_0_cov\tCount_0_cov_trimmed\tProportionTrim_zero_coverage\tSkewness (~0)\tKurtosis {~0)\tMean-based_Regions_of_Interest\tNotes\n"]
	handle.write("\n".join(titleline))
	handle.close()

count=0


print "Checking for files..."
for s in sras:
	found=False	
	for f in files:
		if s in f:
			found=True
			break
	assert found==True, "A file for %s was not found. Cancelled processing" % s
for s in sras:
	summary=[]
	count+=1
	print "NUMBER %s" % count
	for f in files:
		if s in f:
			dh=open(jn(cov_dir,f),"r")
			break
	print s
	fileout=f.strip("_depth.txt")
	loci=[]
	cov=[]
	writeout=[]
	while True:
		x=dh.readline().strip("\n").split("\t")
		if x[0]=="":
			break
		else:	
			loci.append(int(x[1]))
			cov.append(int(x[2]))
#	print len(loci), len(cov)
	print "Read in data for %s" % f
	loci_trim, cov_trim, scaffold_avg_cov, scaffold_stdev,length,num_scaffolds,count_0_cov=fmt_data(loci,cov)
#	print len(loci_trim), len(cov_trim), scaffold_avg_cov, scaffold_stdev
	for i in range(len(loci_trim)):  #this will make a line for each scaffold in the summary file
		notes=[]
		print "Count 0: ", cov_trim[i].count(0)
		print "Length: ", len(cov_trim[i])
		perc_zero_cov=str(round(cov_trim[i].count(0)/len(cov_trim[i]),3))
		skewness=str(round(skew(cov_trim[i]),2))
		kurt=str(round(kurtosis(cov_trim[i]),2))
#		skewness=""  #FIX THIS AND KURT IF GET SCIPY.STATS TO WORK AGAIN
 #		kurt=""
		print perc_zero_cov, skewness, kurt
		writeout.append("Scaffold: %s" % i)
		if len(loci_trim[i])<5000:
			notes.append("Too short to analyze (<5000)")  

#"SRA\tNum_Scaffolds\tScaffold\tLength\tAverage_Coverage\tStDev\tCount_0_cov\tCount_0_cov_trimmed\tProportionTrim_zero_coverage\tSkewness (~0)\tKurtosis {~0)\tMean-based_Regions_of_Interest\tNotes"
			summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i],str(count_0_cov[i]),str(cov_trim[i].count(0)),perc_zero_cov,skewness,kurt,"",";".join(notes)])))
			writeout.append("Too short to analyze (<5000)\n")
		elif cov_trim[i].count(0)==len(cov_trim[i]):  # none of scaffold has coverage
			notes.append("No coverage of scaffold")  
			summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i],str(count_0_cov[i]),str(cov_trim[i].count(0)),perc_zero_cov,skewness,kurt,"",";".join(notes)])))
			writeout.append("No coverage of scaffold\n")
		else:
			loci_step=Window_Step(loci_trim[i],window_size,step_size)  
			cov_step=Window_Step(cov_trim[i],window_size,step_size)
			print "Number of Windows: ", len(cov_step)
			mean_based=Find_Discontinuities2(cov_step,loci_step,scaffold_avg_cov[i],scaffold_stdev[i])					
			if scaffold_avg_cov[i]<50:
				notes.append("Low coverage")
			if scaffold_stdev[i]>0.5*scaffold_avg_cov[i]:
				notes.append("High standard deviation")
			if len(mean_based)>0:  #found differences
				mean_based_out=combine3(mean_based,step_size)
				cutoff=scaffold_avg_cov[i]*1.5
				out=Websters_Discontinuities(mean_based_out,loci_trim[i],cov_trim[i],cutoff)
				summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i],str(count_0_cov[i]),str(cov_trim[i].count(0)),perc_zero_cov,skewness,kurt,len(out)-1,";".join(notes)])))
				print "Mean-based regions of interest: ", len(out)-1  #1st is title line
				writeout.append("\n".join(out)+"\n")
			else:  
				print "No regions of interest"
				summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i], str(count_0_cov[i]),str(cov_trim.count(0)),perc_zero_cov,skewness,kurt,"0",""])))
	handle=open(jn(out_dir,fileout+"_mean_based.txt"),"w")
	handle.write("\n".join(writeout))
	handle.close()	
	handle2=open(summary_out,"a")
	handle2.write("\n".join(summary)+"\n")
	handle2.close()



