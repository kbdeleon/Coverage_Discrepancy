# Map sra to assembly.
# Short reads done with Bowtie2, long reads done with Blasr (from h5 files)
# Output is bam
# samtools depth -aa <mapped.bam> <output_depth.txt>

# Input here is currently depth.txt file
# This program was written for batch input and so
# output is a summary file containing stats of all depth files analyzed as well as a small file of the disrepancy information for each sample

#consider output graph


#NEED TO CONVERT TO PYTHON3!!!


import sys
import subprocess
import os
from os.path import isfile
from os.path import join as jn
from os import listdir
from numpy import mean,median
import numpy
from scipy.stats import kurtosis,skew  #2018-10-24 scipy.stats throwing error with conda


def fmt_data(loci,cov):
	scaffolds_loci = []
	scaffolds_cov = []
	temp_loci = []
	temp_cov = []
	loci_trim = []
	cov_trim = []
	scaffold_avg_cov = []
	scaffold_stdev = []
	length = []
	num_scaffolds = []
	count_0_cov = []
	ct = 0
	for i in range(len(loci)):
		if loci[i] > ct:
			temp_loci.append(loci[i])
			temp_cov.append(cov[i])
			ct = loci[i]
		else:  #new scaffold
			scaffolds_loci.append(temp_loci)
			scaffolds_cov.append(temp_cov)
			temp_loci = [loci[i]]
			temp_cov = [cov[i]]
			ct = loci[i]
	scaffolds_loci.append(temp_loci)
	scaffolds_cov.append(temp_cov)
	num_scaffolds = len(scaffolds_loci)
	for i in range(len(scaffolds_loci)):
		length.append(len(scaffolds_loci[i]))
		if len(scaffolds_loci[i])<401:
			print "Scaffold <400 bp"  #DO SOMETHING HERE TO ALERT. would result in a scaffold of len 0 after trimming below
		else:
			l_trim = scaffolds_loci[i][199:-200]  #remove the 1st and last 200bp because lower coverage at ends of assembly
			c_trim = scaffolds_cov[i][199:-200]
			count_0_cov.append(scaffolds_cov[i].count(0))
			loci_trim.append(l_trim)
			cov_trim.append(c_trim)
			cov_avg = mean(c_trim)
			cov_stdev = numpy.std(c_trim)
			scaffold_avg_cov.append(cov_avg)
			scaffold_stdev.append(cov_stdev)
			print "Average coverage: ", round(cov_avg,2)
			print "Standard Deviation: ", round(cov_stdev,2)
	print count_0_cov
	return loci_trim, cov_trim, scaffold_avg_cov, scaffold_stdev,length, num_scaffolds, count_0_cov

def Window_Step(data,window,step):
	start = 0
	out=[]
	while start <= len(data)-window:
		out.append(data[start:start+window])
		start += step
	out.append(data[start:])  # pick up last fragment
	return out


def Find_Discontinuities(cov_step,loci_step, cov_avg,cov_stdev,times_above_avg):
	avg_cov_step = []
	mean_based = []
	for i in cov_step:
		avg_cov_step.append(mean(i))
		if mean(i) > (times_above_avg * cov_avg):   #this is still the window size set below
			ind = cov_step.index(i)
			mean_based.append([min(loci_step[ind]),max(loci_step[ind])])
	return mean_based

def combine(mean_based,step_size):  #combine any overlapping windows; mean_based is [start,end,avg_cov]
	merged = [mean_based[0]]
	for current in mean_based:
		previous = merged[-1]
		if current[0] <= previous[1]:
			previous[1] = max(previous[1],current[1])
		else:  #start new range
			merged.append(current)
	mean_based_out = merged
	print "Mean_based_out", mean_based_out
	return mean_based_out


def Websters_Discontinuities(mean_based_out,loci_trim,cov_trim,cutoff,webster_window,webster_step):
	out=["Start_Discontinuity\tEnd_Discontinutity\tLength\tAvg_Cov"]
	for t in mean_based_out:
		index_start = loci_trim.index(t[0])
		index_end = loci_trim.index(t[1])
		region_cov = cov_trim[index_start:(index_end+1)]  #add 1 so includes last value
		region_loc = loci_trim[index_start:(index_end+1)]
		#now find discontinutities within region_cov
		window = webster_window
		step = webster_step
		st = 0
		cov_step = []
		loci_step = []
		avg_step = []
		while st <= len(region_cov)-window:
			cov_step.append(region_cov[st:st+window])
			loci_step.append(region_loc[st:st+window])
			avg_step.append(mean(region_cov[st:st+window]))
			st+=step
		if len(region_cov[st:st+window]) > 0.5*window:  #will only be >0 if change window and step size; this is not currently being used, so need to test if change
			cov_step.append(region_cov[st:])  #should pick up last fragment
			loci_step.append(region_loc[st:])
			avg_step.append(mean(region_cov[st:st+window]))
		for a in avg_step:
			if a > cutoff:  #start of region of interest
				s_ind = min(loc for loc, val in enumerate(avg_step) if val==a)
				break
		for b in reversed(avg_step):
			if b > cutoff:
				e_ind = max(loc for loc, val in enumerate(avg_step) if val==b)
				break
		min_loci = min(loci_step[s_ind])
		max_loci = max(loci_step[e_ind])
		len_discrep = max_loci-min_loci
		if len_discrep < 1500:  #CAN PLAY WITH THIS NUMBER.  THIS WOULD BE A 500BP GENE
			print "Too short: (%s - %s)" % (min_loci,max_loci)
			continue  #don't want it
		else:
			cov_discrep = mean([cov_trim[loci_trim.index(min_loci):loci_trim.index(max_loci)+1]])
#			print "Coverage Discrepancy:", cov_discrep
			out.append("\t".join(map(str,[min_loci,max_loci,len_discrep,cov_discrep])))
	return out

window_size = 3000
step_size = 100


loci_trim, cov_trim, scaffold_avg_cov, scaffold_stdev,length,num_scaffolds,count_0_cov = fmt_data(loci,cov)
for i in range(len(loci_trim)):  #this will make a line for each scaffold in the summary file
	notes = []
	print "Count 0: ", cov_trim[i].count(0)
	print "Length: ", len(cov_trim[i])
	perc_zero_cov = str(round(cov_trim[i].count(0)/len(cov_trim[i]),3))
	skewness = str(round(skew(cov_trim[i]),2))
	kurt = str(round(kurtosis(cov_trim[i]),2))
	print perc_zero_cov, skewness, kurt
	writeout.append("Scaffold: %s" % i)
	if len(loci_trim[i])<5000:
		notes.append("Too short to analyze (<5000)")

#"SRA\tNum_Scaffolds\tScaffold\tLength\tAverage_Coverage\tStDev\tCount_0_cov\tCount_0_cov_trimmed\tProportionTrim_zero_coverage\tSkewness (~0)\tKurtosis {~0)\tMean-based_Regions_of_Interest\tNotes"
		summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i],str(count_0_cov[i]),str(cov_trim[i].count(0)),perc_zero_cov,skewness,kurt,"",";".join(notes)])))
		writeout.append("Too short to analyze (<5000)\n")
	elif cov_trim[i].count(0) == len(cov_trim[i]):  # none of scaffold has coverage
		notes.append("No coverage of scaffold")
		summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i],str(count_0_cov[i]),str(cov_trim[i].count(0)),perc_zero_cov,skewness,kurt,"",";".join(notes)])))
		writeout.append("No coverage of scaffold\n")
	else:
		loci_step = Window_Step(loci_trim[i],window_size,step_size)
		cov_step = Window_Step(cov_trim[i],window_size,step_size)
		print "Number of Windows: ", len(cov_step)
        times_above_avg = 1.5
		mean_based = Find_Discontinuities(cov_step,loci_step,scaffold_avg_cov[i],scaffold_stdev[i], times_above_avg)
		if scaffold_avg_cov[i] < 50:
			notes.append("Low coverage")
		if scaffold_stdev[i] > 0.5*scaffold_avg_cov[i]:
			notes.append("High standard deviation")
		if len(mean_based) > 0:  #found differences
			mean_based_out = combine(mean_based,step_size)
			cutoff = scaffold_avg_cov[i]*1.5
            webster_window = 50 #consider making an option in KBase App
            webster_step = 50  #consider making an option in KBase App
			out = Websters_Discontinuities(mean_based_out,loci_trim[i],cov_trim[i],cutoff, webster_window, webster_step)
			summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i],str(count_0_cov[i]),str(cov_trim[i].count(0)),perc_zero_cov,skewness,kurt,len(out)-1,";".join(notes)])))
			print "Mean-based regions of interest: ", len(out)-1  #1st is title line
			writeout.append("\n".join(out)+"\n")
		else:
			print "No regions of interest"
			summary.append("\t".join(map(str,[s,num_scaffolds,"scaffold_"+str(i+1),length[i],scaffold_avg_cov[i],scaffold_stdev[i], str(count_0_cov[i]),str(cov_trim.count(0)),perc_zero_cov,skewness,kurt,"0",""])))
handle = open(jn(out_dir,fileout+"_mean_based.txt"),"w")
handle.write("\n".join(writeout))
handle.close()
handle2 = open(summary_out,"a")
handle2.write("\n".join(summary)+"\n")
handle2.close()
