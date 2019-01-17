import subprocess
#from numpy import mean


#make a class (learn how to do this); pipeline these self.something as input_mapping
#an example is jgi_mg_assembly

def build_Bowtie2_index(genome_fasta):
    idx="genome_index"
    subprocess.call(["bowtie2-build", genome_fasta, idx])

def run_Bowtie2(reads_input,idx):
    alignout_sam="alignout_sam.sam"
	proc=subprocess.Popen(["bowtie2", "-t", "-p", "3", "-x", idx, reads_input, "-S", alignout_sam],stderr=subprocess.PIPE)
	output_align_stats=proc.communicate()[1]
    print output_align_stats

def samtools_convert_depth(alignout_sam,bamout,sortout,depthout):
	#samtools conversion sam to bam
	f=open(bamout,"w")
	subprocess.call(["samtools", "view", "-S", "-b", alignout_sam], stdout=f)
	f.close()
	#samtools bam sort
	subprocess.call(["samtools", "sort",bamout, sortout])
    #samtools create depth file
	f=open(depthout,"w")
	subprocess.call(["samtools","depth", sortout], stdout=f)  #Check! sorting may add another .bam to end of file name
	f.close()
