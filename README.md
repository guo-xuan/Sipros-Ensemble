# Sipros Ensemble

Sipros is a database-searching algorithm for peptide and protein identification in shotgun meta/proteomics. The detailed user manual of the database-searching and how to use it to achieve best results is provided here: [http://siprosensemble.omicsbio.org/user-manual](http://siprosensemble.omicsbio.org/user-manual). This is a quick start guide generally for developers and testers. Users with limited experience with MS-based database-searching are advised to use the user manual.

### Current Version
* v1.0

### Setup and Installation

#### Basic Dependencies

1. GNU GCC or Intel C++  with C++11 support i.e. gcc4.9+ or above, icpc15.0+ or above.
2. MPI Library with MPI-3 support i.e. OpenMPI 1.8 and above or cray-mpich/7.4.0 and above. By default the mpic++ wrapper is needed. If you are on a Cray cluster and the wrapper is "CC". You will need to edit the compiler.mk file. Uncomment the line "MCC := CC" and comment out "MCC := mpic++".   
3. Python 2.7.2 or above Python 2 versions, numpy 1.11.2 or above, scipy 0.13.3 or above, [scikit-learn](http://scikit-learn.org/) 0.17.1 or above.
 
#### Installation Steps
1. Download the tarball with compiled executables for Linux with GCC 4.9 and above from  [https://github.com/guo-xuan/Sipros-Ensemble/releases](https://github.com/guo-xuan/Sipros-Ensemble/releases). The code has been tested only on Linux.
2. If you decide to download the source code, use the following commands to build:
  1. OpenMP version "make openmp".
  2. MPI version version "make mpi" 
  3. All the versions can be built with "make all"
If compiled successfully, the required executables will be in `bin` directory and the various `runSipros...` scripts can be used to run the database-searching. 

#### Configure File Setting

\# is for comments.

[] is used for section name, e.g., [Section Name].

= is used for assigning features, e.g., Search_Type = Regular

{} is used for specifying key value, e.g., PTM{!} = NQR

Currently, there are 35 symbols available for specifying ptms, which are

~ ! @ $ % ^ & * ( ) _ + ` - | \ : " ; ' < > ? . / 1 2 3 4 5 6 7 8 9 0

Please don't use these reserved symbols: { } # [ ] = ,

Neutral loss can be specified by PTM{1to2}, e.g., PTM{>to|} = ST. If symbol2 is nothing, it can be specified by PTM{1to}, e.g., PTM{>to} = ST.

#### Generate Reverse Sequences
```
python reverseseq.py -i original_database_file -o output_database_file
```
The step will generate a new database file with reverse sequences. Update the path of `FASTA_Database` in the configuration file.

#### Quickly Running The Database-searching and Filtering/Assembling

There are two basic versions of the database-searching: one for running on a single machine and another for running with MPI on a cluster.  

* __Single Machine Version:__ This version of the assembler should be used if you are going to run the database-searching on a single machine with one or more cores. The searching is invoked through a run script `./runSipros.sh`. The quick start command as shown below will be used in a batch job submission script or directly typed on the command line terminal.   

```
#!/bin/bash

# Single MS2 file
runSipros.sh -o ${output_dir} -f data.ms2 -c SiprosConfig.cfg

# Multiple MS2 files in a folder
runSipros.sh -o ${output_dir} -w workingdirectory -c SiprosConfig.cfg

```
Use `./runSipros.sh -h` for help information.

* __MPI Version:__ This version of the assembler should be used if you are going to run the assembler with MPI support on a cluster. The run script to invoke the assembler depends on the cluster management and job scheduling system.

	1. If you have ORTE i.e. __mpirun__ is available, invoke the assembler using the run script `runDisco-MPI.sh`. 
	2. If you have SLRUM i.e. __srun__ is available, invoke the assembler using the run script `runDisco-MPI-SLRUM.sh`.
	3. If you have ALPS i.e. __aprun__ is available, invoke the assembler using the run script `runDisco-MPI-ALPS.sh`.
 
For the basic MPI version make sure the RAM on the nodes is more than the disk space size of the reads. If you have a large dataset, then use the Remote Memory Access (RMA) version. The RMA version of the assembler will equally distribute about 70% of the memory usage across all the MPI nodes. The quick start commands are:
```
#!/bin/bash

### MPI Verion 
### Seperated paired end reads
runDisco-MPI.sh -d ${output_dir} -in1 {read_1.fastq} -in2 ${read2_2.fastq} -o ${OP_PREFIX} 

### MPI Remote Memory Access(RMA) Verion 
### Seperated paired end reads
runDisco-MPI.sh -d ${output_directory} -in1 {read_1.fastq}  -in2 ${read2_2.fastq} -o ${OP_PREFIX} -rma 

```
Use `runDisco-MPI.sh -h` for help information.

### Guide to Assembly of Raw Metagenomic Illumina data

The raw Illumina sequences need to be preprocessed before assembly with Disco. Disco provides wrapper scripts to perform preprocessing with BBTools. Please see user manual for more details: [http://disco.omicsbio.org/user-manual](http://disco.omicsbio.org/user-manual)

#### Preprocessing of the Illumina data

Since Disco works best with reads without errors, preprocessing plays an important role in deciding the quality of the assembly results. The 3 basic pre-processing steps are trimming, filtering and eror correction.

##### Trimming, filtering, (merging), and eror correction

We have tested Brian Bushnell's suite of tools [BBTools](http://sourceforge.net/projects/bbmap/files/) extensively on Illumina data and have obtained good results. Suppose the Illumina reads data set is called `$reads`, the steps we recommend are following:

```
#!sh

# Use bbduk.sh to quality and length trim the Illumina reads and remove adapter sequences
# 1. ftm = 5, right trim read length to a multiple of 5
# 2. k = 23, Kmer length used for finding contaminants
# 3. ktrim=r, Trim reads to remove bases matching reference kmers to the right
# 4. mink=11, look for shorter kmers at read tips down to 11 bps
# 5. qhdist=1, hamming distance for query kmers
# 6. tbo, trim adapters based on where paired reads overlap
# 7. tpe, when kmer right-trimming, trim both reads to the minimum length of either
# 8. qtrim=r, trim read right ends to remove bases with low quality
# 9. trimq=10, regions with average quality below 10 will be trimmed.
# 10. minlength=70, reads shorter than 70bps after trimming will be discarded.
# 11. ref=$adapters, adapters shipped with bbnorm tools
# 12. –Xmx8g, use 8G memory
# 13. 1>trim.o 2>&1, redirect stderr to stdout, and save both to file *trim.o*
adapters= bbmap_dir/resources/adapters.fa
phiX_adapters= bbmap_dir/resources/phix174_ill.ref.fa.gz
bbduk.sh in=$reads out=trim.fq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=${adapters} ftm=5 qtrim=r trimq=10
bbduk.sh in=trim.fq.gz out=filter.fq.gz ref=$phiX_adapters hdist=1 k=31
```

##### Error correction with Tadpole

Tarpole is a memory efficient error correction tool from the bbtools package that runs within reasonable time. We suggest using this for error correction. 
```
#!bash
# 1. mode=correct, use tadpole for correction
# 2. ecc=t, error correct via kmer counts
tadpole.sh in=filter.fq.gz out=ecc.fq.gz mode=correct prefilter=1 prealloc k=31
#If the above goes out of memory, try
tadpole.sh in=filter.fq.gz out=ecc.fq.gz mode=correct prefilter=2 prealloc k=31
```

### Assembly of Error Corrected Data

#### Assembly on a Single Node

The Disco assembler is invoked through the run script `./runDisco.sh`. The basic quick start commands with default parameters are as follows. The default parameters are based on empherical tests on real metagenomic datasets.     

```
#!/bin/bash

# Separated paired end reads
runDisco.sh -d ${output_directory} -in1 {read_1.fastq}  -in2 ${read2_2.fastq} -n ${num_threads} -m {max_mem_usage} -o ${64gen} 

# Interleaved paired end reads
runDisco.sh -d ${output_directory} -inP {read_P.fastq} -n ${num_threads} -m {max_mem_usage} -o ${64gen}

# Single end reads
runDisco.sh -d ${output_directory} -inS {read.fastq} -n ${num_threads} -m {max_mem_usage} -o ${64gen} 
```
For all the options of Disco, use `./runDisco.sh -h`

In case the program crashes due to exceeding wall clock time, the assembler can be restarted with the same command. 

#### Assembly on a Distributed Nodes

The assembler can be run on a distributed machine using the three distributed assembly scripts. 

#### Assembly Run Script Options

Usage:

   runDisco.sh [OPTION]...<PARAM>...


<PARAMS>

   -inS	 single read filenames (comma seperated fasta/fastq/fastq.gz file).

   -in1	 forward paired read filename (single fasta/fastq/fastq.gz file).

   -in2	 reverse paired read filename (single fasta/fastq/fastq.gz file).

   -inP	 interleaved paired read filenames (comma seperated fasta/fastq/fastq.gz file).

   -d	 output directory path (DEFAULT: current directory).

   -o	 output filename prefix (DEFAULT: disco).

<OPTIONS>

   -h	 help.

   -m	 maximum memory to be used (DEFAULT: 125 GB).

   -n	 number of threads (DEFAULT: 32).

   -obg	 only build overlap graph (DEFAULT: False).

   -osg	 only simplify existing overlap graph (DEFAULT: False).
   
   -p	 assembly parameter file for 1st assembly iteration.

   -p2	 assembly parameter file for 2nd assembly iteration.

   -p3	 assembly parameter file for 3rd assembly iteration.


The assembly script has basic options to specify required parameters. 

#### Controlling memory usage

The memory usage of Disco can be controlled using the `-m` option to the run script as shown above. The default memory usage is to take all the system resources. In case that has to be avoided or the program crashes ot is too slow due to memory page swapping, the user can set a ubber bound on the memory. The minumum memory to assemble a dataset is:

```
Min Required Memory (GB) = (Disk Space of Reads) + (1GB * num_threads)
``` 
The program will run faster if more memory is made available.

#### Restarting Disco for repeat assembly and handling assembly crashes

Disco assembler can be restarted with changed assembly and scaffolding parameters using the `-osg` option. Setting this option while invoking `runDisco.sh` will reuse the overlap graph constructed earlier and only perform the graph simplification step. This will significantly reduce executime time of assemblies on the same dataset with different parameters.    

Disco assembler can also be restarted after a crash caused due to exceeding wall clock time or out of memory errors. The job must be restarted with the same command as before and Disco will attempt to continue the assembly. Do not set the `-osg` option in this case.   

#### Setting assembly parameters

The assembly parameters can be modified to attempt better assembly. This can be done through a parameter file passed using the `-p` parameter to the run script. 

The configurable parameters are described in the user manual [http://disco.omicsbio.org/user-manual](http://disco.omicsbio.org/user-manual).

The default configuration parameters are in disco.cfg, disco_2.cfg, and disco_3.cfg. 

#### Disco Assembler Output

Please see the OUTPUT.md file for description of the output files.  

### Questions?

* [Xuan Guo](mailto:xuan_guo@outlook.com)
* [Chongle Pan](mailto:chongle.pan@gmail.com)
