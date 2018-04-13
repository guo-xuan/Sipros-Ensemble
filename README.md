# Sipros Ensemble

Sipros is a database-searching algorithm for peptide and protein identification in shotgun meta/proteomics. To run Sipors, you need one or more input spectral files in mzML, FT2, or ms2 formats and a SiprosConfig.cfg file. Then issue a command such as:
```
Sipros_Openmp -f input.mzML -c SiprosConfig.cfg -o destination_folder
Sipros_Openmp -f input.FT2 -c SiprosConfig.cfg -o destination_folder
Sipros_Openmp -f input.ms2 -c SiprosConfig.cfg -o destination_folder
Sipros_MPI -f input.mzML -c SiprosConfig.cfg -o destination_folder
Sipros_Openmp -w input_folder -c SiprosConfig.cfg -o destination_folder
Sipros_MPI -w input_folder -c SiprosConfig.cfg -o destination_folder
```

The detailed user manual of the database-searching and how to use it to achieve best results is provided here: [http://siprosensemble.omicsbio.org/user-manual](http://siprosensemble.omicsbio.org/user-manual). This is a quick start guide generally for developers and testers. Users with limited experience with MS-based database-searching are advised to use the user manual.

### Current Version
* v1.1

### Setup and Installation

#### Basic Dependencies

1. GNU GCC or Intel C++  with C++11 support i.e. gcc4.9+ or above, icpc15.0+ or above.
2. MPI Library with MPI-3 support i.e. OpenMPI 1.8 and above or cray-mpich/7.4.0 and above. By default the mpic++ wrapper is needed. If you are on a Cray cluster and the wrapper is "CC". You will need to edit the compiler.mk file. Uncomment the line "MCC := CC" and comment out "MCC := mpic++".   
3. Python 2.7.2 or above Python 2 versions, numpy 1.11.2 or above, scipy 0.13.3 or above, [scikit-learn](http://scikit-learn.org/) 0.17.1 or above. A guide to setup Python environment is [here](#PythonSetup).
 
#### Installation Steps
1. Download the tarball with compiled executables for Linux with GCC 4.9 and above from  [https://github.com/guo-xuan/Sipros-Ensemble/releases](https://github.com/guo-xuan/Sipros-Ensemble/releases). The code has been tested only on Linux.
2. If you decide to download the source code, use the following commands to build:
  1. OpenMP version "make openmp".
  2. MPI version version "make mpi" 
  3. All the versions can be built with "make all"
If compiled successfully, the required executables will be in `bin` directory and the various `runSipros...` scripts can be used to run the database-searching. 

#### <a name="config"></a>Configure File Setting

A sample configuration file, `SiprosConfig.cfg`, is available in `configs` directory. The configuration settings used for benchmarking is in `SiprosConfigBenchmark.cfg`.

`#` is for comments.

`[]` is used for section name, e.g., `[Section Name]`.

`=` is used for assigning features, e.g., `Search_Type = Regular`

`{}` is used for specifying key value, e.g., `PTM{!} = NQR`

Currently, there are 35 symbols available for specifying ptms, which are
```
~ ! @ $ % ^ & * ( ) _ + ` - | \ : " ; ' < > ? . / 1 2 3 4 5 6 7 8 9 0
```
Please don't use these reserved symbols: `{ } # [ ] = ,`

Neutral loss can be specified by `PTM{1to2}`, e.g., `PTM{>to|} = ST`. If symbol2 is nothing, it can be specified by `PTM{1to}`, e.g. `PTM{>to} = ST`.

#### Generate Reverse Sequences
```
python sipros_prepare_protein_database.py -i original_database_file -o output_database_file -c config_file
```
The step will generate a new database file with reverse sequences. Update the path of `FASTA_Database` in the configuration file.

#### <a name="labelds"></a>Running The Database-searching

There are two basic versions of the database-searching: one for running on a single machine and another for running with MPI on a cluster.  

* __Single Machine Version:__ This version of the assembler should be used if you are going to run the database-searching on a single machine with one or more cores. The searching is invoked through a run `Sipros_OpemMP` in `bin` directory. The quick start command as shown below will be used in a batch job submission script or directly typed on the command line terminal.   

```
#!/bin/bash

# Single MS2 file
Sipros_OpemMP -o output_dir -f ms_data -c SiprosConfig.cfg

# Multiple MS2 files in a working directory
Sipros_OpemMP -o output_dir -w workingdirectory -c SiprosConfig.cfg

```
Results (`.Spe2Pep` files) will be saved on the output directory. if you have many configure files, specify `-g`, like `Sipros_OpemMP -o output_dir -w workingdirectory -g configurefiledirectory`. Use `./Sipros_OpemMP -h` for help information. 

* __MPI Version:__ This version of the database-searching should be used if you are going to run on a cluster with MPI support. The run script to invoke `Sipros_MPI` depends on the cluster management and job scheduling system. An example bash script `submit_job.pbs` is provide in `configs` directory.
 
The quick start commands are:
```
### MPI Verion 
Sipros_MPI -o output_dir -w workingdirectory -c SiprosConfig.cfg

```
Results (`.Spe2Pep` files) will be saved on the output directory. if you have many configure files, specify `-g`, like `Sipros_MPI -o output_dir -w workingdirectory -g configurefiledirectory`.

### Guide to Regular Search

#### Create configure file

Please refer to [Configure File Setting](#config) for technical details. An example is available at [SiprosConfig.cfg](SiprosConfig.cfg).

#### Run Sipros

Please refer to [Running The Database-searching](#labelds).

#### Post-processing

The current version of scripts has been tested using Python 2.7.2, so if you are using different versions of Python (2.6.X or 3.X), you are encouraged to try with Python 2.7.2.

```
#!/bin/bash

cd Scripts
runSiprosFiltering.sh -in Spe2Pep_dir -o workingdirectory -c SiprosConfig.cfg

```

This step will generate related `tab`, `psm.txt`, `pep.txt`, `pro.txt`, `pro2pep.txt`, and `pro2psm.txt` files. Please see the [OUTPUT.md](OUTPUT.md) file for description of the output files.

### Miscellany

#### <a name="PythonSetup"></a>Python Environment Setup

It is recommended to use `anaconda` to setup necessary Python libraries. Take Linux based system as an example:

1. Download `anaconda` at [https://www.continuum.io/downloads](https://www.continuum.io/downloads).

2. In your terminal window type the following instructions:
```
bash Anaconda2-4.3.1-Linux-x86_64.sh
```
3. Create an environment named sipros-env and activate the new environment to use it:
```
conda create --prefix ~/sipros-env
source activate ~/sipros-env
```

4. Install a new package (numpy, scipy, scikit-learn, lxml) in a this environment (~/sipros-env):
```
conda install --prefix ~/sipros-env numpy
```

```
conda install --prefix ~/sipros-env scipy
```

```
conda install --prefix ~/sipros-env scikit-learn
```

```
conda install --prefix ~/sipros-env lxml
```

5. You are good to go.

### Questions?

* [Xuan Guo](mailto:xuan_guo@outlook.com)
* [Chongle Pan](mailto:chongle.pan@gmail.com)
