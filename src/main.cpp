#include <iostream>
#include <string>
#include <cstdlib>

#include "directoryStructure.h"
#include "proNovoConfig.h"
#include "ms2scanvector.h"
// #include <gperftools/profiler.h>

using namespace std;

void searchFT2Files(vector<string> & vsFT2Filenames, const string & sWorkingDirectory, bool bScreenOutput) {
	int i, iFileNum;
	DirectoryStructure working_dir(sWorkingDirectory);
	working_dir.setPattern(".ft2");
	working_dir.getFiles(vsFT2Filenames);
	working_dir.setPattern(".FT2");
	working_dir.getFiles(vsFT2Filenames);
	working_dir.setPattern(".ms2");
	working_dir.getFiles(vsFT2Filenames);
	working_dir.setPattern(".MS2");
	working_dir.getFiles(vsFT2Filenames);
	iFileNum = (int) vsFT2Filenames.size();
	if (iFileNum == 0) {
		cerr << "no scan file in the working directory" << endl;
		exit(1);
	}
	for (i = 0; i < iFileNum; i++)
		vsFT2Filenames.at(i) = sWorkingDirectory + ProNovoConfig::getSeparator() + vsFT2Filenames.at(i);
}

void searchConfigureFiles(vector<string> & vsConfigureFilenames, const string & sConfigFileDirectory,
		bool bScreenOutput) {
	int i, iFileNum;
	DirectoryStructure working_dir(sConfigFileDirectory);
	working_dir.setPattern(".cfg");
	working_dir.getFiles(vsConfigureFilenames);
	working_dir.setPattern(".CFG");
	working_dir.getFiles(vsConfigureFilenames);

	iFileNum = (int) vsConfigureFilenames.size();
	if (iFileNum == 0) {
		cerr << "no configure file in the directory" << endl;
		exit(1);
	}
	for (i = 0; i < iFileNum; i++)
		vsConfigureFilenames.at(i) = sConfigFileDirectory + ProNovoConfig::getSeparator() + vsConfigureFilenames.at(i);
}

/* 
 * Parse command line arguments
 * Populate vsFT2Filenames
 * Set up SiprosConfig
 */

void initializeArguments(int argc, char **argv, vector<string> & vsFT2Filenames, string & sWorkingDirectory,
		vector<string> & vsConfigureFilenames, string & sSingleWorkingFile, string & sOutputDirectory,
		bool & bScreenOutput)
// Under MPI mode, a user can specify configure file directory by specifying -g
// If no configre file or directory is specified, configure file directory is working directory
		{
	int i;

	string sConfigFileDirectory, sConfigFilename;
	// Grab command line arguments
	vector<string> vsArguments;

	sWorkingDirectory = "";
	sConfigFilename = "";
	sSingleWorkingFile = "";
	sOutputDirectory = "";
	sConfigFileDirectory = "";
	bScreenOutput = true;

	while (argc--)
		vsArguments.push_back(*argv++);
	for (i = 1; i <= (int) vsArguments.size() - 1; i++)
		if (vsArguments[i] == "-w")
			sWorkingDirectory = vsArguments[++i];
		else if (vsArguments[i] == "-c")
			sConfigFilename = vsArguments[++i];
		else if (vsArguments[i] == "-f")
			sSingleWorkingFile = vsArguments[++i];
		else if (vsArguments[i] == "-o")
			sOutputDirectory = vsArguments[++i];
		else if (vsArguments[i] == "-g")
			sConfigFileDirectory = vsArguments[++i];
		else if (vsArguments[i] == "-s")
			bScreenOutput = false;
		else if ((vsArguments[i] == "-h") || (vsArguments[i] == "--help")) {
			cout << "Usage: -w WorkingDirectory -c ConfigurationFile, -f: A single MS2 or FT2 file to be processed"
					<< endl;
			cout
					<< "If configuration file is not specified, Sipros will look for SiprosConfig.cfg in the directory of FT2 files"
					<< endl;
			cout << "-o output directory. If not specified, it is the same as that of the input scan file," << endl;
			cout << "-g configure file directory. -s silence all standard output." << endl;
			exit(0);
		} else {
			cerr << "Unknown option " << vsArguments[i] << endl << endl;
			exit(1);
		}
	if ((sWorkingDirectory == "") && (sSingleWorkingFile == ""))
		sWorkingDirectory = ".";

	if ((sWorkingDirectory != "") && (sSingleWorkingFile != "")) {
		cerr << "Either a input scan file or the directory of input scan files needs to be specified" << endl;
		exit(1);
	}
	if ((sConfigFilename == "") && (sConfigFileDirectory == ""))
		//sConfigFilename = sWorkingDirectory + ProNovoConfig::getSeparator() + "SiprosConfig.cfg";
		// Without specifying configure file and configure file directory, the default configure file directory
		// is working directory. In the openmp only version, the default configure file is SiprosConfig.cfg
		sConfigFileDirectory = sWorkingDirectory;

	if (sConfigFileDirectory != "")
		searchConfigureFiles(vsConfigureFilenames, sConfigFileDirectory, bScreenOutput);
	else
		vsConfigureFilenames.push_back(sConfigFilename);

	if (sSingleWorkingFile != "")
		vsFT2Filenames.push_back(sSingleWorkingFile);
	else
		searchFT2Files(vsFT2Filenames, sWorkingDirectory, bScreenOutput);
	if ((sOutputDirectory == "") && (sWorkingDirectory != ""))
		sOutputDirectory = sWorkingDirectory;

}

void handleScan(const string & sFT2filename, const string & sOutputDirectory, const string & sConfigFilename,
		bool bScreenOutput) {
	MS2ScanVector * pMainMS2ScanVector = new MS2ScanVector(sFT2filename, sOutputDirectory, sConfigFilename,
			bScreenOutput);

	if (bScreenOutput) {
		cout << "Reading MS2 scan file " << sFT2filename << endl;
	}

	if (!pMainMS2ScanVector->loadFile()) {
		cerr << "Error: Failed to load file: " << sFT2filename << endl;
	} else {
		// search all MS2 scans and write output to a file
		pMainMS2ScanVector->startProcessing();
	}

	delete pMainMS2ScanVector; //free memory of vpAllMS2Scans
}

int main(int argc, char **argv) {
	// A list of FT2/MS2 files to be searched
	double begin = omp_get_wtime();
	vector<string> vsFT2Filenames;
	// A list of configure files
	vector<string> vsConfigureFilenames;
	vsConfigureFilenames.clear();
	bool bScreenOutput;
	string sWorkingDirectory, sConfigFilename, sSingleWorkingFile, sOutputDirectory;
	initializeArguments(argc, argv, vsFT2Filenames, sWorkingDirectory, vsConfigureFilenames, sSingleWorkingFile,
			sOutputDirectory, bScreenOutput);
	if (vsConfigureFilenames.empty()) {
		vsConfigureFilenames.push_back(sConfigFilename);
	}
	// ProfilerStart("low_res.log");
	for (size_t j = 0; j < vsConfigureFilenames.size(); j++) {
		// Load config file
		if (!ProNovoConfig::setFilename(vsConfigureFilenames.at(j))) {
			cerr << "Could not load config file " << sConfigFilename << endl;
			return 0;
		}

		// Process one FT2 file at a time
		for (size_t i = 0; i < vsFT2Filenames.size(); i++) {
			cout << "MS file name:\t" << vsFT2Filenames.at(i) << "\nCfg:\t" << vsConfigureFilenames.at(j) << endl;
			handleScan(vsFT2Filenames.at(i), sOutputDirectory, vsConfigureFilenames.at(j), bScreenOutput);
		}
	}
	// ProfilerStop();

	double end = omp_get_wtime();
	cout << endl;
	cout << "Total time:\t" << double(end - begin) << " Seconds." << endl << endl;
	return 0;
}
