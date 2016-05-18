#include <iostream>
#include <string>
#include <cstdlib>

#include "directoryStructure.h"
#include "proNovoConfig.h"
#include "ms2scanvector.h"
#include <gperftools/profiler.h>

using namespace std;

void debug();

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

/* 
 * Parse command line arguments
 * Populate vsFT2Filenames
 * Set up SiprosConfig
 */

void initializeArguments(int argc, char **argv, vector<string> & vsFT2Filenames, string & sWorkingDirectory,
		string & sConfigFilename, string & sSingleWorkingFile, string & sOutputDirectory, bool & bScreenOutput) {
	int i;
	// Grab command line arguments
	vector<string> vsArguments;

	sWorkingDirectory = "";
	sConfigFilename = "";
	sSingleWorkingFile = "";
	sOutputDirectory = "";
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
		else if (vsArguments[i] == "-s")
			bScreenOutput = false;
		else if ((vsArguments[i] == "-h") || (vsArguments[i] == "--help")) {
			cout << "Usage: -w WorkingDirectory -c ConfigurationFile, -f: A single MS2 or FT2 file to be processed"
					<< endl;
			cout
					<< "If configuration file is not specified, Sipros will look for SiprosConfig.cfg in the directory of FT2 files"
					<< endl;
			cout << "-o output directory. If not specified, it is the same as that of the input scan file," << endl;
			cout << "-s silence all standard output." << endl;
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
	if (sConfigFilename == "")
		sConfigFilename = sWorkingDirectory + ProNovoConfig::getSeparator() + "SiprosConfig.cfg";
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

	if (bScreenOutput)
		cout << "Reading MS2 scan file " << sFT2filename << endl;

	if (!pMainMS2ScanVector->loadFT2file())
		cerr << "Error: Failed to load file: " << sFT2filename << endl;
	else {
		// search all MS2 scans and write output to a file
		pMainMS2ScanVector->startProcessing();
	}
	delete pMainMS2ScanVector; //free memory of vpAllMS2Scans
}

int main(int argc, char **argv) {
	// A list of FT2/MS2 files to be searched
	double begin = omp_get_wtime();
	vector<string> vsFT2Filenames;
	bool bScreenOutput;
	string sWorkingDirectory, sConfigFilename, sSingleWorkingFile, sOutputDirectory;
	initializeArguments(argc, argv, vsFT2Filenames, sWorkingDirectory, sConfigFilename, sSingleWorkingFile,
			sOutputDirectory, bScreenOutput);
	// Load config file
	if (!ProNovoConfig::setFilename(sConfigFilename)) {
		cerr << "Could not load config file " << sConfigFilename << endl;
		return 0;
	}

	if (0) {
		ProfilerStart("Sipros5_openmp.prof");
		debug();
		ProfilerStop();
		return 0;
	}

	// Process one FT2 file at a time

	for (int i = 0; i < (int) vsFT2Filenames.size(); i++) {
		handleScan(vsFT2Filenames[i], sOutputDirectory, sConfigFilename, bScreenOutput);
	}
	double end = omp_get_wtime();
	cout << endl;
	cout << "Total time:\t" << double(end - begin) << " Seconds." << endl << endl;
	return 0;
}

void debug() {
	string xxx = "N!CYDDFESSIGTTSQHM~RYR";
	string zzz = "CHONPS";
	vector<int> viComposition(6);
	map<string, vector<int> >::iterator ResidueIter;
	for (size_t i = 0; i < xxx.length(); i++) {
		ResidueIter = ProNovoConfig::configIsotopologue.mResidueAtomicComposition.find(xxx.substr(i, 1));
		if (ResidueIter == ProNovoConfig::configIsotopologue.mResidueAtomicComposition.end()) {
			cout << "error" << endl;
		} else {
			for (size_t j = 0; j < 6; j++) {
				viComposition.at(j) += ResidueIter->second.at(j);
			}
		}

	}
	ResidueIter = ProNovoConfig::configIsotopologue.mResidueAtomicComposition.find("Nterm");
	if (ResidueIter == ProNovoConfig::configIsotopologue.mResidueAtomicComposition.end()) {
		cout << "error" << endl;
	} else {
		for (size_t j = 0; j < 6; j++) {
			viComposition.at(j) += ResidueIter->second.at(j);
		}
	}
	ResidueIter = ProNovoConfig::configIsotopologue.mResidueAtomicComposition.find("Cterm");
	if (ResidueIter == ProNovoConfig::configIsotopologue.mResidueAtomicComposition.end()) {
		cout << "error" << endl;
	} else {
		for (size_t j = 0; j < 6; j++) {
			viComposition.at(j) += ResidueIter->second.at(j);
		}
	}
	for (size_t i = 0; i < 6; i++) {
		cout << zzz.at(i) << viComposition.at(i);
	}
	cout << endl;
	IsotopeDistribution yyy;
	double begin = omp_get_wtime();
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < 1000000; i++) {
		ProNovoConfig::configIsotopologue.computeIsotopicDistribution2(xxx, yyy);
	}
	double end = omp_get_wtime();
	cout << "1:\t" << double(end - begin) << " Seconds." << endl << endl;
	begin = omp_get_wtime();
	for (int i = 0; i < 1000000; i++) {
		ProNovoConfig::configIsotopologue.computeIsotopicDistribution3(xxx, yyy);
	}
	end = omp_get_wtime();
	cout << "2:\t" << double(end - begin) << " Seconds." << endl << endl;
	//yyy.print();
}
