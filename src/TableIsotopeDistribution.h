/*
 * TableIsotopeDistribution.h
 *
 *  Created on: Oct 1, 2015
 *      Author: xgo
 */

#ifndef TABLEISOTOPEDISTRIBUTION_H_
#define TABLEISOTOPEDISTRIBUTION_H_

#include "proNovoConfig.h"

class IsotopeDistribution;

using namespace std;

typedef long long INT64;

#define CLOCKSTART double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) << " Seconds." << endl << endl;

#define CURTIME time_t curr=time(0);cout << "current time is: " << ctime(&curr) <<endl;
#define CURMEM INT64 mem = checkMemoryUsage();cout<<"Memory used: " << mem <<  " MB."<< endl;

// To keep time information of functions.
#define MEMORYSTART INT64 mem_start = checkMemoryUsage(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define MEMORYSTOP INT64 mem_end = checkMemoryUsage(); cout << "Function " << __FUNCTION__ << "() finished . " << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl;
// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage() {
	// get KB memory into count
	unsigned int count = 0;

#if defined(__linux__)
	ifstream f("/proc/self/status"); // read the linux file
	while (!f.eof()) {
		string key;
		f >> key;
		if (key == "VmData:") {     // size of data
			f >> count;
			break;
		}

	}
	f.close();
#endif

	// return MBs memory (size of data)
	return (count / 1024);
}
;

IsotopeDistribution * sum(const IsotopeDistribution * distribution0, const IsotopeDistribution * distribution1, IsotopeDistribution * sumDistribution);

class TableIsotopeDistribution {
private:
	vector<vector<int> > vviIndexDictionary;
	map<char, int> mResidueToIndex;
	map<int, string> mIndexToResidue;
	int iNumberOfUniqueResidue;
	int iNumberForPrecalculatedResidue;
public:
	TableIsotopeDistribution();
	~TableIsotopeDistribution();

	vector<vector<IsotopeDistribution*> *> vpIsotopeDistributionNterm;
	vector<vector<IsotopeDistribution*> *> vpIsotopeDistributionCterm;

	bool computeProductIon(string _sSequence, vector<vector<double> > & _vvdYionMass, vector<vector<double> > & _vvdYionProb, vector<vector<double> > & _vvdBionMass, vector<vector<double> > & _vvdBionProb);
	IsotopeDistribution* getDistributionByGivenResidues(const string & _sSequence, const int & _positionBegin, const int & _positionEnd);
	int getIndexByGivenResidues(const string & _sSequence, const int & _positionBegin, const int & _positionEnd);
	int getIndexByGivenResidues(const vector<int> & _viSequence);
	int getIndexByGivenResidues(const vector<int> & _viSequence, const int & _positionBegin, const int & _positionEnd);
	void decodeIndexToResidues(const int & _index, const int & _size, const int & _uniqueResidue, vector<int> & _viSequence);
	void initializeEverything(const int & _numberOfPrecalculatedResidues);
	void initializeIndexDictionary(const int & _nResidues, const int & _nMemoryResidues);
	void precalculateIsotopicDistribution(const vector<int> & _viSequence, IsotopeDistribution * _sumDistribution);
	bool setResidueIndex();
	void setVectorSize(const int & _nResidues, const int & _nPrecalculatedResidues);

	void startPrecalculation();

	void debug();

};

#endif /* TABLEISOTOPEDISTRIBUTION_H_ */

/*
 * Parse command line arguments
 * Populate vsFT2Filenames
 * Set up SiprosConfig
 */

/*void initializeArguments2(int argc, char **argv, vector<string> & vsFT2Filenames, string & sWorkingDirectory, string & sConfigFilename, string & sSingleWorkingFile, string & sOutputDirectory, bool & bScreenOutput) {
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
			cout << "Usage: -w WorkingDirectory -c ConfigurationFile, -f: A single MS2 or FT2 file to be processed" << endl;
			cout << "If configuration file is not specified, Sipros will look for SiprosConfig.cfg in the directory of FT2 files" << endl;
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
	//searchFT2Files(vsFT2Filenames, sWorkingDirectory, bScreenOutput);
	if ((sOutputDirectory == "") && (sWorkingDirectory != ""))
		sOutputDirectory = sWorkingDirectory;

};

void print(vector<vector<double> > & _v) {
	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(6);
	for (size_t i = 0; i < _v.size(); i++) {
		for (size_t j = 0; j < _v.at(i).size(); j++) {
			cout << _v.at(i).at(j) << endl;
		}
	}
};


int xxx(int argc, char **argv) {

	vector<string> vsFT2Filenames;
	bool bScreenOutput;
	string sWorkingDirectory, sConfigFilename, sSingleWorkingFile, sOutputDirectory;
	initializeArguments2(argc, argv, vsFT2Filenames, sWorkingDirectory, sConfigFilename, sSingleWorkingFile, sOutputDirectory, bScreenOutput);
	// Load config file
	if (!ProNovoConfig::setFilename(sConfigFilename)) {
		cerr << "Could not load config file " << sConfigFilename << endl;
		return 0;
	}
	TableIsotopeDistribution * tid = new TableIsotopeDistribution();
	{
		//CLOCKSTART
		//MEMORYSTART
		tid->initializeEverything(4);
		//MEMORYSTOP
		//CLOCKSTOP
	}
	{
		//CLOCKSTART
		//MEMORYSTART
		//tid->debug();
		tid->startPrecalculation();
		//MEMORYSTOP
		//CLOCKSTOP
	}
	string _sequence = "[HGDIEYR]";
	vector<vector<double> > _vvdYionMass;
	vector<vector<double> > _vvdYionProb;
	vector<vector<double> > _vvdBionMass;
	vector<vector<double> > _vvdBionProb;
	tid->computeProductIon(_sequence, _vvdYionMass, _vvdYionProb, _vvdBionMass, _vvdBionProb);
	cout << "BMass" << endl;
	print(_vvdBionMass);
	cout << "BProb" << endl;
	print(_vvdBionProb);
	cout << "YMass" << endl;
	print(_vvdYionMass);
	cout << "YMass" << endl;
	print(_vvdBionProb);
	delete tid;
	cout << "Test" << endl;

}*/
