#ifndef MS2SCANVECTOR_H
#define MS2SCANVECTOR_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "ms2scan.h"
#include "tokenvector.h"
#include "proNovoConfig.h"
#include "proteindatabase.h"
#include "./Scores/CometSearch.h"
#include "./Scores/MVH.h"
#include "ProteinDbParser.h"

#define ZERO            0.00000001
#define PEPTIDE_ARRAY_SIZE  10000
#define TASKWAIT_SIZE 300

using namespace std;

class MS2ScanVector {
	// All MS2 scans to be scored
	// the MS2 scans are sorted by their precursor masses
	vector<MS2Scan *> vpAllMS2Scans;
	// the precursor mass of these MS2 scans
	// this is used for quickly inserting a peptide into correct MS2 scans
	// vpAllMS2Scans and vpPrecursorMasses are in the same order
	vector<double> vpPrecursorMasses;
	string sFT2Filename;    // the FT2 filename
	string sOutputFile; // the output file name
	string sConfigFile; // the configure file name
	// mass except N and C termini;
	map<char, double> mapResidueMass;
	// if true, allows standard output
	bool bScreenOutput;

	// find every MS2 scan whose precursor mass matches peptide mass
	bool assignPeptides2Scans(Peptide * currentPeptide);
	//return true if parent_charge should be 1
	bool ChargeDetermination(const vector<double> & vdAllmz, double pmz);
	bool isMS1HighRes(const string & target);
	static bool mygreater(double i, double j);
	static bool myless(MS2Scan * pMS2Scan1, MS2Scan * pMS2Scan2);
	static bool mylessScanId(MS2Scan * pMS2Scan1, MS2Scan * pMS2Scan2);
	// postprocessing all MS2 scans' results by multi-threading
	void postProcessAllMS2();
	// preprocessing all MS2 scans by multi-threading
	void preProcessAllMS2();
	void processPeptideArray(vector<Peptide*>& vpPeptideArray);
	// task version of processing peptides
	void processPeptideArrayTask(vector<Peptide*>& vpPeptideArray, omp_lock_t * pLck);
	void saveScan(MS2Scan * pMS2Scan);
	// search all MS2 scans against the protein list by multi-threading
	void searchDatabase();
	void searchDatabaseSnp();
	// task version of search database
	void searchDatabaseTask();
	pair<int, int> GetRangeFromMass(double lb, double ub);
	// write results to a SIP file
	// the SIP file will the same base file name as sFT2Filename
	// change the extension filename to ".SIP"
	void setOutputFile(const string & sFT2FilenameInput, const string & sOutputDirectory);
	void writeOutput();
	void writeOutputMultiScoresPin();
	void writeOutputMultiScoresSip();
	void writeOutputMultiScoresSpectrum2MutiPep();
	//calculate mean and standard deviation of scores of a ms2 scan
	void calculateMeanAndDeviation(int inumberScore, double dScoreSum, double dScoreSquareSum, double & dMean,
			double & dDeviation);
	void GetAllRangeFromMass(double dPeptideMass, vector<pair<int, int> > & vpPeptideMassRanges);
	string ParsePath(string sPath);

public:
	MS2ScanVector(const string & sFT2FilenameInput, const string & sOutputDirectory, const string & sConfigFilename,
			bool bScreenOutput);
	~MS2ScanVector();

	// Populate vpAllMS2Scans from the input FT2 file
	// Determine the charge state of every scan by calling the function MS2Scan::isSinglyCharged().
	// Create a +2 scan and a +3 scan for an unknown multiple charged scan.
	// Return false if there is a problem with the file
	bool loadFile();
	bool loadFT2file();
	bool loadMs2File();
	//Read FT2 files
	bool ReadFT2File();
	//Read Ms2 files
	bool ReadMs2File();
	// start functions to process the loaded FT2 file
	void startProcessing();

	size_t iMaxNumProteins;

	// variables for the MVH thread
	vector<double> ** _ppdAAforward;
	vector<double> ** _ppdAAreverse;
	vector<double> ** psequenceIonMasses;
	vector<char> ** pSeqs;
	int num_max_threads;
	void preMvh();
	void postMvh();
};

#endif // MS2SCANVECTOR_H
