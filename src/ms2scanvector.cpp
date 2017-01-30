#include "ms2scanvector.h"

MS2ScanVector::MS2ScanVector(const string & sFT2FilenameInput, const string & sOutputDirectory,
		const string & sConfigFilename, bool bScreenOutput) {
	unsigned int n;
	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;
	sFT2Filename = sFT2FilenameInput;
	sConfigFile = sConfigFilename;
	//mass_w = ProNovoConfig::getParentMassWindows();
	setOutputFile(sFT2FilenameInput, sOutputDirectory);
	this->bScreenOutput = bScreenOutput;
	for (n = 0; n < vsSingleResidueNames.size(); ++n) {
		mapResidueMass[vsSingleResidueNames[n][0]] = vdSingleResidueMasses[n];
	}
	iMaxNumProteins = 0;
}

MS2ScanVector::~MS2ScanVector() {
	// the destructor will free memory from vpAllMS2Scans
	for (size_t i = 0; i < vpAllMS2Scans.size(); i++) {
		delete vpAllMS2Scans.at(i);
	}
	/*vector<MS2Scan *>::iterator it;
	 for (it = vpAllMS2Scans.begin(); it < vpAllMS2Scans.end(); it++){
	 delete *it;
	 }*/
	vpAllMS2Scans.clear();
	vpPrecursorMasses.clear();
}

void MS2ScanVector::setOutputFile(const string& sFT2FilenameInput, const string& sOutputDirectory) {
	string sLeftFileName, sMiddleFileName, sRightFileName = ".sip";
	size_t pos;
	//cout<<ProNovoConfig::getSearchName()<<endl;
	if ((ProNovoConfig::getSearchName() == "") || (ProNovoConfig::getSearchName() == "Null")
			|| (ProNovoConfig::getSearchName() == "NULL") || (ProNovoConfig::getSearchName() == "null"))
		sRightFileName = ".sip";
	else
		sRightFileName = "." + ProNovoConfig::getSearchName() + ".sip";

	sLeftFileName = sFT2FilenameInput.substr(0, sFT2FilenameInput.length() - 4);
	if (sOutputDirectory == "")
		// if output directory is not specified, it is working directory by default
		sOutputFile = sLeftFileName + sRightFileName;
	else {
		pos = sLeftFileName.rfind(ProNovoConfig::getSeparator());
		if (pos == string::npos)
			sMiddleFileName = sLeftFileName;
		else
			sMiddleFileName = sLeftFileName.substr(pos + 1);
		sOutputFile = sOutputDirectory + ProNovoConfig::getSeparator() + sMiddleFileName + sRightFileName;
	}
}

bool MS2ScanVector::ReadFT2File() {
	bool bReVal, flag_1stScan = true; //flag_1stScan true indicates pMS2Scan is empty
	string sline;
	istringstream input;
	MS2Scan * pMS2Scan = NULL;
	ifstream ft2_stream(sFT2Filename.c_str());
	int tmp_charge;
	double tmp_mz, tmp_intensity;

	bReVal = ft2_stream.is_open();
	if (bReVal) {
		while (!ft2_stream.eof()) {
			sline.clear();
			// string tmp_sline = sline;
			getline(ft2_stream, sline);
			// cout << sline <<"wyf"<<endl;
			if (sline == "")
				continue;
			if ((sline.at(0) >= '0') && (sline.at(0) <= '9')) {
				TokenVector words(sline, " \t\n\r");
				// vector<string> words;
				// TokenVector::Parse(sline, " \r\t\n", words);
				if (words.size() < 6) {
					// The judgement of MS scan resolution here will be overwritten,
					// we still keep the related codes here for future possible usuage
					pMS2Scan->isMS2HighRes = false;
					pMS2Scan->viCharge.push_back(0);
				} else {
					pMS2Scan->isMS2HighRes = true;
					input.clear();
					input.str(words[5]);
					input >> tmp_charge;
					input.clear();
					pMS2Scan->viCharge.push_back(tmp_charge);
				}
				input.clear();
				input.str(words[0]);
				input >> tmp_mz;
				input.clear();
				pMS2Scan->vdMZ.push_back(tmp_mz);
				input.str(words[1]);
				input >> tmp_intensity;
				input.clear();
				pMS2Scan->vdIntensity.push_back(tmp_intensity);
				words.clear();
			} else if (sline.at(0) == 'S') {
				if (flag_1stScan)
					flag_1stScan = false;
				else if (pMS2Scan->vdIntensity.empty())
					delete pMS2Scan;
				else {
					if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
						pMS2Scan->isMS2HighRes = true;
					else
						pMS2Scan->isMS2HighRes = false;
					saveScan(pMS2Scan);
				}
				pMS2Scan = new MS2Scan;
				pMS2Scan->sFT2Filename = sFT2Filename;
				// vector<string> words;
				// TokenVector::Parse(sline, " \r\t\n", words);
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words.at(1));
				input >> pMS2Scan->iScanId;
				input.clear();
				input.str(words.at(3));
				input >> pMS2Scan->dParentMZ;
				input.clear();
				pMS2Scan->isMS1HighRes = isMS1HighRes(words.at(3));
				// in case no Z Line
				pMS2Scan->iParentChargeState = 0;
				pMS2Scan->dParentNeutralMass = 0;
				words.clear();
			} else if (sline.at(0) == 'Z') {
				// vector<string> words;
				// TokenVector::Parse(sline, " \r\t\n", words);
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan->iParentChargeState;
				input.clear();
				//input.str(words[2]);
				//input >> pMS2Scan->dParentNeutralMass;
				//input.clear();
				words.clear();
			} else if (sline.at(0) == 'I') {
				// vector<string> words;
				// TokenVector::Parse(sline, " \r\t\n", words);
				TokenVector words(sline, " \r\t\n");
				if (words.at(1) == "ScanType") {
					pMS2Scan->setScanType(words.at(2) + words.at(5) + words.at(6));
				}
				if (words.at(1) == "RetentionTime") {
					pMS2Scan->setRTime(words.at(2));
				}
				words.clear();
			}
		}
		ft2_stream.clear();
		ft2_stream.close();

		//recognition high or low MS2
//	cout<<ProNovoConfig::getMassAccuracyFragmentIon()<<endl;
		if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1)
			pMS2Scan->isMS2HighRes = true;
		else
			pMS2Scan->isMS2HighRes = false;
		if (!flag_1stScan) // To avoid empty file
			saveScan(pMS2Scan);
	}
	return bReVal;
}

bool MS2ScanVector::ReadMs2File() {
	bool bReVal, flag_1stScan = true; //flag_1stScan true indicates pMS2Scan is empty
	string sline;
	istringstream input;
	MS2Scan * pMS2Scan = NULL;
	ifstream ft2_stream(sFT2Filename.c_str());
	int tmp_charge;
	double tmp_mz, tmp_intensity;

	bReVal = ft2_stream.is_open();
	if (bReVal) {
		while (!ft2_stream.eof()) {
			sline.clear();
			// string tmp_sline = sline;
			getline(ft2_stream, sline);
			// cout << sline <<"wyf"<<endl;
			if (sline == "") {
				continue;
			}
			if ((sline.at(0) >= '0') && (sline.at(0) <= '9')) {
				TokenVector words(sline, " \t\n\r");
				input.clear();
				input.str(words[0]);
				input >> tmp_mz;
				input.clear();
				pMS2Scan->vdMZ.push_back(tmp_mz);
				input.str(words[1]);
				input >> tmp_intensity;
				input.clear();
				pMS2Scan->vdIntensity.push_back(tmp_intensity);
				pMS2Scan->viCharge.push_back(0);
				/*if (words.size() >= 2) {
					input.clear();
					input.str(words[2]);
					input >> tmp_charge;
					input.clear();
					pMS2Scan->viCharge.push_back(tmp_charge);
				} else {
					pMS2Scan->viCharge.push_back(0);
				}*/
			} else if (sline.at(0) == 'S') {
				if (flag_1stScan) {
					flag_1stScan = false;
				} else if (pMS2Scan->vdIntensity.empty()) {
					delete pMS2Scan;
				} else {
					if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1) {
						pMS2Scan->isMS2HighRes = true;
					} else {
						pMS2Scan->isMS2HighRes = false;
					}
					saveScan(pMS2Scan);
				}
				pMS2Scan = new MS2Scan;
				pMS2Scan->sFT2Filename = sFT2Filename;
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan->iScanId;
				input.clear();
				input.str(words[3]);
				input >> pMS2Scan->dParentMZ;
				input.clear();
				pMS2Scan->isMS1HighRes = isMS1HighRes(words[3]);
				// in case no Z Line
				pMS2Scan->iParentChargeState = 0;
				pMS2Scan->dParentNeutralMass = 0;
			} else if (sline.at(0) == 'Z') {
				TokenVector words(sline, " \r\t\n");
				input.clear();
				input.str(words[1]);
				input >> pMS2Scan->iParentChargeState;
				input.clear();
				//input.str(words[2]);
				//input >> pMS2Scan->dParentNeutralMass;
				//input.clear();
			} else if (sline.at(0) == 'I') {
				TokenVector words(sline, " \r\t\n");
				if (words[1] == "ScanType") {
					pMS2Scan->setScanType(words[2] + words[5] + words[6]);
				}
				if (words[1] == "RetTime") {
					pMS2Scan->setRTime(words[2]);
				}
			}
		}
		ft2_stream.clear();
		ft2_stream.close();

		// recognition high or low MS2
		// cout<<ProNovoConfig::getMassAccuracyFragmentIon()<<endl;
		if (ProNovoConfig::getMassAccuracyFragmentIon() < 0.1) {
			pMS2Scan->isMS2HighRes = true;
		} else {
			pMS2Scan->isMS2HighRes = false;
		}
		if (!flag_1stScan) { // To avoid empty file
			saveScan(pMS2Scan);
		}
	}
	return bReVal;
}

bool matchPattern(const string & sPattern, const string & sFilename) {
	// if the filename ends with sPattern
	// then it is matched
	if ((sFilename.compare(sFilename.length() - sPattern.length(), sPattern.length(), sPattern) == 0))
		return true;
	else
		return false;

}

bool MS2ScanVector::loadFile() {
	bool bReVal;
	if (matchPattern("ms2", sFT2Filename) || matchPattern("MS2", sFT2Filename)) {
		bReVal = loadMs2File();
	} else {
		bReVal = loadFT2file();
	}
	return bReVal;
}

bool MS2ScanVector::loadFT2file() {
	// read all MS2 scans from the file and populate vpAllMS2Scans
	// sort all MS2 scans in vpAllProteins by ascending order of their precursor masses
	// save their precursor masses in vpPrecursorMasses to quick look-up in assignPeptides2Scans()

	bool bReVal; //false when the file fails to be opened.
	vector<MS2Scan *>::iterator it;
	bReVal = ReadFT2File();
	if (bReVal) {
		sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), myless);
		for (it = vpAllMS2Scans.begin(); it < vpAllMS2Scans.end(); it++) {
			vpPrecursorMasses.push_back((*it)->dParentNeutralMass);
		}
	}

	return bReVal;
}

bool MS2ScanVector::loadMs2File() {
	// read all MS2 scans from the file and populate vpAllMS2Scans
	// sort all MS2 scans in vpAllProteins by ascending order of their precursor masses
	// save their precursor masses in vpPrecursorMasses to quick look-up in assignPeptides2Scans()

	bool bReVal; //false when the file fails to be opened.
	vector<MS2Scan *>::iterator it;
	bReVal = ReadMs2File();
	if (bReVal) {
		sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), myless);
		for (it = vpAllMS2Scans.begin(); it < vpAllMS2Scans.end(); it++)
			vpPrecursorMasses.push_back((*it)->dParentNeutralMass);
	}
	return bReVal;
}

bool MS2ScanVector::isMS1HighRes(const std::string& target) {
	bool bReVal = true;
	if (target.at(target.size() - 2) == '.')
		bReVal = false;
	return bReVal;
}

bool MS2ScanVector::ChargeDetermination(const std::vector<double>& vdAllmz, double pmz)
// return true, if the charge state is decided to be one

		{
	bool bReVal = false;
	int iPosi;
	vector<double> vdtempMZ = vdAllmz;
	iPosi = max(0, (int) (((int) vdtempMZ.size()) * 0.05 - 1));
	nth_element(vdtempMZ.begin(), vdtempMZ.begin() + iPosi, vdtempMZ.end(), mygreater);
	if (vdtempMZ.at(iPosi) < pmz)
		bReVal = true;
	return bReVal;
}

bool MS2ScanVector::mygreater(double i, double j) {
	return (i > j);
}

bool MS2ScanVector::myless(MS2Scan* pMS2Scan1, MS2Scan* pMS2Scan2) {
	return (pMS2Scan1->dParentNeutralMass < pMS2Scan2->dParentNeutralMass);
}

void MS2ScanVector::saveScan(MS2Scan * pMS2Scan) { //parentChargeState > 0, save, otherwise, try 1, or 2 and 3
	int j;
	bool bchargeOne;
	MS2Scan * pMS2newScan;
	if (pMS2Scan->iParentChargeState == 0) {
		bchargeOne = ChargeDetermination(pMS2Scan->vdMZ, pMS2Scan->dParentMZ);
		if (bchargeOne) {
			j = 1;
			pMS2newScan = new MS2Scan;
			*pMS2newScan = *pMS2Scan;
			pMS2newScan->iParentChargeState = j;
			pMS2newScan->dParentNeutralMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState)
					- pMS2newScan->iParentChargeState * ProNovoConfig::getProtonMass();
			pMS2newScan->dParentMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState);
			vpAllMS2Scans.push_back(pMS2newScan);
		} else {
			// if the charge state is not +1, it could be greater than +3,
			// but current we just try +2 and +3.
			for (j = 2; j <= 3; j++) {
				pMS2newScan = new MS2Scan;
				*pMS2newScan = *pMS2Scan;
				pMS2newScan->iParentChargeState = j;
				pMS2newScan->dParentNeutralMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState)
						- pMS2newScan->iParentChargeState * ProNovoConfig::getProtonMass();
				pMS2newScan->dParentMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState);
				vpAllMS2Scans.push_back(pMS2newScan);
			}
		}
	} else {
		pMS2newScan = new MS2Scan;
		*pMS2newScan = *pMS2Scan;
		pMS2newScan->dParentNeutralMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState)
				- (pMS2newScan->iParentChargeState) * ProNovoConfig::getProtonMass();
		pMS2newScan->dParentMass = (pMS2newScan->dParentMZ) * (pMS2newScan->iParentChargeState);
		vpAllMS2Scans.push_back(pMS2newScan);
	}
	if (pMS2newScan->dParentMass > ProNovoConfig::dMaxMS2ScanMass) {
		ProNovoConfig::dMaxMS2ScanMass = pMS2newScan->dParentMass;
	}
	if (pMS2newScan->iParentChargeState > ProNovoConfig::options.iMaxFragmentCharge) {
		ProNovoConfig::options.iMaxFragmentCharge = pMS2newScan->iParentChargeState;
	}
	delete pMS2Scan;
}

void MS2ScanVector::preProcessAllMS2() {
	//MEMORYSTART
	size_t i, iScanSize;
	iScanSize = vpAllMS2Scans.size();
	if (bScreenOutput) {
		cout << "Preprocessing " << vpAllMS2Scans.size() << " scans. " << endl;
	}

	if (ProNovoConfig::bXcorrEnable) {
#pragma omp parallel
		{
			//MH: Must be equal to largest possible array
			int iArraySize = (int) ((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0) * ProNovoConfig::dHighResInverseBinWidth);
			double *pdTmpRawData;
			double *pdTmpFastXcorrData;
			double *pdTmpCorrelationData;
			double *pdTmpSmoothedSpectrum;
			double *pdTmpPeakExtracted;
			try {
				pdTmpRawData = new double[iArraySize]();
				pdTmpFastXcorrData = new double[iArraySize]();
				pdTmpCorrelationData = new double[iArraySize]();
				pdTmpSmoothedSpectrum = new double[iArraySize]();
				pdTmpPeakExtracted = new double[iArraySize]();
			} catch (std::bad_alloc& ba) {
				char szErrorMsg[256];
				sprintf(szErrorMsg, " Error - new(pd[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
				sprintf(szErrorMsg + strlen(szErrorMsg),
						"Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
				sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
				string strErrorMsg(szErrorMsg);
				logerr(szErrorMsg);
				exit(1);
			}
#pragma omp for schedule(dynamic)
			for (size_t ii = 0; ii < iScanSize; ii++) {
				struct Query * pQuery = new Query();
				if (!CometSearch::Preprocess(pQuery, vpAllMS2Scans.at(ii), pdTmpRawData, pdTmpFastXcorrData,
						pdTmpCorrelationData, pdTmpSmoothedSpectrum, pdTmpPeakExtracted)) {
					cout << "Error" << endl;
					exit(1);
				}
				vpAllMS2Scans.at(ii)->pQuery = pQuery;
			}
			delete[] pdTmpRawData;
			delete[] pdTmpFastXcorrData;
			delete[] pdTmpCorrelationData;
			delete[] pdTmpSmoothedSpectrum;
			delete[] pdTmpPeakExtracted;
		}
	} //if (ProNovoConfig::bXcorrEnable)

	cout << "Rank " << ProNovoConfig::iRank << ":\tMvh preprocess begin" << endl;
	if (ProNovoConfig::bMvhEnable) {
		cout << "Rank " << ProNovoConfig::iRank << ":\tMvh preprocess sorting begin" << endl;
		{
#pragma omp for schedule(dynamic)
			for (size_t ii = 0; ii < iScanSize; ii++) {
				vpAllMS2Scans.at(ii)->sortPeakList();
			}
		}
		cout << "Rank " << ProNovoConfig::iRank << ":\tMvh preprocess sorting end" << endl;
		double minObservedMz = numeric_limits<double>::max();
		double maxObservedMz = 0;
		for (size_t ii = 0; ii < iScanSize; ii++) {
			minObservedMz = min(minObservedMz, vpAllMS2Scans.at(ii)->vdMZ.front());
			maxObservedMz = max(maxObservedMz, vpAllMS2Scans.at(ii)->vdMZ.back());
		}

#pragma omp parallel
		{
			multimap<double, double> * IntenSortedPeakPreData = new multimap<double, double>();
#pragma omp for schedule(dynamic)
			for (size_t ii = 0; ii < iScanSize; ii++) {
				vpAllMS2Scans.at(ii)->peakData = new map<double, char>();
				vpAllMS2Scans.at(ii)->peakData->clear();
				vpAllMS2Scans.at(ii)->intenClassCounts = new vector<int>();
				vpAllMS2Scans.at(ii)->intenClassCounts->clear();
				MVH::Preprocess(vpAllMS2Scans.at(ii), IntenSortedPeakPreData, minObservedMz, maxObservedMz);
			}
			delete IntenSortedPeakPreData;
			IntenSortedPeakPreData = NULL;
		}
		int maxPeakBins = 0;
		int iNumSkippedScans = 0;
		for (size_t ii = 0; ii < iScanSize; ii++) {
			if (vpAllMS2Scans.at(ii)->totalPeakBins > maxPeakBins) {
				maxPeakBins = vpAllMS2Scans.at(ii)->totalPeakBins;
			}
			if (vpAllMS2Scans.at(ii)->bSkip) {
				iNumSkippedScans++;
			}
		}
		cout << "Total Scans: " << iScanSize << "\t Skipped Scans: " << iNumSkippedScans << endl;
		cout << "Rank " << ProNovoConfig::iRank << ":\tMvh preprocess creating table begin" << endl;
		MVH::initialLnTable(((size_t) maxPeakBins));
		cout << "Rank " << ProNovoConfig::iRank << ":\tMvh preprocess creating table end" << endl;
	} // if (ProNovoConfig::bMvhEnable) {
	cout << "Rank " << ProNovoConfig::iRank << ":\tMvh preprocess end" << endl;

	cout << "Rank " << ProNovoConfig::iRank << ":\tScan preprocess begin" << endl;
#pragma omp parallel for schedule(guided)
	for (i = 0; i < iScanSize; i++) {
		vpAllMS2Scans.at(i)->preprocess();
	}
	cout << "Rank " << ProNovoConfig::iRank << ":\tScan preprocess end" << endl;

	cout << "Rank " << ProNovoConfig::iRank << ":\tMultiscore preprocess begin" << endl;
	if (ProNovoConfig::bMultiScores) {
		for (i = 0; i < iScanSize; i++) {
			if (vpAllMS2Scans.at(i)->iParentChargeState > ProNovoConfig::iMaxPercusorCharge) {
				ProNovoConfig::iMaxPercusorCharge = vpAllMS2Scans.at(i)->iParentChargeState;
			}
		}
		ProNovoConfig::TOP_N = 50;
#pragma omp parallel for schedule(guided)
		for (i = 0; i < iScanSize; i++) {
			vpAllMS2Scans.at(i)->sumIntensity();
		}
	}
	cout << "Rank " << ProNovoConfig::iRank << ":\tMultiscore preprocess end" << endl;
	//MEMORYSTOP
}

void MS2ScanVector::GetAllRangeFromMass(double dPeptideMass, vector<std::pair<int, int> >& vpPeptideMassRanges)
// all ranges of MS2 scans are stored in  vpPeptideMassWindows
		{
	int i;
	pair<int, int> pairMS2Range;
	pair<int, int> lastPairRange(-100, -100);
	vector<pair<double, double> > vpPeptideMassWindows;
	vpPeptideMassWindows.clear();
	vpPeptideMassRanges.clear();
	ProNovoConfig::getPeptideMassWindows(dPeptideMass, vpPeptideMassWindows);
	for (i = 0; i < (int) vpPeptideMassWindows.size(); i++) {
		pairMS2Range = GetRangeFromMass(vpPeptideMassWindows.at(i).first, vpPeptideMassWindows.at(i).second);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1)) {
			if ((lastPairRange.first < 0) || (lastPairRange.second < 0))
				lastPairRange = pairMS2Range;
			else {
				if (lastPairRange.second > pairMS2Range.first)
					lastPairRange.second = pairMS2Range.second;
				else {
					vpPeptideMassRanges.push_back(lastPairRange);
					lastPairRange = pairMS2Range;
				}
			}
		}
	}
	if ((lastPairRange.first > -1) && (lastPairRange.second > -1))
		vpPeptideMassRanges.push_back(lastPairRange);
}

pair<int, int> MS2ScanVector::GetRangeFromMass(double lb, double ub)
//lb and ub are lower and upper bounds of acceptable parent mass values
		{
	pair<int, int> p;
	int low = 0, high = vpPrecursorMasses.size() - 1, mid;
	double target;
	target = (lb + ub) / 2.0;
	//double lb, ub;  // lower and upper bounds on acceptable parent mass values
	//ub = target + error;
	//lb = target - error;
	while ((high - low) > 1) {
		mid = (high + low) / 2;
		if (vpPrecursorMasses[mid] > target)
			high = mid;
		else
			low = mid;
	}

	// Iterate till we get to the first element > than the lower bound
	int ndx = low;
	//cout<<scan_mass_list_.size()<<endl;
	if (vpPrecursorMasses[ndx] >= lb) {
		while (ndx >= 0 && vpPrecursorMasses[ndx] >= lb)
			ndx--;
		ndx++;
	} else
		while (ndx < (int) vpPrecursorMasses.size() && vpPrecursorMasses[ndx] < lb)
			ndx++;
	if (ndx == (int) vpPrecursorMasses.size() || vpPrecursorMasses[ndx] > ub)
		p = make_pair(-1, -1);
	else {
		low = ndx;
		while (ndx < (int) vpPrecursorMasses.size() && vpPrecursorMasses[ndx] <= ub)
			ndx++;
		high = ndx - 1;
		p = make_pair(low, high);
	}

	if (p.first == -1)
		return p;
	if (vpPrecursorMasses[p.first] < lb)
		cerr << "ERROR L " << vpPrecursorMasses[p.first] << " " << lb << endl;
	if (vpPrecursorMasses[p.second] > ub)
		cerr << "ERROR U " << vpPrecursorMasses[p.second] << " " << ub << endl;
	return p;
}

/*
 void MS2ScanVector::assignPeptides2Scans(Peptide * currentPeptide)
 {
 int i;
 pair<int, int> pairMS2Range;
 //   cout<<currentPeptide->getPeptideMass()<<endl;
 pairMS2Range = GetRangeFromMass(currentPeptide->getPeptideMass(),
 ProNovoConfig::getMassAccuracyParentIon());
 if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1))
 for (i= pairMS2Range.first; i<= pairMS2Range.second; i++)
 vpAllMS2Scans.at(i)->vpPeptides.push_back(currentPeptide);
 }*/

bool MS2ScanVector::assignPeptides2Scans(Peptide * currentPeptide) {
	int i, j;
	bool bGetAssigned = false;
	vector<pair<int, int> > vpPeptideMassRanges;
	pair<int, int> pairMS2Range;

	GetAllRangeFromMass(currentPeptide->getPeptideMass(), vpPeptideMassRanges);

	for (j = 0; j < (int) vpPeptideMassRanges.size(); j++) {
		pairMS2Range = vpPeptideMassRanges.at(j);
		if ((pairMS2Range.first > -1) && (pairMS2Range.second > -1)) {
			for (i = pairMS2Range.first; i <= pairMS2Range.second; i++) {
				vpAllMS2Scans.at(i)->vpPeptides.push_back(currentPeptide);
				vpAllMS2Scans.at(i)->iNumPeptideAssigned++;
			}
			bGetAssigned = true;
		}
	}
	return bGetAssigned;
}

void MS2ScanVector::processPeptideArray(vector<Peptide*>& vpPeptideArray) {
	size_t i, iPeptideArraySize, iScanSize;
	iPeptideArraySize = vpPeptideArray.size();
	iScanSize = vpAllMS2Scans.size();

#pragma omp parallel for shared(vpPeptideArray) private(i) schedule(guided)
	for (i = 0; i < iPeptideArraySize; i++) {
		vpPeptideArray.at(i)->preprocessing(vpAllMS2Scans.at(0)->isMS2HighRes, mapResidueMass);
	}

	if (ProNovoConfig::bMvhEnable) {
#pragma omp parallel
		{
			vector<double> * _pdAAforward = new vector<double>();
			vector<double> * _pdAAreverse = new vector<double>();
			vector<double> * sequenceIonMasses = new vector<double>();
#pragma omp for schedule(guided)
			for (size_t ii = 0; ii < iScanSize; ii++) {
				if (vpAllMS2Scans.at(ii)->bSkip) {
					vpAllMS2Scans.at(ii)->vpPeptides.clear();
					continue;
				}
				size_t j = 0, iPeptideArrayPerScan = vpAllMS2Scans.at(ii)->vpPeptides.size();
				for (j = 0; j < iPeptideArrayPerScan; j++) {
					double dMvh = 0;
					if (!vpAllMS2Scans.at(ii)->mergePeptide(vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides,
							vpAllMS2Scans.at(ii)->vpPeptides.at(j)->getPeptideForScoring(),
							vpAllMS2Scans.at(ii)->vpPeptides.at(j)->getProteinName())) {
						if (MVH::ScoreSequenceVsSpectrum(vpAllMS2Scans.at(ii)->vpPeptides.at(j), vpAllMS2Scans.at(ii),
								sequenceIonMasses, _pdAAforward, _pdAAreverse, dMvh)) {
							vpAllMS2Scans.at(ii)->saveScore(dMvh, vpAllMS2Scans.at(ii)->vpPeptides.at(j),
									vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides,
									vpAllMS2Scans.at(ii)->vdWeightSumAllScores, "MVH");
						}
					}
				}
				vpAllMS2Scans.at(ii)->vpPeptides.clear();
			}
			_pdAAforward->clear();
			delete _pdAAforward;
			_pdAAforward = NULL;
			_pdAAreverse->clear();
			delete _pdAAreverse;
			_pdAAreverse = NULL;
			sequenceIonMasses->clear();
			delete sequenceIonMasses;
			sequenceIonMasses = NULL;
		}
	} // if (ProNovoConfig::bMvhEnable) {

	if (ProNovoConfig::bWeightDotSumEnable) {
#pragma omp parallel for schedule(guided)
		for (i = 0; i < iScanSize; i++) {
			vpAllMS2Scans[i]->scorePeptides();
		}
	}

	// free memory of all peptide objects
	for (i = 0; i < vpPeptideArray.size(); i++) {
		delete vpPeptideArray.at(i);
	}

	// empty peptide array
	vpPeptideArray.clear();
}

void MS2ScanVector::searchDatabaseSnp() {
	cout << "Rank " << ProNovoConfig::iRank << ":\tSearch database begin" << endl;
	ProteinDbParser * db = new ProteinDbParser(bScreenOutput);
	size_t iNumPeptides;
	size_t iNumPtmPeptides;
	Peptide * currentPeptide = NULL;
	vector<Peptide *> vpPeptideArray;
	while (db->getNextProtein()) {
		db->proteinDigest();
		iNumPeptides = db->getNumMutatedPeptide();
		for (size_t i = 0; i < iNumPeptides; i++) {
			db->peptideModify(i);
			iNumPtmPeptides = db->getNumPtmPeptide();
			for (size_t j = 0; j < iNumPtmPeptides; j++) {
				currentPeptide = new Peptide();
				db->setPeptide(j, currentPeptide);
				if (assignPeptides2Scans(currentPeptide)) {
					// cout << currentPeptide->sPeptide << ":\t" << currentPeptide->getPeptideMass() <<endl;
					// save the new peptide to the array
					vpPeptideArray.push_back(currentPeptide);
					if (currentPeptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass) {
						ProNovoConfig::dMaxPeptideMass = currentPeptide->getPeptideMass();
					}
				} else {
					delete currentPeptide;
					currentPeptide = NULL;
				}
				// when the vpPeptideArray is full
				if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE) {
					//cout << "Protein #" << myProteinDatabase.iProteinId << "." << endl;
					processPeptideArray(vpPeptideArray);
				}
			}
		}
	}
	if (!vpPeptideArray.empty()) {
		//cout << "Protein #" << myProteinDatabase.iProteinId << "." << endl;
		processPeptideArray(vpPeptideArray);
	}
	/*// the last peptide object is an empty object and need to be deleted
	 if (currentPeptide != NULL) {
	 delete currentPeptide;
	 }*/
	delete db;
	cout << "Rank " << ProNovoConfig::iRank << ":\tSearch database end" << endl;
	cout << "Rank " << ProNovoConfig::iRank << ":\tdestroy mvh table begin" << endl;
	if (ProNovoConfig::bMvhEnable) {
		MVH::destroyLnTable();
	}
	cout << "Rank " << ProNovoConfig::iRank << ":\tdestroy mvh table end" << endl;
}

void MS2ScanVector::searchDatabase() {
	//MEMORYSTART
	ProteinDatabase myProteinDatabase(bScreenOutput);
	vector<Peptide *> vpPeptideArray;
	Peptide * currentPeptide = NULL;
	myProteinDatabase.loadDatabase();

	cout << "Rank " << ProNovoConfig::iRank << ":\tSearch database begin" << endl;
	if (myProteinDatabase.getFirstProtein()) {
		currentPeptide = new Peptide();
		// get one peptide from the database at a time, until there is no more peptide
		while (myProteinDatabase.getNextPeptide(currentPeptide)) {
			/*if(currentPeptide->getProteinName()=="Rev_UNLACT1_7763G0002"){
			 cout << "check" << endl;
			 }*/
			// cout << currentPeptide->sPeptide << endl;
			/*if (currentPeptide->sPeptide == "[ALLS<LAYCK]") {
				cout << "check" << endl;
			}*/
			// save the new peptide to the array
			if (assignPeptides2Scans(currentPeptide)) {
				// cout << currentPeptide->sPeptide << ":\t" << currentPeptide->getPeptideMass() <<endl;
				// save the new peptide to the array
				vpPeptideArray.push_back(currentPeptide);
				if (currentPeptide->getPeptideMass() > ProNovoConfig::dMaxPeptideMass) {
					ProNovoConfig::dMaxPeptideMass = currentPeptide->getPeptideMass();
				}
			} else {
				delete currentPeptide;
				currentPeptide = NULL;
			}
			// create a new peptide for the next iteration
			currentPeptide = new Peptide;
			// when the vpPeptideArray is full
			if (vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE) {
				//cout << "Protein #" << myProteinDatabase.iProteinId << "." << endl;
				processPeptideArray(vpPeptideArray);
			}
		}

		// there are still unprocessed peptides in the vpPeptideArray
		// need to process them in the same manner
		// the following code is the same as inside if(vpPeptideArray.size() >= PEPTIDE_ARRAY_SIZE )
		// cout<<vpPeptideArray.size()<<endl;
		if (!vpPeptideArray.empty()) {
			//cout << "Protein #" << myProteinDatabase.iProteinId << "." << endl;
			processPeptideArray(vpPeptideArray);
		}
	}
	// the last peptide object is an empty object and need to be deleted
	if (currentPeptide != NULL) {
		delete currentPeptide;
	}
	cout << "Rank " << ProNovoConfig::iRank << ":\tSearch database end" << endl;
	cout << "Rank " << ProNovoConfig::iRank << ":\tdestroy mvh table begin" << endl;
	if (ProNovoConfig::bMvhEnable) {
		MVH::destroyLnTable();
	}
	cout << "Rank " << ProNovoConfig::iRank << ":\tdestroy mvh table end" << endl;
	//MEMORYSTOP
}

void MS2ScanVector::startProcessing() {

	// Preprocessing all MS2 scans by multi-threading
	preProcessAllMS2();

	// Search all MS2 scans against the database by mult-threading
	// searchDatabaseSnp();

	searchDatabase();

	// Postprocessing all MS2 scans' results by mult-threading
	postProcessAllMS2();

	// write results to a SIP file
	if (ProNovoConfig::bMultiScores) {
		//writeOutputMultiScoresPin();
		//writeOutputMultiScoresSip();
		writeOutputMultiScoresSpectrum2MutiPep();
	} else {
		writeOutput();
	}
}

void MS2ScanVector::postProcessAllMS2() {
	//MEMORYSTART
	size_t iScanSize = vpAllMS2Scans.size();
	if (ProNovoConfig::bMultiScores) {
		PeptideUnit::iNumScores = 1;

		cout << "Rank " << ProNovoConfig::iRank << ":\tXcorr begin" << endl;
		if (!ProNovoConfig::bXcorrEnable) {
			//if (false) {
			//MH: Must be equal to largest possible array
			int iArraySizePreprocess = (int) ((ProNovoConfig::dMaxMS2ScanMass + 3 + 2.0)
					* ProNovoConfig::dHighResInverseBinWidth);
			//MH: Must be equal to largest possible array
			int iArraySizeScore =
					(int) ((ProNovoConfig::dMaxPeptideMass + 100) * ProNovoConfig::dHighResInverseBinWidth);
			{
				CometSearch::iArraySizePreprocess = iArraySizePreprocess;
				CometSearch::iArraySizeScore = iArraySizeScore;
				CometSearch::iDimesion2 = 9;
				CometSearch::iMAX_PEPTIDE_LEN = MAX_PEPTIDE_LEN;
				CometSearch::iMaxPercusorCharge = ProNovoConfig::iMaxPercusorCharge + 1;
			}
#pragma omp parallel
			{
				double * pdTmpRawData;
				double * pdTmpFastXcorrData;
				double * pdTmpCorrelationData;
				double * pdTmpSmoothedSpectrum;
				double * pdTmpPeakExtracted;
				bool * pbDuplFragment;
				double * _pdAAforward;
				double * _pdAAreverse;
				unsigned int *** _uiBinnedIonMasses;
				try {
					pdTmpRawData = new double[iArraySizePreprocess]();
					pdTmpFastXcorrData = new double[iArraySizePreprocess]();
					pdTmpCorrelationData = new double[iArraySizePreprocess]();
					pdTmpSmoothedSpectrum = new double[iArraySizePreprocess]();
					pdTmpPeakExtracted = new double[iArraySizePreprocess]();
					pbDuplFragment = new bool[iArraySizeScore]();
					_pdAAforward = new double[MAX_PEPTIDE_LEN]();
					_pdAAreverse = new double[MAX_PEPTIDE_LEN]();
					_uiBinnedIonMasses = new unsigned int**[ProNovoConfig::iMaxPercusorCharge + 1]();
					for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
						_uiBinnedIonMasses[ii] = new unsigned int*[9]();
						for (int j = 0; j < 9; j++) {
							_uiBinnedIonMasses[ii][j] = new unsigned int[MAX_PEPTIDE_LEN]();
						}
					}
				} catch (std::bad_alloc& ba) {
					char szErrorMsg[256];
					sprintf(szErrorMsg, " Error - new(pd[%d]). bad_alloc: %s.\n", iArraySizePreprocess, ba.what());
					sprintf(szErrorMsg + strlen(szErrorMsg),
							"Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
					sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
					string strErrorMsg(szErrorMsg);
					logerr(szErrorMsg);
					exit(1);
				}
#pragma omp for schedule(guided)
				for (size_t ii = 0; ii < iScanSize; ii++) {
					struct Query * pQuery = new Query();
					if (!CometSearch::Preprocess(pQuery, vpAllMS2Scans.at(ii), pdTmpRawData, pdTmpFastXcorrData,
							pdTmpCorrelationData, pdTmpSmoothedSpectrum, pdTmpPeakExtracted)) {
						cout << "Error" << endl;
						exit(1);
					} else {
						vpAllMS2Scans.at(ii)->pQuery = pQuery;
						for (size_t j = 0; j < vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.size(); j++) {
							double dXcorr = 0;
							CometSearch::ScorePeptides(
									&(vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
									pbDuplFragment, _pdAAforward, _pdAAreverse, vpAllMS2Scans.at(ii),
									_uiBinnedIonMasses, dXcorr, iArraySizeScore);
							vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.at(j)->addScores(dXcorr);
						}
						delete pQuery;
						vpAllMS2Scans.at(ii)->pQuery = NULL;
					}
				}
				delete[] pdTmpRawData;
				delete[] pdTmpFastXcorrData;
				delete[] pdTmpCorrelationData;
				delete[] pdTmpSmoothedSpectrum;
				delete[] pdTmpPeakExtracted;
				delete[] pbDuplFragment;
				delete[] _pdAAforward;
				delete[] _pdAAreverse;
				for (int ii = 0; ii < ProNovoConfig::iMaxPercusorCharge + 1; ii++) {
					for (int j = 0; j < 9; j++) {
						delete[] _uiBinnedIonMasses[ii][j];
					}
					delete[] _uiBinnedIonMasses[ii];
				}
				delete[] _uiBinnedIonMasses;
				pdTmpRawData = NULL;
				pdTmpFastXcorrData = NULL;
				pdTmpCorrelationData = NULL;
				pdTmpSmoothedSpectrum = NULL;
				pdTmpPeakExtracted = NULL;
				pbDuplFragment = NULL;
				_pdAAforward = NULL;
				_pdAAreverse = NULL;
				_uiBinnedIonMasses = NULL;
			}
			++PeptideUnit::iNumScores;
		}		//end comet score
		cout << "Rank " << ProNovoConfig::iRank << ":\tXcorr end" << endl;

		cout << "Rank " << ProNovoConfig::iRank << ":\tWDP begin" << endl;
		if (!ProNovoConfig::bWeightDotSumEnable) {
			//if (false) {
			bool bHighRes = vpAllMS2Scans.at(0)->isMS2HighRes;
#pragma omp parallel
			{
				vector<vector<double> > * pvvdYionMass = new vector<vector<double> >();
				vector<vector<double> > * pvvdYionProb = new vector<vector<double> >();
				vector<vector<double> > * pvvdBionMass = new vector<vector<double> >();
				vector<vector<double> > * pvvdBionProb = new vector<vector<double> >();
				vector<double> * pvdYionMass = new vector<double>();
				vector<double> * pvdBionMass = new vector<double>();
				//preprocess peptide
#pragma omp for schedule(guided)
				for (size_t ii = 0; ii < iScanSize; ii++) {
					for (size_t j = 0; j < vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.size(); j++) {
						Peptide::preprocessing(vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring,
								bHighRes, mapResidueMass, pvvdYionMass, pvvdYionProb, pvvdBionMass, pvvdBionProb,
								pvdYionMass, pvdBionMass);
						double dWeightSum = 0;
						if (bHighRes) {
							dWeightSum = vpAllMS2Scans.at(ii)->scoreWeightSumHighMS2(
									&(vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
									pvvdYionMass, pvvdYionProb, pvvdBionMass, pvvdBionProb);
						} else {
							dWeightSum = vpAllMS2Scans.at(ii)->scoreWeightSum(
									&(vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.at(j)->sPeptideForScoring),
									pvdYionMass, pvdBionMass);
						}
						vpAllMS2Scans.at(ii)->vpWeightSumTopPeptides.at(j)->addScores(dWeightSum);
					}
				}

				delete pvvdYionMass;
				delete pvvdYionProb;
				delete pvvdBionMass;
				delete pvvdBionProb;
				delete pvdYionMass;
				delete pvdBionMass;
				pvvdYionMass = NULL;
				pvvdYionProb = NULL;
				pvvdBionMass = NULL;
				pvvdBionProb = NULL;
				pvdYionMass = NULL;
				pvdBionMass = NULL;
			}
			++PeptideUnit::iNumScores;
		}		//end sipros score

#pragma omp for schedule(guided)
		for (size_t ii = 0; ii < iScanSize; ii++) {
			vpAllMS2Scans.at(ii)->scoreFeatureCalculation();
		}
		cout << "Rank " << ProNovoConfig::iRank << ":\tWDP end" << endl;
	}		// end multiple scores
	//MEMORYSTOP
}

bool MS2ScanVector::mylessScanId(MS2Scan* pMS2Scan1, MS2Scan* pMS2Scan2) {
	return (pMS2Scan1->iScanId < pMS2Scan2->iScanId);
}

string MS2ScanVector::ParsePath(string sPath)
// return the filename without path
		{
	string sTailFileName;
	size_t iPosition;
	iPosition = sPath.rfind(ProNovoConfig::getSeparator());
	if (iPosition == string::npos)
		sTailFileName = sPath;
	else
		sTailFileName = sPath.substr(iPosition + 1);
	return sTailFileName;
}

void MS2ScanVector::writeOutput() {
	int i, j;
	ofstream outputFile;
	ifstream configStream;
	string sline, sTailFT2FileName;

	double dcurrentMeanWeightSum, dcurrentDeviationWeightSum;

	sTailFT2FileName = ParsePath(sFT2Filename);

	sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), mylessScanId);
	outputFile.open(sOutputFile.c_str());

	configStream.open(sConfigFile.c_str());
	while (!configStream.eof()) {
		sline.clear();
		getline(configStream, sline);
		outputFile << "#\t" << sline << endl;
	}

	outputFile << "Filename\tScanNumber\tParentCharge\tMeasuredParentMass\t";
	outputFile << "CalculatedParentMass\tScanType\tSearchName\tScoringFunction\tRank\tScore\t";
	outputFile << "IdentifiedPeptide\tOriginalPeptide\tProteinNames" << endl;
	for (i = 0; i < (int) vpAllMS2Scans.size(); i++)
		if (!vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.empty()) {
			if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(0)->sScoringFunction == "WeightSum") {
				calculateMeanAndDeviation(vpAllMS2Scans.at(i)->inumberofWeightSumScore,
						vpAllMS2Scans.at(i)->dsumofWeightScore, vpAllMS2Scans.at(i)->dsumofSquareWeightSumScore,
						dcurrentMeanWeightSum, dcurrentDeviationWeightSum);
			}
			for (j = 0; j < (int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				/*		cout<<sFT2Filename<<"\t"<<vpAllMS2Scans.at(i)->iScanId<<"\t";
				 cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getPeptideScore()<<"\t";
				 cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getPeptideSeq()<<"\t";
				 cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getOriginalPeptideSeq()<<"\t";
				 cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getProteinName()<<"\t";
				 cout<<vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->getBeginPosProtein()<<"\t";
				 cout<<vpAllMS2Scans.at(i)->iParentChargeState<<"\t";
				 cout<<i<<"\t";
				 cout<<vpAllMS2Scans.at(i)->vdWeightSumAllScores.size()<<"\t0\thit"<<endl;
				 */
				outputFile << sTailFT2FileName << "\t";
				outputFile << vpAllMS2Scans.at(i)->iScanId << "\t";
				outputFile << vpAllMS2Scans.at(i)->iParentChargeState << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(5) << vpAllMS2Scans.at(i)->dParentNeutralMass
						<< "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(5)
						<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dCalculatedParentMass << "\t";

				outputFile << vpAllMS2Scans.at(i)->getScanType() << "\t";
				outputFile << ProNovoConfig::getSearchName() << "\t";

				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sScoringFunction << "\t";
				outputFile << (j + 1) << "\t";

				outputFile << setiosflags(ios::fixed) << setprecision(4)
						<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dScore << "\t";
				//outputFile<<((vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dScore-dcurrentMeanWeightSum)/dcurrentDeviationWeightSum)<<"\t";

				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix;
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix;
				outputFile << "\t";

				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix;
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sOriginalPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix;
				outputFile << "\t";

				outputFile << "{" << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sProteinNames << "}" << endl;
			}
		}
	configStream.clear();
	configStream.close();
	outputFile.close();
}

void MS2ScanVector::writeOutputMultiScoresSip() {
	int i, j;
	ofstream outputFile;
	ifstream configStream;
	string sline, sTailFT2FileName;

	double dcurrentMeanWeightSum, dcurrentDeviationWeightSum;

	sTailFT2FileName = ParsePath(sFT2Filename);

	sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), mylessScanId);
	outputFile.open(sOutputFile.c_str());

	configStream.open(sConfigFile.c_str());
	while (!configStream.eof()) {
		sline.clear();
		getline(configStream, sline);
		outputFile << "#\t" << sline << endl;
	}

	outputFile << "Filename\tScanNumber\tParentCharge\tMeasuredParentMass\t";
	outputFile << "CalculatedParentMass\tScanType\tSearchName\tRank\tMVH\tXcorr\tWeightDotSum\t";
	outputFile << "IdentifiedPeptide\tOriginalPeptide\tProteinNames" << endl;
	for (i = 0; i < (int) vpAllMS2Scans.size(); i++)
		if (!vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.empty()) {
			if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(0)->sScoringFunction == "WeightSum") {
				calculateMeanAndDeviation(vpAllMS2Scans.at(i)->inumberofWeightSumScore,
						vpAllMS2Scans.at(i)->dsumofWeightScore, vpAllMS2Scans.at(i)->dsumofSquareWeightSumScore,
						dcurrentMeanWeightSum, dcurrentDeviationWeightSum);
			}
			for (j = 0; j < (int) vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				outputFile << sTailFT2FileName << "\t";
				outputFile << vpAllMS2Scans.at(i)->iScanId << "\t";
				outputFile << vpAllMS2Scans.at(i)->iParentChargeState << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(5) << vpAllMS2Scans.at(i)->dParentNeutralMass
						<< "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(5)
						<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dCalculatedParentMass << "\t";

				outputFile << vpAllMS2Scans.at(i)->getScanType() << "\t";
				outputFile << ProNovoConfig::getSearchName() << "\t";

				//outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sScoringFunction << "\t";
				outputFile << (j + 1) << "\t";
				for (size_t k = 0; k < vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.size(); k++) {
					outputFile << setiosflags(ios::fixed) << setprecision(4)
							<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.at(k) << "\t";
				}

				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix;
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix;
				outputFile << "\t";

				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix;
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sOriginalPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix != '-')
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix;
				outputFile << "\t";

				outputFile << "{" << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sProteinNames << "}" << endl;
			}
		}
	configStream.clear();
	configStream.close();
	outputFile.close();
}

void MS2ScanVector::writeOutputMultiScoresSpectrum2MutiPep() {
	//MEMORYSTART
	size_t i, j;
	size_t k = 0;
	ofstream outputFile;
	ifstream configStream;
	string sline, sTailFT2FileName;

	sTailFT2FileName = ParsePath(sFT2Filename);

	sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), mylessScanId);

	string sOutputFileTxt = sOutputFile.substr(0, sOutputFile.length() - 4);
	sOutputFileTxt += "_Spe2Pep.txt";
	outputFile.open(sOutputFileTxt.c_str());

	configStream.open(sConfigFile.c_str());
	while (!configStream.eof()) {
		sline.clear();
		getline(configStream, sline);
		outputFile << "#\t" << sline << endl;
	}

	//Spectrum level head
	outputFile << "+\tFilename\tScanNumber\tParentCharge\tMeasuredParentMass"
			<< "\tScanType\tSearchName\tTotalIntensity\tMaxIntensity\tRetentionTime" << endl;
	//PSM level head
	outputFile << "*\tIdentifiedPeptide\tOriginalPeptide\tCalculatedParentMass " << "\tMVH\tXcorr\tWDP\tProteinNames"
			<< endl;
	cout << "Rank " << ProNovoConfig::iRank << ":\tWrite begin" << endl;
	for (i = 0; i < vpAllMS2Scans.size(); i++) {
		if (!vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.empty()) {
			//Spectrum level information
			outputFile << "+";
			outputFile << "\t" << sTailFT2FileName;
			outputFile << "\t" << vpAllMS2Scans.at(i)->iScanId;
			outputFile << "\t" << vpAllMS2Scans.at(i)->iParentChargeState;
			outputFile << "\t" << setiosflags(ios::fixed) << setprecision(5) << vpAllMS2Scans.at(i)->dParentNeutralMass;
			outputFile << "\t" << vpAllMS2Scans.at(i)->getScanType();
			outputFile << "\t" << ProNovoConfig::getSearchName();
			//sInt sum of intensity
			outputFile << "\t" << vpAllMS2Scans.at(i)->dSumIntensity;
			//maxInt
			outputFile << "\t" << vpAllMS2Scans.at(i)->dMaxIntensity;
			// Retention Time
			outputFile << "\t" << vpAllMS2Scans.at(i)->getRTime();
			outputFile << endl;
			//PSM level head
			for (j = 0; j < vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				if (!vpAllMS2Scans.at(i)->isAnyScoreInTopN(j, ProNovoConfig::Log_TOP_N_Output)) {
					continue;
				}
				outputFile << "*\t";
				//IdentifiedPeptide
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix;
				}
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix;
				}
				outputFile << "\t";
				//OriginalPeptide
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalPrefix;
				}
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sOriginalPeptide;
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix != '-') {
					outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cOriginalSuffix;
				}

				//CalculatedParentMass
				outputFile << "\t";
				outputFile << setiosflags(ios::fixed) << setprecision(4)
						<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dPepMass << "\t";
				//MVH, Xcorr, WDP
				for (k = 0; k < vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.size(); k++) {
					outputFile << setiosflags(ios::fixed) << setprecision(4)
							<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.at(k) << "\t";
				}
				//ProteinNames
				outputFile << "{" << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sProteinNames << "}" << endl;
			}
		}
	}
	cout << "Rank " << ProNovoConfig::iRank << ":\tWrite end" << endl;
	configStream.clear();
	configStream.close();
	outputFile.close();
	//MEMORYSTOP
}

void MS2ScanVector::writeOutputMultiScoresPin() {
	size_t i, j, k;
	ofstream outputFile;
	ifstream configStream;
	string sline, sTailFT2FileName;
	double temp = 0;
	sTailFT2FileName = ParsePath(sFT2Filename);

	sort(vpAllMS2Scans.begin(), vpAllMS2Scans.end(), mylessScanId);
	string sOutputFilePin = sOutputFile.substr(0, sOutputFile.length() - 4);
	sOutputFilePin += ".pin";
	outputFile.open(sOutputFilePin.c_str());

	outputFile << "PSMId\tLabel\tScanNr\t";
	for (i = 0; i < PeptideUnit::iNumScores; i++) {
		outputFile << "Score" << (i + 1) << "\t";
		outputFile << "Fraction" << (i + 1) << "\t";
		outputFile << "Rank" << (i + 1) << "\t";
	}
	outputFile << "Mass\t";
	outputFile << "dM\t";
	outputFile << "absdM\t";
	outputFile << "sInt\t";
	outputFile << "maxInt\t";
	outputFile << "pepLen\t";
	outputFile << "charge1\t";
	outputFile << "charge2\t";
	outputFile << "charge3\t";
	outputFile << "enzInt\t";
	outputFile << "lnNumPep\t";
	outputFile << "Peptide\t";
	iMaxNumProteins = 0;
	size_t iTemp = 0;
	for (i = 0; i < vpAllMS2Scans.size(); i++) {
		iTemp = vpAllMS2Scans.at(i)->getMaxNumProteinPsm();
		if (iTemp > iMaxNumProteins) {
			iMaxNumProteins = iTemp;
		}
	}
	for (i = 0; i < iMaxNumProteins; i++) {
		outputFile << "Proteins" << (i + 1) << "\t";
	}
	outputFile << endl;

	for (i = 0; i < vpAllMS2Scans.size(); i++) {
		if (!vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.empty()) {
			for (j = 0; j < vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.size(); j++) {
				outputFile << sTailFT2FileName << "_";
				outputFile << vpAllMS2Scans.at(i)->iScanId << "_" << j << "\t";
				if (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->isDecoy(ProNovoConfig::sDecoyPrefix)) {
					outputFile << "-1\t";
				} else {
					outputFile << "1\t";
				}
				outputFile << vpAllMS2Scans.at(i)->iScanId << "\t";
				//scores and related features
				for (k = 0; k < vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.size(); k++) {
					outputFile << setiosflags(ios::fixed) << setprecision(4)
							<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdScores.at(k) << "\t";
					outputFile << setiosflags(ios::fixed) << setprecision(4)
							<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdFraction.at(k) << "\t";
					outputFile << setiosflags(ios::fixed) << setprecision(4)
							<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->vdRank.at(k) << "\t";
				}
				//Mass
				outputFile << vpAllMS2Scans.at(i)->dParentMass << "\t";
				//dM
				temp = (vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->dPepMass - vpAllMS2Scans.at(i)->dParentMass);
				outputFile << temp << "\t";
				//absdM
				outputFile << fabs(temp) << "\t";
				//sInt sum of intensity
				outputFile << vpAllMS2Scans.at(i)->dSumIntensity << "\t";
				//maxInt
				outputFile << vpAllMS2Scans.at(i)->dMaxIntensity << "\t";
				//pepLen
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->iPepLen << "\t";
				//charge
				if (vpAllMS2Scans.at(i)->iParentChargeState == 1) {
					outputFile << "1" << "\t";
				} else {
					outputFile << "0" << "\t";
				}
				if (vpAllMS2Scans.at(i)->iParentChargeState == 2) {
					outputFile << "1" << "\t";
				} else {
					outputFile << "0" << "\t";
				}
				if (vpAllMS2Scans.at(i)->iParentChargeState >= 3) {
					outputFile << "1" << "\t";
				} else {
					outputFile << "0" << "\t";
				}
				//enzInt
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->iEnzInt << "\t";
				//lnNumPep
				outputFile << (log(vpAllMS2Scans.at(i)->iNumPeptideAssigned)) << "\t";
				//Peptide
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifyPrefix;
				outputFile << ".";
				outputFile
						<< vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide.substr(1,
								vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sIdentifiedPeptide.length() - 2);
				outputFile << ".";
				outputFile << vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->cIdentifySuffix;
				outputFile << "\t";
				//Proteins
				istringstream ss(vpAllMS2Scans.at(i)->vpWeightSumTopPeptides.at(j)->sProteinNames);
				string token;
				while (std::getline(ss, token, ',')) {
					outputFile << token << '\t';
				}
				outputFile << endl;
			}
		}
	}
	configStream.clear();
	configStream.close();
	outputFile.close();
}

void MS2ScanVector::calculateMeanAndDeviation(int inumberScore, double dScoreSum, double dScoreSquareSum, double& dMean,
		double& dDeviation) {
	dMean = dScoreSum / inumberScore;
	//dDeviation is sample standard deviation

	if (inumberScore > 1) {
		dDeviation = sqrt((dScoreSquareSum - dMean * dMean) * inumberScore / (inumberScore - 1));
	}

	if ((inumberScore == 1) || (dDeviation < ZERO)) {
		dDeviation = ZERO;
	}
}

