/*
 * CometSearch.h
 *
 *  Created on: May 13, 2016
 *      Author: xgo
 */

#ifndef SCORES_COMETSEARCH_H_
#define SCORES_COMETSEARCH_H_

class MS2Scan;

#include "../proNovoConfig.h"
#include "../ms2scan.h"
#include <ctype.h>

// Redefined how the bin offset is interpreted and applied.  The valid range for the offset is
// now between 0.0 and 1.0 and scales to the binWidth.
//#define BIN(dMass) (int)(dMass*ProNovoConfig::dInverseBinWidth + ProNovoConfig::dOneMinusBinOffset)
#define BINX(dMass, dInverseBinWidth, dOneMinusBinOffset) (int)(dMass*dInverseBinWidth + dOneMinusBinOffset)
#define PROTON_MASS 1.00727646688
#define logout(szString) fputs(szString, stdout)
#define logerr(szString) fputs(szString, stderr)
#define XCORR_CUTOFF                1E-8   // some near-zero cutoff
struct msdata                    // used in the preprocessing
{
	double dIon;
	double dIntensity;
};

// PreprocessStruct stores information used in preprocessing
// each spectrum.  Information not kept around otherwise
struct PreprocessStruct {
	int iHighestIon;
	double dHighestIntensity;
	struct msdata pTmpSpData[NUM_SP_IONS];
};

struct SpectrumInfoInternal {
	int iArraySize;     // m/z versus intensity array
	int iHighestIon;
	int iChargeState;
	int iMaxFragCharge;
	double dTotalIntensity;
};

struct PepMassInfo {
	double dCalcPepMass;
	double dExpPepMass;
	double dPeptideMassTolerance;
	double dPeptideMassToleranceMinus;
	double dPeptideMassTolerancePlus;
};

// Query stores information for peptide scoring and results
// This struct is allocated for each spectrum/charge combination
struct Query {

	// Sparse matrix representation of data
	int iSpScoreData;    //size of sparse matrix
	int iFastXcorrData;  //MH: I believe these are all the same size now.
	int iFastXcorrDataNL;
	float **ppfSparseSpScoreData;
	float **ppfSparseFastXcorrData;
	float **ppfSparseFastXcorrDataNL;

	// Standard array representation of data
	float *pfSpScoreData;
	float *pfFastXcorrData;
	float *pfFastXcorrDataNL;  // pfFastXcorrData with NH3, H2O contributions

	SpectrumInfoInternal _spectrumInfoInternal;

	Query() {

		iSpScoreData = 0;
		iFastXcorrData = 0;
		iFastXcorrDataNL = 0;

		ppfSparseSpScoreData = NULL;
		ppfSparseFastXcorrData = NULL;
		ppfSparseFastXcorrDataNL = NULL;          // pfFastXcorrData with NH3, H2O contributions

		pfSpScoreData = NULL;
		pfFastXcorrData = NULL;
		pfFastXcorrDataNL = NULL;              // pfFastXcorrData with NH3, H2O contributions

		_spectrumInfoInternal.dTotalIntensity = 0.0;
		_spectrumInfoInternal.iArraySize = 0;
		_spectrumInfoInternal.iHighestIon = 0;
		_spectrumInfoInternal.dTotalIntensity = 0.0;
	}

	~Query() {
		int i;
		for (i = 0; i < iSpScoreData; i++) {
			if (ppfSparseSpScoreData[i] != NULL)
				delete[] ppfSparseSpScoreData[i];
		}
		delete[] ppfSparseSpScoreData;
		ppfSparseSpScoreData = NULL;

		for (i = 0; i < iFastXcorrData; i++) {
			if (ppfSparseFastXcorrData[i] != NULL)
				delete[] ppfSparseFastXcorrData[i];
		}
		delete[] ppfSparseFastXcorrData;
		ppfSparseFastXcorrData = NULL;

		if (ProNovoConfig::ionInformation.bUseNeutralLoss
				&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A]
						|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
						|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
			for (i = 0; i < iFastXcorrDataNL; i++) {
				if (ppfSparseFastXcorrDataNL[i] != NULL)
					delete[] ppfSparseFastXcorrDataNL[i];
			}
			delete[] ppfSparseFastXcorrDataNL;
			ppfSparseFastXcorrDataNL = NULL;
		}
	}
};

struct BinnedIonMasses {
	unsigned int _uiBinnedIonMasses[MAX_FRAGMENT_CHARGE + 1][9][MAX_PEPTIDE_LEN];
};

struct IonSeriesStruct         // defines which fragment ion series are considered
{
	int bPreviousMatch[8];
};

class CometSearch {
public:
	CometSearch();
	~CometSearch();

	static bool Preprocess(struct Query *pScoring, MS2Scan * mstSpectrum, double *pdTmpRawData,
			double *pdTmpFastXcorrData, double *pdTmpCorrelationData, double *pdSmoothedSpectrum,
			double *pdTmpPeakExtracted);
	static bool LoadIons(struct Query *pScoring, double *pdTmpRawData, MS2Scan * mstSpectrum,
			struct PreprocessStruct *pPre);
	static void MakeCorrData(double *pdTmpRawData, double *pdTmpCorrelationData, struct Query *pScoring,
			struct PreprocessStruct *pPre);
	static bool Smooth(double *data, int iArraySize, double *pdSmoothedSpectrum);
	static bool PeakExtract(double *data, int iArraySize, double *pdTmpPeakExtracted);
	static void GetTopIons(double *pdTmpRawData, struct msdata *pTmpSpData, int iArraySize);
	static int QsortByIon(const void *p0, const void *p1);
	static void StairStep(struct msdata *pTmpSpData, double dFragmentBinSize);
	static void print(struct Query *pScoring);

	static bool ScorePeptides(Peptide * currentPeptide, bool *pbDuplFragment, double* _pdAAforward,
			double * _pdAAreverse, MS2Scan * mstSpectrum, unsigned int *** _uiBinnedIonMasses, double & dXcorr,
			int test);
	static double GetFragmentIonMass(int iWhichIonSeries, int i, int ctCharge, double *_pdAAforward,
			double *_pdAAreverse);

	static bool CalculateSP(double & fScoreSp, double* _pdAAforward, double * _pdAAreverse, MS2Scan * mstSpectrum, int iLenPeptide);
	static double FindSpScore(Query *pQuery, int bin);
};

#endif /* SCORES_COMETSEARCH_H_ */
