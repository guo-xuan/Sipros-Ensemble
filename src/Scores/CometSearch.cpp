/*
 * CometSearch.cpp
 *
 *  Created on: May 13, 2016
 *      Author: xgo
 */

#include "CometSearch.h"

CometSearch::CometSearch() {

}

CometSearch::~CometSearch() {

}

bool CometSearch::Preprocess(struct Query *pScoring, MS2Scan * mstSpectrum, double *pdTmpRawData,
		double *pdTmpFastXcorrData, double *pdTmpCorrelationData, double *pdTmpSmoothedSpectrum,
		double *pdTmpPeakExtracted) {
	int i;
	int x;
	int y;
	struct msdata pTmpSpData[NUM_SP_IONS];
	struct PreprocessStruct pPre;
	double dInverseBinWidth = 0, iMinus17 = 0, iMinus18 = 0, dFragmentBinSize = 0;
	//mstSpectrum->isMS2HighRes = false;
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		iMinus17 = ProNovoConfig::precalcMasses.iMinus17HighRes;
		iMinus18 = ProNovoConfig::precalcMasses.iMinus18HighRes;
		dFragmentBinSize = ProNovoConfig::dHighResFragmentBinSize;
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		iMinus17 = ProNovoConfig::precalcMasses.iMinus17LowRes;
		iMinus18 = ProNovoConfig::precalcMasses.iMinus18LowRes;
		dFragmentBinSize = ProNovoConfig::dLowResFragmentBinSize;
	}
	pPre.iHighestIon = 0;
	pPre.dHighestIntensity = 0;

	//MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
	// very wide mass tolerance searches (i.e. 500 Da).
	double dCushion = 3.0;
	pScoring->_spectrumInfoInternal.iArraySize = (int) ((mstSpectrum->dParentMass + dCushion + 2.0) * dInverseBinWidth);
	if (mstSpectrum->iParentChargeState == 1) {
		pScoring->_spectrumInfoInternal.iMaxFragCharge = 1;
	} else {
		pScoring->_spectrumInfoInternal.iMaxFragCharge = mstSpectrum->iParentChargeState - 1;
	}
	if (mstSpectrum->iParentChargeState > ProNovoConfig::iMaxPercusorCharge) {
		ProNovoConfig::iMaxPercusorCharge = mstSpectrum->iParentChargeState;
	}

	// initialize these temporary arrays before re-using
	size_t iTmp = (size_t) ((ProNovoConfig::dMaxMS2ScanMass + dCushion + 2.0) * dInverseBinWidth) * sizeof(double);
	memset(pdTmpRawData, 0, iTmp);
	memset(pdTmpFastXcorrData, 0, iTmp);
	memset(pdTmpCorrelationData, 0, iTmp);
	memset(pdTmpSmoothedSpectrum, 0, iTmp);
	memset(pdTmpPeakExtracted, 0, iTmp);

	// pdTmpRawData is a binned array holding raw data
	if (!LoadIons(pScoring, pdTmpRawData, mstSpectrum, &pPre)) {
		return false;
	}

	try {
		pScoring->pfFastXcorrData = new float[pScoring->_spectrumInfoInternal.iArraySize]();
	} catch (std::bad_alloc& ba) {
		char szErrorMsg[256];
		sprintf(szErrorMsg, " Error - new(pfFastXcorrData[%d]). bad_alloc: %s.\n",
				pScoring->_spectrumInfoInternal.iArraySize, ba.what());
		sprintf(szErrorMsg + strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
		sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
		string strErrorMsg(szErrorMsg);
		logerr(szErrorMsg);
		return false;
	}

	if (ProNovoConfig::ionInformation.bUseNeutralLoss
			&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A]
					|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
					|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
		try {
			pScoring->pfFastXcorrDataNL = new float[pScoring->_spectrumInfoInternal.iArraySize]();
		} catch (std::bad_alloc& ba) {
			char szErrorMsg[256];
			sprintf(szErrorMsg, " Error - new(pfFastXcorrDataNL[%d]). bad_alloc: %s.\n",
					pScoring->_spectrumInfoInternal.iArraySize, ba.what());
			sprintf(szErrorMsg + strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
			sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
			string strErrorMsg(szErrorMsg);
			logerr(szErrorMsg);
			return false;
		}
	}

	// Create data for correlation analysis.
	// pdTmpRawData intensities are normalized to 100; pdTmpCorrelationData is windowed
	MakeCorrData(pdTmpRawData, pdTmpCorrelationData, pScoring, &pPre);

	// Make fast xcorr spectrum.
	double dSum = 0.0;
	int iTmpRange = 2 * ProNovoConfig::iXcorrProcessingOffset + 1;
	double dTmp = 1.0 / (double) (iTmpRange - 1);

	dSum = 0.0;
	for (i = 0; i < ProNovoConfig::iXcorrProcessingOffset; i++)
		dSum += pdTmpCorrelationData[i];
	for (i = ProNovoConfig::iXcorrProcessingOffset;
			i < pScoring->_spectrumInfoInternal.iArraySize + ProNovoConfig::iXcorrProcessingOffset; i++) {
		if (i < pScoring->_spectrumInfoInternal.iArraySize)
			dSum += pdTmpCorrelationData[i];
		if (i >= iTmpRange)
			dSum -= pdTmpCorrelationData[i - iTmpRange];
		pdTmpFastXcorrData[i - ProNovoConfig::iXcorrProcessingOffset] = (dSum
				- pdTmpCorrelationData[i - ProNovoConfig::iXcorrProcessingOffset]) * dTmp;
	}

	pScoring->pfFastXcorrData[0] = 0.0;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		double dTmp = pdTmpCorrelationData[i] - pdTmpFastXcorrData[i];

		pScoring->pfFastXcorrData[i] = (float) dTmp;

		// Add flanking peaks if used if it is high resolution
		if (mstSpectrum->isMS2HighRes) {
			int iTmp;

			iTmp = i - 1;
			pScoring->pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.5);

			iTmp = i + 1;
			if (iTmp < pScoring->_spectrumInfoInternal.iArraySize)
				pScoring->pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.5);
		}

		// If A, B or Y ions and their neutral loss selected, roll in -17/-18 contributions to pfFastXcorrDataNL
		if (ProNovoConfig::ionInformation.bUseNeutralLoss
				&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A]
						|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
						|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
			int iTmp;

			pScoring->pfFastXcorrDataNL[i] = pScoring->pfFastXcorrData[i];

			iTmp = i - iMinus17;
			if (iTmp >= 0) {
				pScoring->pfFastXcorrDataNL[i] +=
						(float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.2);
			}

			iTmp = i - iMinus18;
			if (iTmp >= 0) {
				pScoring->pfFastXcorrDataNL[i] +=
						(float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.2);
			}

		}
	}

	// Using sparse matrix which means we free pScoring->pfFastXcorrData, ->pfFastXcorrDataNL here
	// If A, B or Y ions and their neutral loss selected, roll in -17/-18 contributions to pfFastXcorrDataNL.
	if (ProNovoConfig::ionInformation.bUseNeutralLoss
			&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A]
					|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
					|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
		pScoring->iFastXcorrDataNL = pScoring->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

		try {
			pScoring->ppfSparseFastXcorrDataNL = new float*[pScoring->iFastXcorrDataNL]();
		} catch (std::bad_alloc& ba) {
			char szErrorMsg[256];
			sprintf(szErrorMsg, " Error - new(pScoring->ppfSparseFastXcorrDataNL[%d]). bad_alloc: %s.",
					pScoring->iFastXcorrDataNL, ba.what());
			sprintf(szErrorMsg + strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
			sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
			string strErrorMsg(szErrorMsg);
			logerr(szErrorMsg);
			return false;
		}

		for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
			if (pScoring->pfFastXcorrDataNL[i] > FLOAT_ZERO || pScoring->pfFastXcorrDataNL[i] < -FLOAT_ZERO) {
				x = i / SPARSE_MATRIX_SIZE;
				if (pScoring->ppfSparseFastXcorrDataNL[x] == NULL) {
					try {
						pScoring->ppfSparseFastXcorrDataNL[x] = new float[SPARSE_MATRIX_SIZE]();
					} catch (std::bad_alloc& ba) {
						char szErrorMsg[256];
						sprintf(szErrorMsg,
								" Error - new(pScoring->ppfSparseFastXcorrDataNL[%d][%d]). bad_alloc: %s.\n", x,
								SPARSE_MATRIX_SIZE, ba.what());
						sprintf(szErrorMsg + strlen(szErrorMsg),
								"Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
						sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
						string strErrorMsg(szErrorMsg);
						logerr(szErrorMsg);
						return false;
					}
					for (y = 0; y < SPARSE_MATRIX_SIZE; y++)
						pScoring->ppfSparseFastXcorrDataNL[x][y] = 0;
				}
				y = i - (x * SPARSE_MATRIX_SIZE);
				pScoring->ppfSparseFastXcorrDataNL[x][y] = pScoring->pfFastXcorrDataNL[i];
			}
		}

		delete[] pScoring->pfFastXcorrDataNL;
		pScoring->pfFastXcorrDataNL = NULL;

	}

	pScoring->iFastXcorrData = pScoring->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

	//MH: Fill sparse matrix
	try {
		pScoring->ppfSparseFastXcorrData = new float*[pScoring->iFastXcorrData]();
	} catch (std::bad_alloc& ba) {
		char szErrorMsg[256];
		sprintf(szErrorMsg, " Error - new(pScoring->ppfSparseFastXcorrData[%d]). bad_alloc: %s.\n",
				pScoring->iFastXcorrData, ba.what());
		sprintf(szErrorMsg + strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
		sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
		string strErrorMsg(szErrorMsg);
		logerr(szErrorMsg);
		return false;
	}

	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		if (pScoring->pfFastXcorrData[i] > FLOAT_ZERO || pScoring->pfFastXcorrData[i] < -FLOAT_ZERO) {
			x = i / SPARSE_MATRIX_SIZE;
			if (pScoring->ppfSparseFastXcorrData[x] == NULL) {
				try {
					pScoring->ppfSparseFastXcorrData[x] = new float[SPARSE_MATRIX_SIZE]();
				} catch (std::bad_alloc& ba) {
					char szErrorMsg[256];
					sprintf(szErrorMsg, " Error - new(pScoring->ppfSparseFastXcorrData[%d][%d]). bad_alloc: %s.\n", x,
					SPARSE_MATRIX_SIZE, ba.what());
					sprintf(szErrorMsg + strlen(szErrorMsg),
							"Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
					sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
					string strErrorMsg(szErrorMsg);
					logerr(szErrorMsg);
					return false;
				}
				for (y = 0; y < SPARSE_MATRIX_SIZE; y++)
					pScoring->ppfSparseFastXcorrData[x][y] = 0;
			}
			y = i - (x * SPARSE_MATRIX_SIZE);
			pScoring->ppfSparseFastXcorrData[x][y] = pScoring->pfFastXcorrData[i];
		}
	}

	delete[] pScoring->pfFastXcorrData;
	pScoring->pfFastXcorrData = NULL;

	// Create data for sp scoring.
	// Arbitrary bin size cutoff to do smoothing, peak extraction.
	if (dFragmentBinSize >= 0.10) {
		if (!Smooth(pdTmpRawData, pScoring->_spectrumInfoInternal.iArraySize, pdTmpSmoothedSpectrum)) {
			return false;
		}

		if (!PeakExtract(pdTmpRawData, pScoring->_spectrumInfoInternal.iArraySize, pdTmpPeakExtracted)) {
			return false;
		}
	}

	for (i = 0; i < NUM_SP_IONS; i++) {
		pTmpSpData[i].dIon = 0.0;
		pTmpSpData[i].dIntensity = 0.0;
	}

	GetTopIons(pdTmpRawData, &(pTmpSpData[0]), pScoring->_spectrumInfoInternal.iArraySize);

	qsort(pTmpSpData, NUM_SP_IONS, sizeof(struct msdata), QsortByIon);

	// Modify for Sp data.
	StairStep(pTmpSpData, dFragmentBinSize);

	try {
		pScoring->pfSpScoreData = new float[pScoring->_spectrumInfoInternal.iArraySize]();
	} catch (std::bad_alloc& ba) {
		char szErrorMsg[256];
		sprintf(szErrorMsg, " Error - new(pfSpScoreData[%d]). bad_alloc: %s.\n",
				pScoring->_spectrumInfoInternal.iArraySize, ba.what());
		sprintf(szErrorMsg + strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
		sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
		string strErrorMsg(szErrorMsg);
		logerr(szErrorMsg);
		return false;
	}

	// note that pTmpSpData[].dIon values are already BIN'd
	for (i = 0; i < NUM_SP_IONS; i++)
		pScoring->pfSpScoreData[(int) (pTmpSpData[i].dIon)] = (float) pTmpSpData[i].dIntensity;

	// MH: Fill sparse matrix for SpScore
	pScoring->iSpScoreData = pScoring->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

	try {
		pScoring->ppfSparseSpScoreData = new float*[pScoring->iSpScoreData]();
	} catch (std::bad_alloc& ba) {
		char szErrorMsg[256];
		sprintf(szErrorMsg, " Error - new(pScoring->ppfSparseSpScoreData[%d]). bad_alloc: %s.\n",
				pScoring->iSpScoreData, ba.what());
		sprintf(szErrorMsg + strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
		sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
		string strErrorMsg(szErrorMsg);
		logerr(szErrorMsg);
		return false;
	}

	for (i = 0; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		if (pScoring->pfSpScoreData[i] > FLOAT_ZERO) {
			x = i / SPARSE_MATRIX_SIZE;
			if (pScoring->ppfSparseSpScoreData[x] == NULL) {
				try {
					pScoring->ppfSparseSpScoreData[x] = new float[SPARSE_MATRIX_SIZE]();
				} catch (std::bad_alloc& ba) {
					char szErrorMsg[256];
					sprintf(szErrorMsg, " Error - new(pScoring->ppfSparseSpScoreData[%d][%d]). bad_alloc: %s.\n", x,
					SPARSE_MATRIX_SIZE, ba.what());
					sprintf(szErrorMsg + strlen(szErrorMsg),
							"Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
					sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
					string strErrorMsg(szErrorMsg);
					logerr(szErrorMsg);
					return false;
				}
				for (y = 0; y < SPARSE_MATRIX_SIZE; y++)
					pScoring->ppfSparseSpScoreData[x][y] = 0;
			}
			y = i - (x * SPARSE_MATRIX_SIZE);
			pScoring->ppfSparseSpScoreData[x][y] = pScoring->pfSpScoreData[i];
		}
	}

	delete[] pScoring->pfSpScoreData;
	pScoring->pfSpScoreData = NULL;
	//print(pScoring);
	return true;
}

//  Reads MSMS data file as ASCII mass/intensity pairs.
bool CometSearch::LoadIons(struct Query *pScoring, double *pdTmpRawData, MS2Scan * mstSpectrum,
		struct PreprocessStruct *pPre) {
	double dIon, dIntensity;
	double dInverseBinWidth =
			mstSpectrum->isMS2HighRes ?
					(ProNovoConfig::dHighResInverseBinWidth) : (ProNovoConfig::dLowResInverseBinWidth);
	double dOneMinusBinOffset =
			mstSpectrum->isMS2HighRes ?
					(ProNovoConfig::dHighResOneMinusBinOffset) : (ProNovoConfig::dLowResOneMinusBinOffset);
	size_t i = 0;
	while (true) {
		if (i >= mstSpectrum->vdMZ.size()) {
			break;
		}
		dIon = mstSpectrum->vdMZ.at(i);
		dIntensity = mstSpectrum->vdIntensity.at(i);
		i++;

		pScoring->_spectrumInfoInternal.dTotalIntensity += dIntensity;

		if ((dIntensity >= ProNovoConfig::options.dMinIntensity) && (dIntensity > 0.0)) {
			if (dIon < (mstSpectrum->dParentMass + 50.0)) {
				int iBinIon = BINX(dIon, dInverseBinWidth, dOneMinusBinOffset);

				dIntensity = sqrt(dIntensity);

				if (iBinIon > pPre->iHighestIon)
					pPre->iHighestIon = iBinIon;

				if ((iBinIon < pScoring->_spectrumInfoInternal.iArraySize) && (dIntensity > pdTmpRawData[iBinIon])) {
					if (ProNovoConfig::options.iRemovePrecursor == 1) {
						double dMZ = (mstSpectrum->dParentMass
								+ (pScoring->_spectrumInfoInternal.iChargeState - 1) * PROTON_MASS)
								/ (double) (pScoring->_spectrumInfoInternal.iChargeState);
						if (fabs(dIon - dMZ) > ProNovoConfig::options.dRemovePrecursorTol) {
							if (dIntensity > pdTmpRawData[iBinIon]) {
								pdTmpRawData[iBinIon] = dIntensity;
							}
							if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity) {
								pPre->dHighestIntensity = pdTmpRawData[iBinIon];
							}
						}
					} else if (ProNovoConfig::options.iRemovePrecursor == 2) {
						int j;
						int bNotPrec = 1;
						for (j = 1; j <= pScoring->_spectrumInfoInternal.iChargeState; j++) {
							double dMZ;
							dMZ = (mstSpectrum->dParentMass + (j - 1) * PROTON_MASS) / (double) (j);
							if (fabs(dIon - dMZ) < ProNovoConfig::options.dRemovePrecursorTol) {
								bNotPrec = 0;
								break;
							}
						}
						if (bNotPrec) {
							if (dIntensity > pdTmpRawData[iBinIon])
								pdTmpRawData[iBinIon] = dIntensity;

							if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
								pPre->dHighestIntensity = pdTmpRawData[iBinIon];
						}
					} else // iRemovePrecursor==0
					{
						if (dIntensity > pdTmpRawData[iBinIon])
							pdTmpRawData[iBinIon] = dIntensity;

						if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
							pPre->dHighestIntensity = pdTmpRawData[iBinIon];
					}
				}
			}
		}
	}
	return true;
}

// pdTmpRawData now holds raw data, pdTmpCorrelationData is windowed data after this function
void CometSearch::MakeCorrData(double *pdTmpRawData, double *pdTmpCorrelationData, struct Query *pScoring,
		struct PreprocessStruct *pPre) {
	int i, ii, iBin, iWindowSize, iNumWindows = 10;
	double dMaxWindowInten, dTmp1, dTmp2;

	iWindowSize = (int) ((pPre->iHighestIon) / iNumWindows) + 1;

	for (i = 0; i < iNumWindows; i++) {
		dMaxWindowInten = 0.0;

		for (ii = 0; ii < iWindowSize; ii++) { // Find max inten. in window.
			iBin = i * iWindowSize + ii;
			if (iBin < pScoring->_spectrumInfoInternal.iArraySize) {
				if (pdTmpRawData[iBin] > dMaxWindowInten) {
					dMaxWindowInten = pdTmpRawData[iBin];
				}
			}
		}

		if (dMaxWindowInten > 0.0) {
			dTmp1 = 50.0 / dMaxWindowInten;
			dTmp2 = 0.05 * pPre->dHighestIntensity;

			for (ii = 0; ii < iWindowSize; ii++) { // Normalize to max inten. in window.
				iBin = i * iWindowSize + ii;
				if (iBin < pScoring->_spectrumInfoInternal.iArraySize) {
					if (pdTmpRawData[iBin] > dTmp2) {
						pdTmpCorrelationData[iBin] = pdTmpRawData[iBin] * dTmp1;
					}
				}
			}
		}
	}
}

// Smooth input data over 5 points.
bool CometSearch::Smooth(double *data, int iArraySize, double *pdTmpSmoothedSpectrum) {
	int i;

	data[0] = 0.0;
	data[1] = 0.0;
	data[iArraySize - 1] = 0.0;
	data[iArraySize - 2] = 0.0;

	for (i = 2; i < iArraySize - 2; i++) {
		// *0.0625 is same as divide by 16.
		pdTmpSmoothedSpectrum[i] = (data[i - 2] + 4.0 * data[i - 1] + 6.0 * data[i] + 4.0 * data[i + 1] + data[i + 2])
				* 0.0625;
	}

	memcpy(data, pdTmpSmoothedSpectrum, iArraySize * sizeof(double));

	return true;
}

// Run 2 passes through to pull out peaks.
bool CometSearch::PeakExtract(double *data, int iArraySize, double *pdTmpPeakExtracted) {
	int i, ii, iStartIndex, iEndIndex;
	double dStdDev, dAvgInten;

	// 1st pass, choose only peak greater than avg + dStdDev.
	for (i = 0; i < iArraySize; i++) {
		pdTmpPeakExtracted[i] = 0.0;
		dAvgInten = 0.0;

		iStartIndex = i - 50;
		if (i - 50 < 0)
			iStartIndex = 0;

		iEndIndex = i + 50;
		if (i + 50 > iArraySize - 1)
			iEndIndex = iArraySize - 1;

		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dAvgInten += (double) data[ii];
		dAvgInten /= iEndIndex - iStartIndex;

		dStdDev = 0.0;
		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dStdDev += (data[ii] - dAvgInten) * (data[ii] - dAvgInten);
		dStdDev = sqrt(dStdDev / (iEndIndex - iStartIndex + 1));

		if ((i > 0) && (i < iArraySize - 1)) {
			if (data[i] > (dAvgInten + dStdDev)) {
				pdTmpPeakExtracted[i] = data[i] - dAvgInten + dStdDev;
				data[i] = 0;     // Remove the peak before 2nd pass.
			}
		}
	}

	// 2nd pass, choose only peak greater than avg + 2*dStdDev.
	for (i = 0; i < iArraySize; i++) {
		dAvgInten = 0.0;

		iStartIndex = i - 50;
		if (i - 50 < 0)
			iStartIndex = 0;

		iEndIndex = i + 50;
		if (i + 50 > iArraySize - 1)
			iEndIndex = iArraySize - 1;

		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dAvgInten += (double) data[ii];
		dAvgInten /= iEndIndex - iStartIndex;

		dStdDev = 0.0;
		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dStdDev += (data[ii] - dAvgInten) * (data[ii] - dAvgInten);
		dStdDev = sqrt(dStdDev / (iEndIndex - iStartIndex + 1));

		if ((i > 0) && (i < iArraySize - 1)) {
			if (data[i] > (dAvgInten + 2 * dStdDev))
				pdTmpPeakExtracted[i] = data[i] - dAvgInten + dStdDev;
		}
	}

	memcpy(data, pdTmpPeakExtracted, (size_t) iArraySize * sizeof(double));

	return true;
}

// Pull out top # ions for intensity matching in search.
void CometSearch::GetTopIons(double *pdTmpRawData, struct msdata *pTmpSpData, int iArraySize) {
	int i, ii, iLowestIntenIndex = 0;
	double dLowestInten = 0.0, dMaxInten = 0.0;

	for (i = 0; i < iArraySize; i++) {
		if (pdTmpRawData[i] > dLowestInten) {
			(pTmpSpData + iLowestIntenIndex)->dIntensity = (double) pdTmpRawData[i];
			(pTmpSpData + iLowestIntenIndex)->dIon = (double) i;

			if ((pTmpSpData + iLowestIntenIndex)->dIntensity > dMaxInten)
				dMaxInten = (pTmpSpData + iLowestIntenIndex)->dIntensity;

			dLowestInten = (pTmpSpData + 0)->dIntensity;
			iLowestIntenIndex = 0;

			for (ii = 1; ii < NUM_SP_IONS; ii++) {
				if ((pTmpSpData + ii)->dIntensity < dLowestInten) {
					dLowestInten = (pTmpSpData + ii)->dIntensity;
					iLowestIntenIndex = ii;
				}
			}
		}
	}

	if (dMaxInten > FLOAT_ZERO) {
		for (i = 0; i < NUM_SP_IONS; i++)
			(pTmpSpData + i)->dIntensity = (((pTmpSpData + i)->dIntensity) / dMaxInten) * 100.0;
	}
}

int CometSearch::QsortByIon(const void *p0, const void *p1) {
	if (((struct msdata *) p1)->dIon < ((struct msdata *) p0)->dIon)
		return (1);
	else if (((struct msdata *) p1)->dIon > ((struct msdata *) p0)->dIon)
		return (-1);
	else
		return (0);
}

// Works on Sp data.
void CometSearch::StairStep(struct msdata *pTmpSpData, double dFragmentBinSize) {
	int i, ii, iii;
	double dMaxInten, dGap;

	i = 0;
	while (i < NUM_SP_IONS - 1) {
		ii = i;
		dMaxInten = (pTmpSpData + i)->dIntensity;
		dGap = 0.0;

		while (dGap <= dFragmentBinSize && ii < NUM_SP_IONS - 1) {
			ii++;
			dGap = (pTmpSpData + ii)->dIon - (pTmpSpData + ii - 1)->dIon;

			// Finds the max intensity for adjacent points.
			if (dGap <= dFragmentBinSize) {
				if ((pTmpSpData + ii)->dIntensity > dMaxInten)
					dMaxInten = (pTmpSpData + ii)->dIntensity;
			}
		}

		// Sets the adjacent points to the dMaxInten.
		for (iii = i; iii < ii; iii++)
			(pTmpSpData + iii)->dIntensity = dMaxInten;

		i = ii;
	}
}

void CometSearch::print(Query * pScoring) {
	int i = 0, x = 0, y = 0;
	cout << "\n" << endl;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		x = i / SPARSE_MATRIX_SIZE;
		if (pScoring->ppfSparseFastXcorrData[x] == NULL) {
			cout << "0\t";
		} else {
			y = i - (x * SPARSE_MATRIX_SIZE);
			cout << pScoring->ppfSparseFastXcorrData[x][y] << "\t";
		}
	}
	cout << endl;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		x = i / SPARSE_MATRIX_SIZE;
		if (pScoring->ppfSparseFastXcorrDataNL[x] == NULL) {
			cout << "0\t";
		} else {
			y = i - (x * SPARSE_MATRIX_SIZE);
			cout << pScoring->ppfSparseFastXcorrDataNL[x][y] << "\t";
		}
	}
	cout << endl;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		x = i / SPARSE_MATRIX_SIZE;
		if (pScoring->ppfSparseSpScoreData[x] == NULL) {
			cout << "0\t";
		} else {
			y = i - (x * SPARSE_MATRIX_SIZE);
			cout << pScoring->ppfSparseSpScoreData[x][y] << "\t";
		}
	}
	cout << "\n" << endl;
}

bool CometSearch::ScorePeptides(Peptide * currentPeptide, bool *pbDuplFragment, double * _pdAAforward,
		double * _pdAAreverse, MS2Scan * mstSpectrum, unsigned int *** _uiBinnedIonMasses, double & dXcorr, int test) {
	double dInverseBinWidth = 0, dOneMinusBinOffset = 0;
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dHighResOneMinusBinOffset;
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dLowResOneMinusBinOffset;
	}

	string * sSequence = &(currentPeptide->sPeptide);
	size_t i = 0;
	size_t j = 0;
	size_t k = sSequence->length() - 1;
	string currentPTM;
	size_t iPeptideLength = 0;
	size_t iPos = 0;
	map<string, double>::iterator iterResidueMonoMass;
	for (i = 0; i <= k; ++i) {
		if (isalpha(sSequence->at(i))) {
			iPeptideLength = iPeptideLength + 1;
		}
	}
	int iLenMinus1 = iPeptideLength - 1;
	if ((int) iPeptideLength < ProNovoConfig::getMinPeptideLength()) {
		cerr << "ERROR: Peptide sequence is too short " << sSequence << endl;
		return false;
	}

	if (sSequence->at(j) != '[') {
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}
	double dBion = ProNovoConfig::precalcMasses.dNtermProton;
	double dYion = ProNovoConfig::precalcMasses.dCtermOH2Proton;
	j++;
	if (!isalpha(sSequence->at(j))) {
		currentPTM = sSequence->at(j);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
			return false;
		}
		dBion += iterResidueMonoMass->second;
		j++;
	}
	if (sSequence->at(k) == ']') {
		k--;
	} else {
		currentPTM = sSequence->at(k);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
			return false;
		}
		dYion += iterResidueMonoMass->second;
		k--;
		if (sSequence->at(k) != ']') {
			cerr << "ERROR: (second) Last character in a peptide sequence must be ]." << endl;
			return false;
		}
		k--;
	}
	for (i = 0; i < (size_t) iLenMinus1; i++) {
		//First forward
		if (!isalpha(sSequence->at(j))) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			return false;
		}
		currentPTM = sSequence->at(j);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			return false;
		}
		dBion += iterResidueMonoMass->second;
		j++;
		if (!isalpha(sSequence->at(j))) {
			currentPTM = sSequence->at(j);
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
				return false;
			}
			dBion += iterResidueMonoMass->second;
			j++;
		}
		_pdAAforward[iPos] = dBion;
		//Now reverse
		if (!isalpha(sSequence->at(k))) {
			currentPTM = sSequence->at(k);
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find this PTM in the config file " << currentPTM << endl;
				return false;
			}
			dYion += iterResidueMonoMass->second;
			k--;
		}
		if (!isalpha(sSequence->at(k))) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			return false;
		}
		currentPTM = sSequence->at(k);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find this PTM in the config file" << currentPTM << endl;
			return false;
		}
		dYion += iterResidueMonoMass->second;
		k--;
		_pdAAreverse[iPos] = dYion;
		iPos++;
	}
	// Now get the set of binned fragment ions once to compare this peptide against all matching spectra.
	int ctCharge = 0;
	int ctIonSeries = 0;
	int iWhichIonSeries = 0;
	int ctLen = 0;
	int iVal = 0;
	Query* pQuery = mstSpectrum->pQuery;
	int iMaxFragCharge = pQuery->_spectrumInfoInternal.iMaxFragCharge;
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ctIonSeries = 0; ctIonSeries < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ctIonSeries++) {
			iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ctIonSeries];
			for (ctLen = 0; ctLen < iLenMinus1; ctLen++) {
				iVal = BINX(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse),
						dInverseBinWidth, dOneMinusBinOffset);
				pbDuplFragment[iVal] = false;
			}
		}
	}
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ctIonSeries = 0; ctIonSeries < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ctIonSeries++) {
			iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ctIonSeries];
			for (ctLen = 0; ctLen < iLenMinus1; ctLen++) {
				iVal = BINX(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse),
						dInverseBinWidth, dOneMinusBinOffset);
				if (pbDuplFragment[iVal] == false) {
					_uiBinnedIonMasses[ctCharge][ctIonSeries][ctLen] = iVal;
					pbDuplFragment[iVal] = true;
				} else {
					_uiBinnedIonMasses[ctCharge][ctIonSeries][ctLen] = 0;
				}
			}
		}
	}
	dXcorr = 0;
	bool bUseNLPeaks = false;
	float **ppSparseFastXcorrData;              // use this if bSparseMatrix
	int bin, x, y;
	int iMax = pQuery->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ctIonSeries = 0; ctIonSeries < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ctIonSeries++) {
			iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ctIonSeries];
			if (ProNovoConfig::ionInformation.bUseNeutralLoss
					&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A]
							|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
							|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
				bUseNLPeaks = true;
			} else {
				bUseNLPeaks = false;
			}
			if (ctCharge == 1 && bUseNLPeaks) {
				ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrDataNL;
			} else {
				ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrData;
			}
			for (ctLen = 0; ctLen < iLenMinus1; ctLen++) {
				//MH: newer sparse matrix converts bin to sparse matrix bin
				bin = _uiBinnedIonMasses[ctCharge][ctIonSeries][ctLen];
				x = bin / SPARSE_MATRIX_SIZE;
				if (ppSparseFastXcorrData[x] == NULL || x > iMax) // x should never be > iMax so this is just a safety check
					continue;
				y = bin - (x * SPARSE_MATRIX_SIZE);
				dXcorr += ppSparseFastXcorrData[x][y];
				//cout << x << "," << y << endl;
			}
		}
	}
	if (dXcorr < XCORR_CUTOFF) {
		dXcorr = XCORR_CUTOFF;
	} else {
		dXcorr *= 0.005;  // Scale intensities to 50 and divide score by 1E4.
	}

	/*cout << "\n" << endl;
	 for (i = 1; (int) i < pQuery->_spectrumInfoInternal.iArraySize; i++) {
	 x = i / SPARSE_MATRIX_SIZE;
	 if (pQuery->ppfSparseFastXcorrDataNL[x] == NULL) {
	 //cout << "0\t";
	 } else {
	 y = i - (x * SPARSE_MATRIX_SIZE);
	 if (pQuery->ppfSparseFastXcorrDataNL[x][y] != 0) {
	 cout << pQuery->ppfSparseFastXcorrDataNL[x][y] << "\t";
	 }
	 }
	 }
	 cout << endl;
	 for (i = 1; (int) i < pQuery->_spectrumInfoInternal.iArraySize; i++) {
	 x = i / SPARSE_MATRIX_SIZE;
	 if (pQuery->ppfSparseFastXcorrData[x] == NULL) {
	 //cout << "0\t";
	 } else {
	 y = i - (x * SPARSE_MATRIX_SIZE);
	 if (pQuery->ppfSparseFastXcorrData[x][y] != 0) {
	 cout << pQuery->ppfSparseFastXcorrData[x][y] << "\t";
	 }
	 }
	 }
	 cout << endl;*/

	return true;
}

double CometSearch::GetFragmentIonMass(int iWhichIonSeries, int i, int ctCharge, double *_pdAAforward,
		double *_pdAAreverse) {
	double dFragmentIonMass = 0.0;

	switch (iWhichIonSeries) {
	case ION_SERIES_B:
		dFragmentIonMass = _pdAAforward[i];
		break;
	case ION_SERIES_Y:
		dFragmentIonMass = _pdAAreverse[i];
		break;
	case ION_SERIES_A:
		dFragmentIonMass = _pdAAforward[i] - ProNovoConfig::precalcMasses.dCO;
		break;
	case ION_SERIES_C:
		dFragmentIonMass = _pdAAforward[i] + ProNovoConfig::precalcMasses.dNH3;
		break;
	case ION_SERIES_Z:
		dFragmentIonMass = _pdAAreverse[i] - ProNovoConfig::precalcMasses.dNH2;
		break;
	case ION_SERIES_X:
		dFragmentIonMass = _pdAAreverse[i] + ProNovoConfig::precalcMasses.dCOminusH2;
		break;
	}

	return (dFragmentIonMass + (ctCharge - 1) * PROTON_MASS) / ctCharge;
}

bool CometSearch::CalculateSP(double & fScoreSp, double* _pdAAforward, double * _pdAAreverse, MS2Scan * mstSpectrum,
		int iLenPeptide) {
	double dInverseBinWidth = 0, dOneMinusBinOffset = 0;
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dHighResOneMinusBinOffset;
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dLowResOneMinusBinOffset;
	}
	int iMaxFragCharge = mstSpectrum->pQuery->_spectrumInfoInternal.iMaxFragCharge;
	int iMatchedFragmentIonCt = 0;
	IonSeriesStruct ionSeries[9];
	double dTmpIntenMatch = 0.0;
	double dConsec = 0.0;
	int ii;
	for (ii = 0; ii < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ii++) {
		int iii;
		for (iii = 1; iii <= iMaxFragCharge; iii++)
			ionSeries[ProNovoConfig::ionInformation.piSelectedIonSeries[ii]].bPreviousMatch[iii] = 0;
	}
	int ctCharge;
	double dFragmentIonMass = 0.0;
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ii = 0; ii < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ii++) {
			int iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ii];

			// As both _pdAAforward and _pdAAreverse are increasing, loop through
			// iLenPeptide-1 to complete set of internal fragment ions.
			for (int iii = 0; iii < iLenPeptide - 1; iii++) {
				// Gets fragment ion mass.
				dFragmentIonMass = GetFragmentIonMass(iWhichIonSeries, iii, ctCharge, _pdAAforward, _pdAAreverse);

				if (!(dFragmentIonMass <= FLOAT_ZERO)) {
					int iFragmentIonMass = BINX(dFragmentIonMass, dInverseBinWidth, dOneMinusBinOffset);
					float fSpScore;

					fSpScore = FindSpScore(mstSpectrum->pQuery, iFragmentIonMass);

					if (fSpScore > FLOAT_ZERO) {
						iMatchedFragmentIonCt++;

						// Simple sum intensity.
						dTmpIntenMatch += fSpScore;

						// Increase score for consecutive fragment ion series matches.
						if (ionSeries[iWhichIonSeries].bPreviousMatch[ctCharge])
							dConsec += 0.075;

						ionSeries[iWhichIonSeries].bPreviousMatch[ctCharge] = 1;
					} else {
						ionSeries[iWhichIonSeries].bPreviousMatch[ctCharge] = 0;
					}
				}
			}
		}
	}
	fScoreSp = (double) ((dTmpIntenMatch * iMatchedFragmentIonCt * (1.0 + dConsec))
			/ ((iLenPeptide) * iMaxFragCharge * ProNovoConfig::ionInformation.iNumIonSeriesUsed));

	return true;
}

double CometSearch::FindSpScore(Query *pQuery, int bin) {
	int x = bin / SPARSE_MATRIX_SIZE;
	if (pQuery->ppfSparseSpScoreData[x] == NULL)
		return 0.0f;
	int y = bin - (x * SPARSE_MATRIX_SIZE);
	return pQuery->ppfSparseSpScoreData[x][y];
}
