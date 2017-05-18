/*
 * SiprosReader.cpp
 *
 *  Created on: May 16, 2017
 *      Author: xgo
 */

#include "Common.h"
#include "CometData.h"
#include "CometInterfaces.h"
#include "CometDataInternal.h"
#include <string>
#include <iostream>
#include "SiprosReader.h"

using namespace CometInterfaces;

bool SiprosReader::MzmlReader(string & filename_str, vector<Spectrum> *  _vSpectra){

	// std::string filename = "/home/xgo/Temp/OSU_D2_FASP_Elite_02262014_10.mzML";

	ICometSearchManager* pCometSearchMgr = GetCometSearchManager();

	bool bSearchSucceeded = pCometSearchMgr->ReadSpectrumData(filename_str.c_str(), _vSpectra);

	ReleaseCometSearchManager();

	return bSearchSucceeded;
}

