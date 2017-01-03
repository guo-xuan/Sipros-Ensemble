/*
 * ProteinDbParser.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: xgo
 */

#include "ProteinDbParser.h"

string ProteinDbParser::sOpenCharMutation = "((((";
string ProteinDbParser::sEndCharMutation = "))))";

void LookupTable::clean() {
	for (size_t i = 0; i < vSet.size(); i++) {
		vb.at(vSet.at(i)) = false;
	}
	vSet.clear();
}

bool LookupTable::exist(size_t _index) {
	if (_index >= vb.size()) {
		return false;
	}
	return vb.at(_index);
}

bool LookupTable::push(size_t _index) {
	if (_index >= vb.size()) {
		vb.resize(_index + 1, false);
	}
	if (vb.at(_index) == false) {
		vSet.push_back(_index);
		vb.at(_index) = true;
		return true;
	}
	return false;
}

ProteinDbParser::ProteinDbParser(bool _bScreenOutput) {
	imaxPTM = ProNovoConfig::getMaxPTMcount();
	iMaxMissedCleavageSites = ProNovoConfig::getMaxMissedCleavages();
	bScreenOutput = _bScreenOutput;
	sDBFilename = ProNovoConfig::getFASTAfilename();
	sResiduesBeforeCleavageSite = ProNovoConfig::getCleavageAfterResidues();
	sResiduesAfterCleavageSite = ProNovoConfig::getCleavageBeforeResidues();
	vector<string> vsallNames = ProNovoConfig::vsSingleResidueNames;
	slegalChar.clear();
	for (size_t i = 0; i < vsallNames.size(); i++) {
		if (isalpha(vsallNames.at(i).at(0))) {
			slegalChar.push_back(vsallNames.at(i).at(0));
		}
	}
	db_stream.open(sDBFilename.c_str());
	if (!(db_stream.is_open())) {
		cerr << "fail to open Database file " << sDBFilename << endl;
		exit(1);
	}
	getline(db_stream, sLine);
	if (sLine.at(0) == '>') {
		snextProteinName = sLine.substr(1, sLine.find_first_of(" \t\f\v\n\r") - 1);
		sMutationInfo.clear();
		iMutationInfoOpenPos = sLine.rfind(sOpenCharMutation);
		iMutationInfoClosePos = sLine.rfind(sEndCharMutation);
		if (iMutationInfoOpenPos != string::npos && iMutationInfoClosePos != string::npos) {
			sMutationInfo = sLine.substr(iMutationInfoOpenPos + sOpenCharMutation.length(), iMutationInfoClosePos - iMutationInfoOpenPos - sOpenCharMutation.length());
		}
	}
	iProteinId = 0;
	iMutationInfoOpenPos = 0;
	iMutationInfoClosePos = 0;

	orderstring = "[]";
	for (size_t i = 0; i < ProNovoConfig::vsSingleResidueNames.size(); i++) {
		if ((ProNovoConfig::vsSingleResidueNames.at(i).length() == 1)
				&& (isalpha(ProNovoConfig::vsSingleResidueNames[i][0]))) {
			orderstring += ProNovoConfig::vsSingleResidueNames[i];
		}
	}
	// Currently, in the peptide generating part, SIP and Regular modes are same.
	if (!ptmlist.populate_from_xml_config()) {
		// the current configure file is not xml format, but we still use this function
		cerr << "Error in parsing PTM rules from config " << endl;
		exit(1);
	}
	Initial_PTM_Map();

	dTerminusMassN = ProNovoConfig::getTerminusMassN();
	dTerminusMassC = ProNovoConfig::getTerminusMassC();

	for (size_t i = 0; i < orderstring.size(); i++) {
		if (orderstring.at(i) == '[' || orderstring.at(i) == ']') {
			mResiduleMass[orderstring.at(i)] = 0;
		} else {
			mResiduleMass[orderstring.at(i)] = ProNovoConfig::getResidueMass(orderstring.substr(i, 1));
		}
	}

	iMinPeptideLen = ProNovoConfig::getMinPeptideLength();
	iMaxPeptideLen = ProNovoConfig::getMaxPeptideLength();

	bSkipM = false;
}

ProteinDbParser::~ProteinDbParser() {
}

bool ProteinDbParser::getNextProtein() {
	sLine.clear();
	sProteinSeq.clear();
	scurrentProteinName = snextProteinName;
	// no more proteins in the database, if next protein name is empty
	if (snextProteinName.empty()) {
		db_stream.clear();
		db_stream.close();
		return false;
	}
	/*if(scurrentProteinName == "21942184((((105E))))"){
		cout << "check" << endl;
	}*/
	snextProteinName.clear();
	iMutationInfoOpenPos = 0;
	iMutationInfoClosePos = 0;
	while (!db_stream.eof()) {
		sLine.clear();
		getline(db_stream, sLine);
		if (sLine == "") {
			continue;
		}
		if (sLine.at(0) == '>') {
			if ((ProNovoConfig::getTestStartRemoval()) && (sProteinSeq.at(0) == 'M')) {
				sProteinSeq = sProteinSeq.substr(1);
				bSkipM = true;
			} else {
				bSkipM = false;
			}
			ParseMutationInfo(sMutationInfo, vAllMutations, iAllMutationsSize, vviAllMutationLinks,
					iAllMutationLinksSize, miiMutations);
			ParseCleavageSites(sProteinSeq, vAllMutations, iAllMutationsSize, miiMutations, vAllPosibleCleavageSites,
					iAllPosibleCleavageSitesSize, vStaticCleavageSites, vDynamicCleavageSites);
			snextProteinName = sLine.substr(1, sLine.find_first_of(" \t\f\v\n\r") - 1);
			// extract the mutation information
			sMutationInfo.clear();
			iMutationInfoOpenPos = sLine.rfind(sOpenCharMutation);
			iMutationInfoClosePos = sLine.rfind(sEndCharMutation);
			if (iMutationInfoOpenPos != string::npos && iMutationInfoClosePos != string::npos) {
				sMutationInfo = sLine.substr(iMutationInfoOpenPos + sOpenCharMutation.length(),
						iMutationInfoClosePos - iMutationInfoOpenPos - sOpenCharMutation.length());
			}
			break;
		} else {
			RemoveIllegalResidue(sLine);
			sProteinSeq.append(sLine);
		}
	}
	if (db_stream.eof()) {
		if ((ProNovoConfig::getTestStartRemoval()) && (sProteinSeq.at(0) == 'M')) {
			sProteinSeq = sProteinSeq.substr(1);
			bSkipM = true;
		} else {
			bSkipM = false;
		}
		ParseMutationInfo(sMutationInfo, vAllMutations, iAllMutationsSize, vviAllMutationLinks, iAllMutationLinksSize,
				miiMutations);
		ParseCleavageSites(sProteinSeq, vAllMutations, iAllMutationsSize, miiMutations, vAllPosibleCleavageSites,
				iAllPosibleCleavageSitesSize, vStaticCleavageSites, vDynamicCleavageSites);
	}
	iProteinId++;
	if (bScreenOutput) {
		cout << "Processing protein #" << iProteinId << "\r";
	}
	return true;
}

size_t ProteinDbParser::getNumMutatedPeptide() {
	return iMutatedPeptidesSize;
}

size_t ProteinDbParser::getNumPtmPeptide() {
	return iPtmPeptidesSize;
}

const vector<MutatedPeptide> & ProteinDbParser::peptideModify(size_t _iMutatedPeptideIndex) {
	iPtmPeptidesSize = 0;
	generatePtmPeptide(vMutatedPeptides.at(_iMutatedPeptideIndex), vPtmPeptides, iPtmPeptidesSize, ptm_position_all,
			comb_order, ptm_order, ele_num);
	return vPtmPeptides;
}

const vector<MutatedPeptide> & ProteinDbParser::proteinDigest() {
	iMutatedPeptidesSize = 0;
	digest(sProteinSeq, vAllPosibleCleavageSites, iAllPosibleCleavageSitesSize, vStaticCleavageSites,
			vDynamicCleavageSites, vCombinationCleavageSites, vSelectedDynamicCleavageSites, vProteinSegments,
			iProteinSegmentsSize);
	for (size_t i = 0; i < iProteinSegmentsSize; i++) {
		generateMutatedPeptides(vProteinSegments.at(i), vMutatedPeptides, iMutatedPeptidesSize, miiMutations,
				vAllMutations, vviAllMutationLinks, lookupTablePosition, vLinkMasks, vPositionMasks,
				vMutationsInPeptide);
	}
	return vMutatedPeptides;
}

void ProteinDbParser::setPeptide(size_t _index, Peptide * _pPeptide) {
	MutatedPeptide * p = &(vPtmPeptides.at(_index));
	string s;
	s.append("[");
	s.append(p->proteinSegment->sSeq);
	s.append("]");
	_pPeptide->setPeptide(p->sSeq, s, scurrentProteinName, p->proteinSegment->pPostion.first, p->dcurrentMass,
			p->cIdentifyPrefix, p->cIdentifySuffix, p->cOriginalPrefix, p->cOriginalSuffix);
}

// private methods for ProteinDbParser

void ProteinDbParser::applyMutations(ProteinSegment & _proteinsegment, vector<Mutation> & _vAllVariants,
		vector<size_t> & _vIndexDynamicMutations, int64_t _iMasks, vector<MutatedPeptide> & _vMutatedPeptides,
		size_t & _iMutatedPeptidesSize, LookupTable & _lookupTablePosition) {
	/*	if (_proteinsegment.sSeq == "GLEIVKELLKK") {
	 cout << "check" << endl;
	 }*/
	if (_iMutatedPeptidesSize >= _vMutatedPeptides.size()) {
		_vMutatedPeptides.push_back(MutatedPeptide());
	}
	MutatedPeptide * mutatedPeptide = &(_vMutatedPeptides.at(_iMutatedPeptidesSize));
	mutatedPeptide->sSeq.clear();
	mutatedPeptide->sSeq.append("[");
	mutatedPeptide->sSeq.append(_proteinsegment.sSeq);
	mutatedPeptide->sSeq.append("]");
	mutatedPeptide->proteinSegment = &(_proteinsegment);
	size_t beg = _proteinsegment.pPostion.first, end = _proteinsegment.pPostion.second;
	Mutation * mutation;
	int i = 0;
	if (!_vIndexDynamicMutations.empty()) {
		for (i = _vIndexDynamicMutations.size() - 1; i >= 0; i--) {
			if ((_iMasks & 1) == 1) {
				mutation = &(_vAllVariants.at(_vIndexDynamicMutations.at(i)));
				if (mutation->cType == 'S') {
					mutatedPeptide->sSeq.at(mutation->iPos - beg + 1) = mutation->sSeq.at(0);
				} else {
					cout << "Error: mutations other than substitution not allowed right now." << endl;
					exit(1);
				}
			}
			_iMasks = _iMasks >> 1;
		}
	}
	// check the length of peptide
	if (mutatedPeptide->sSeq.length() - 2 < iMinPeptideLen || mutatedPeptide->sSeq.length() - 2 > iMaxPeptideLen) {
		return;
	}
	// check the cleavage sites are all correct
	_lookupTablePosition.clean();
	for (size_t j = 0; j < _proteinsegment.vPositiveDynamicCleavageSites.size(); j++) {
		_lookupTablePosition.push(_proteinsegment.vPositiveDynamicCleavageSites.at(j));
		if (!isCleavageSite(mutatedPeptide->sSeq.at(_proteinsegment.vPositiveDynamicCleavageSites.at(j) - beg),
				mutatedPeptide->sSeq.at(_proteinsegment.vPositiveDynamicCleavageSites.at(j) - beg + 1))) {
			// if cleavage is the begin or the end of a protein, the site must be a cleavage site
			if (_proteinsegment.vPositiveDynamicCleavageSites.at(j) == 0
					|| _proteinsegment.vPositiveDynamicCleavageSites.at(j) == sProteinSeq.length()) {
				continue;
			}
			return;
		}
	}
	map<size_t, CleavageSite*>::iterator it, itLow, itUp;
	itLow = mscCleavageSites.lower_bound(beg);
	itUp = mscCleavageSites.upper_bound(end);
	for (it = itLow; it != itUp; ++it) {
		if (_lookupTablePosition.push(it->second->iPos)) {
			if (isCleavageSite(mutatedPeptide->sSeq.at(it->second->iPos - beg),
					mutatedPeptide->sSeq.at(it->second->iPos - beg + 1))) {
				return;
			}
		}
	}
	/*	if (mutatedPeptide->sSeq == "[ELLKSADVTVENMAPGTIER]") {
	 cout << "check 2" << endl;
	 }*/
	mutatedPeptide->cIdentifyPrefix = _proteinsegment.cIdentifyPrefix;
	mutatedPeptide->cIdentifySuffix = _proteinsegment.cIdentifySuffix;
	mutatedPeptide->cOriginalPrefix = _proteinsegment.cOriginalPrefix;
	mutatedPeptide->cOriginalSuffix = _proteinsegment.cOriginalSuffix;
	// set the mass of this mutatedPeptide
	mutatedPeptide->dcurrentMass = dTerminusMassC + dTerminusMassN;
	for (size_t j = 1; j < mutatedPeptide->sSeq.size() - 1; j++) {
		mutatedPeptide->dcurrentMass += mResiduleMass[mutatedPeptide->sSeq.at(j)];
	}
	// cout << mutatedPeptide->sSeq << endl;
	/*	if (mutatedPeptide->sSeq == "[ELLKSADVTVENMAPGTIER]") {
	 cout << "check 2" << endl;
	 }*/
	++_iMutatedPeptidesSize;
}

bool comparePtrToNode(CleavageSite* a, CleavageSite* b) {
	return (a->iPos < b->iPos);
}

void ProteinDbParser::digest(string & _seq, vector<CleavageSite> & _vAllPosibleCleavageSites,
		size_t _iAllPosibleCleavageSitesSize, vector<size_t> & _vStaticCleavageSites,
		vector<size_t> & _vDynamicCleavageSites, vector<CleavageSite*> & _vCombinationCleavageSites,
		vector<size_t> & _vSelectedDynamicCleavageSites, vector<ProteinSegment> & _vProteinSegments,
		size_t & _iProteinSegmentsSize) {
	int iMissCleavage = 0;
	size_t iNumPermutation = 0;
	size_t mod = 0, remainder = 0;
	_iProteinSegmentsSize = 0;
	// if there is no dynamic cleavage sites
	if (_vDynamicCleavageSites.size() == 0) {
		_vCombinationCleavageSites.clear();
		_vSelectedDynamicCleavageSites.clear();
		for (size_t i = 0; i < _iAllPosibleCleavageSitesSize; i++) {
			_vCombinationCleavageSites.push_back(&(_vAllPosibleCleavageSites.at(i)));
		}
		generateProteinSegment(_seq, _vCombinationCleavageSites, _vSelectedDynamicCleavageSites, _vProteinSegments,
				_iProteinSegmentsSize);
	} else {
		// if there are dynamic cleavage sites
		// try all the static cleavage sites first
		_vCombinationCleavageSites.clear();
		_vSelectedDynamicCleavageSites.clear();
		for (size_t i = 0; i < _vStaticCleavageSites.size(); i++) {
			_vCombinationCleavageSites.push_back(&(_vAllPosibleCleavageSites.at(_vStaticCleavageSites.at(i))));
		}
		generateProteinSegment(_seq, _vCombinationCleavageSites, _vSelectedDynamicCleavageSites, _vProteinSegments,
				_iProteinSegmentsSize);
		// fix the static cleavage sites then try every dynamic cleavage site
		for (int i = 0; i < ((int)_vDynamicCleavageSites.size()); i++) {
			int j = i + 1;
			for (; j < ((int)_vDynamicCleavageSites.size()); j++) {
				iMissCleavage = _vDynamicCleavageSites.at(j) - _vDynamicCleavageSites.at(i) - 1;
				iMissCleavage -= j - i - 1;
				if (iMissCleavage > iMaxMissedCleavageSites) {
					break;
				}
			}
			// permutation
			--j;
			iNumPermutation = pow(2, j - i);
			for (size_t k = 0; k < iNumPermutation; k++) {
				_vCombinationCleavageSites.clear();
				_vSelectedDynamicCleavageSites.clear();
				_vSelectedDynamicCleavageSites.push_back(_vDynamicCleavageSites.at(i));
				_vCombinationCleavageSites.push_back(&(_vAllPosibleCleavageSites.at(_vDynamicCleavageSites.at(i))));
				// decode the permutation
				remainder = k;
				for (int l = 1; l <= j - i; l++) {
					mod = remainder % 2;
					remainder = remainder >> 1;
					if (mod == 1) {
						_vCombinationCleavageSites.push_back(
								(&(_vAllPosibleCleavageSites.at(_vDynamicCleavageSites.at(l + i)))));
						_vSelectedDynamicCleavageSites.push_back(_vDynamicCleavageSites.at(l + i));
					}
				}
				// some dynamics are not selected
				iMissCleavage = _vSelectedDynamicCleavageSites.back() - _vSelectedDynamicCleavageSites.front() - 1;
				iMissCleavage -= _vSelectedDynamicCleavageSites.size() - 2;
				if (iMissCleavage > iMaxMissedCleavageSites) {
					continue;
				}
				// generate protein segments
				for (int l = 0; l < ((int)_vStaticCleavageSites.size()); l++) {
					_vCombinationCleavageSites.push_back(
							(&(_vAllPosibleCleavageSites.at(_vStaticCleavageSites.at(l)))));
				}
				sort(_vCombinationCleavageSites.begin(), _vCombinationCleavageSites.end(), comparePtrToNode);
				_vSelectedDynamicCleavageSites.clear();
				// each position will only have one cleavage site
				for (int l = 0; l < ((int)_vCombinationCleavageSites.size()); l++) {
					if (_vCombinationCleavageSites.at(l)->bIsDynamicCleavageSite) {
						_vSelectedDynamicCleavageSites.push_back(l);
					}
				}
				generateProteinSegment(_seq, _vCombinationCleavageSites, _vSelectedDynamicCleavageSites,
						_vProteinSegments, _iProteinSegmentsSize);
			}
		}
	}
}

/**
 * mutation vector: 2, 3, 4, 5, 6
 * linked mutation: 3, 5
 * mask in binary:  0, 1, 0, 1, 0
 */
void ProteinDbParser::generateLinkMasks(vector<vector<size_t>> & _vviAllMutationLinks,
		const vector<size_t> & _vLinkIndexes, vector<Mutation> & _vAllVariants,
		vector<size_t> & _vIndexDynamicMutations, vector<int64_t> & _vMasks) {
	_vMasks.clear();
	vector<size_t>::iterator first1, first2, last1, last2;
	vector<size_t> * mutationLink;
	for (size_t i = 0; i < _vLinkIndexes.size(); i++) {
		int64_t iMask = 0;
		first1 = _vIndexDynamicMutations.begin();
		last1 = _vIndexDynamicMutations.end();
		mutationLink = &(_vviAllMutationLinks.at(_vLinkIndexes.at(i)));
		first2 = mutationLink->begin();
		last2 = mutationLink->end();
		while (first1 != last1 && first2 != last2) {
			if (*first1 < *first2) {
				++first1;
				iMask = iMask << 1;
			} else if (*first2 < *first1)
				++first2;
			else {
				iMask = iMask << 1;
				iMask += 1;
				++first1;
				++first2;
			}
		}
		while (first1 != last1) {
			iMask = iMask << 1;
			++first1;
		}
		_vMasks.push_back(iMask);
	}

}

void ProteinDbParser::generateMutatedPeptides(ProteinSegment & _proteinsegment,
		vector<MutatedPeptide> & _vMutatedPeptides, size_t & _iMutatedPeptidesSize,
		multimap<size_t, size_t> & _miiMutations, vector<Mutation> & _vAllMutations,
		vector<vector<size_t>> & _vviAllMutationLinks, LookupTable & _lookupTablePosition,
		vector<int64_t> & _vLinkMasks, vector<int64_t> & _vPositionMasks, vector<size_t> & _vDynamicMutations) {
	// get all mutations with static mutations
	vector<size_t> * mutationLink;
	Mutation * variant, *linkedMutation;
	_lookupTablePosition.clean();
	// get all the mutations in the range of this protein segment;
	// considering same position problem (mutation level, link level)
	multimap<size_t, size_t>::iterator it, itlow, itup;
	itlow = _miiMutations.lower_bound(_proteinsegment.pPostion.first);
	itup = _miiMutations.upper_bound(_proteinsegment.pPostion.second - 1);
	_vDynamicMutations.clear();
	for (it = itlow; it != itup; ++it) {
		variant = &(_vAllMutations.at((*it).second));
		_vDynamicMutations.push_back((*it).second);
	}
	bool bConflict = false;
	LookupTable * lookupTableLinkIndex = &(_lookupTablePosition);
	lookupTableLinkIndex->clean();
	int iMutationCount = 0;
	for (size_t i = 0; i < _vDynamicMutations.size(); i++) {
		variant = &(_vAllMutations.at(_vDynamicMutations.at(i)));
		bConflict = false;
		iMutationCount = 0;
		if (variant->iMutationLinkIndex != EMPTY) {
			mutationLink = &(_vviAllMutationLinks.at(variant->iMutationLinkIndex));
			for (size_t j = 0; j < mutationLink->size(); j++) {
				linkedMutation = &(_vAllMutations.at(mutationLink->at(j)));
				if (linkedMutation->iPos >= _proteinsegment.pPostion.first
						&& linkedMutation->iPos < _proteinsegment.pPostion.second) {
					++iMutationCount;
				}
			}
		}
		if (iMutationCount > 1) {
			lookupTableLinkIndex->push(variant->iMutationLinkIndex);
		}
	}
	const vector<size_t> & vLinkIndexes = lookupTableLinkIndex->get();
	sort(_vDynamicMutations.begin(), _vDynamicMutations.end());
	// generate masks
	generateLinkMasks(_vviAllMutationLinks, vLinkIndexes, _vAllMutations, _vDynamicMutations, _vLinkMasks);
	generatePositionMasks(_vAllMutations, _vDynamicMutations, _vPositionMasks);
	// permutation of mutations
	int64_t iMask = 0, iTemp;
	int64_t iTotalMuations = pow(2, _vDynamicMutations.size());
	for (; iMask < iTotalMuations; iMask++) {
		// check if comply with all link constrains
		bConflict = false;
		for (size_t i = 0; i < _vLinkMasks.size(); i++) {
			iTemp = iMask & _vLinkMasks.at(i);
			if (iTemp != 0 && iTemp != _vLinkMasks.at(i)) {
				bConflict = true;
				break;
			}
		}
		if (bConflict) {
			continue;
		}
		// check if comply with position constrains
		for (size_t i = 0; i < _vPositionMasks.size(); i++) {
			iTemp = iMask & _vPositionMasks.at(i);
			if (iTemp == _vPositionMasks.at(i)) {
				bConflict = true;
				break;
			}
		}
		if (bConflict) {
			continue;
		}
		// now generate peptide with mutations
		applyMutations(_proteinsegment, _vAllMutations, _vDynamicMutations, iMask, _vMutatedPeptides,
				_iMutatedPeptidesSize, _lookupTablePosition);
	}
}

// unfinished
/**
 * _seq original sequence
 * _vCleavageSites all the cleavage sites
 * _vCleavageSitesMustHave the protein segment must contain these cleavage sites (index in _vCleavageSites)
 * _vProteinSegments store the protein segments
 */
void ProteinDbParser::generateProteinSegment(string & _seq, vector<CleavageSite*> & _vCleavageSites,
		vector<size_t> & _vCleavageSitesMustHave, vector<ProteinSegment> & _vProteinSegments,
		size_t & _iProteinSegmentsSize) {
	// _iProteinSegmentsSize = 0;
	size_t iBeginCleavePos = -1;
	size_t iEndCleavePos = 0;
	ProteinSegment * pProteinSegment;
	if (_vCleavageSitesMustHave.size() == 0) {
		// enumerate all the peptides
		for (int i = 0; i <= iMaxMissedCleavageSites; i++) {
			for (size_t j = 0; j + 1 + i < _vCleavageSites.size(); j++) {
				iBeginCleavePos = _vCleavageSites.at(j)->iPos;
				iEndCleavePos = _vCleavageSites.at(j + 1 + i)->iPos;
				if (_iProteinSegmentsSize >= _vProteinSegments.size()) {
					_vProteinSegments.push_back(ProteinSegment());
				}
				pProteinSegment = &(_vProteinSegments.at(_iProteinSegmentsSize));
				pProteinSegment->sSeq = _seq.substr(iBeginCleavePos, iEndCleavePos - iBeginCleavePos);
				pProteinSegment->pPostion.first = iBeginCleavePos;
				pProteinSegment->pPostion.second = iEndCleavePos;
				pProteinSegment->cIdentifyPrefix = (iBeginCleavePos == 0 ? '-' : _seq.at(iBeginCleavePos - 1));
				pProteinSegment->cOriginalPrefix = pProteinSegment->cIdentifyPrefix;
				pProteinSegment->cIdentifySuffix = (iEndCleavePos == _seq.size() ? '-' : _seq.at(iEndCleavePos));
				pProteinSegment->cOriginalSuffix = pProteinSegment->cIdentifySuffix;
				pProteinSegment->vPositiveDynamicCleavageSites.clear();
				for (size_t l = j; l <= j + 1 + i; l++) {
					pProteinSegment->vPositiveDynamicCleavageSites.push_back(_vCleavageSites.at(l)->iPos);
				}
				++_iProteinSegmentsSize;
			}
		}
	} else {
		// enumerate the peptides that contain all cleavage sites in _vCleavageSitesMustHave
		for (int i = 0; i <= iMaxMissedCleavageSites; i++) {
			int j = ((int) _vCleavageSitesMustHave.at(0)) - i - 1;
			if (j < 0) {
				j = 0;
			}
			for (; j <= (int) _vCleavageSitesMustHave.at(0) && (j + 1 + i) < ((int) _vCleavageSites.size()); j++) {
				// check if the new peptide contains all cleavage sites that we must have
				if (j + 1 + i < (int) _vCleavageSitesMustHave.back()) {
					continue;
				}
				iBeginCleavePos = _vCleavageSites.at(j)->iPos;
				iEndCleavePos = _vCleavageSites.at(j + 1 + i)->iPos;
				if (_iProteinSegmentsSize >= _vProteinSegments.size()) {
					_vProteinSegments.push_back(ProteinSegment());
				}
				pProteinSegment = &(_vProteinSegments.at(_iProteinSegmentsSize));
				pProteinSegment->sSeq = _seq.substr(iBeginCleavePos, iEndCleavePos - iBeginCleavePos);
				pProteinSegment->pPostion.first = iBeginCleavePos;
				pProteinSegment->pPostion.second = iEndCleavePos;
				pProteinSegment->cOriginalPrefix = (iBeginCleavePos == 0 ? '-' : _seq.at(iBeginCleavePos - 1));
				if (_vCleavageSites.at(j)->bIsDynamicCleavageSite && (!_vCleavageSites.at(j)->vMutations.empty())
						&& (_vCleavageSites.at(j)->iPos - 1 == _vCleavageSites.at(j)->vMutations.at(0)->iPos)) {
					if (_vCleavageSites.at(j)->vMutations.at(0)->cType == SUBSITUE) {
						pProteinSegment->cIdentifyPrefix = _vCleavageSites.at(j)->vMutations.at(0)->sSeq.at(0);
					} else {
						cout << "Error in ProteinDbParser::generateProteinSegment" << endl;
						exit(1);
					}
				} else {
					pProteinSegment->cIdentifyPrefix = pProteinSegment->cOriginalPrefix;
				}
				pProteinSegment->cOriginalSuffix = (iEndCleavePos == _seq.size() ? '-' : _seq.at(iEndCleavePos));
				if (_vCleavageSites.at(j)->bIsDynamicCleavageSite && (!_vCleavageSites.at(j)->vMutations.empty())
						&& (_vCleavageSites.at(j)->iPos == _vCleavageSites.at(j)->vMutations.back()->iPos)) {
					if (_vCleavageSites.at(j)->vMutations.at(0)->cType == SUBSITUE) {
						pProteinSegment->cIdentifySuffix = _vCleavageSites.at(j)->vMutations.back()->sSeq.at(0);
					} else {
						cout << "Error in ProteinDbParser::generateProteinSegment" << endl;
						exit(1);
					}
				} else {
					pProteinSegment->cIdentifySuffix = pProteinSegment->cOriginalSuffix;
				}

				pProteinSegment->vPositiveDynamicCleavageSites.clear();
				for (int l = j; l <= j + 1 + i; l++) {
					pProteinSegment->vPositiveDynamicCleavageSites.push_back(_vCleavageSites.at(l)->iPos);
				}
				++_iProteinSegmentsSize;
			}
		}
	}
}

void ProteinDbParser::generatePtmPeptide(MutatedPeptide & _mutatedPeptide, vector<MutatedPeptide> & _vPtmPeptides,
		size_t & _iPtmPeptidesSize, vector<pair<size_t, vector<pair<string, double> > *> > & _ptm_position_all,
		vector<int> & _comb_order, vector<int> & _ptm_order, vector<int> & _ele_num) {
	// without PTM
	MutatedPeptide * mutatedPeptide;
	if (_iPtmPeptidesSize >= _vPtmPeptides.size()) {
		_vPtmPeptides.push_back(MutatedPeptide());
	}
	mutatedPeptide = &(_vPtmPeptides.at(_iPtmPeptidesSize));
	(*mutatedPeptide) = _mutatedPeptide;
	++_iPtmPeptidesSize;
	// with PTM
	string sPtm;
	int i = 0, ori_val, j = 0, k = 0, l = 0, total_num, iPtmPos;
	// find all positions with PTMs
	size_t residue_id;
	pair<int, vector<pair<string, double> >*> ptm_position;
	_ptm_position_all.clear();
	for (j = 0; j < ((int) (_mutatedPeptide.sSeq.length())); j++) {
		// identify positions on which ptm may happen
		residue_id = orderstring.find(_mutatedPeptide.sSeq.at(j));
		if (!(ptm_map.at(residue_id).empty())) {
			//cout<<residue_id<<endl;
			ptm_position = make_pair(j, &(ptm_map.at(residue_id)));
			_ptm_position_all.push_back(ptm_position);
		}
	}
	total_num = _ptm_position_all.size();
	int icurrentMaxPtm = ((imaxPTM < ((ptm_position_all.size()))) ? imaxPTM : ((ptm_position_all.size())));
	int iPtmCount = 1;
	bool bNextCom = true;
	for (; iPtmCount <= icurrentMaxPtm; iPtmCount++) {
		// Initialize the ptm info
		_comb_order.clear();
		for (i = 0; i < (int) iPtmCount; i++) {
			_comb_order.push_back(i);
		}
		bNextCom = true;
		_comb_order.back() -= 1;
		while (bNextCom) {
			bNextCom = false;
			// comb order update
			for (i = (iPtmCount - 1); i >= 0; i--) {
				if (_comb_order.at(i) < total_num - iPtmCount + i) {
					ori_val = _comb_order[i];
					for (j = i; j < iPtmCount; j++) {
						_comb_order.at(j) = ori_val + 1 + j - i;
					}
					bNextCom = true;
					_ptm_order.clear();
					//ptm_order is for which ptm happen on positions
					_ele_num.clear();
					//ele_num is for numbers of ptm may happen on positions
					for (j = 0; j < ((int) comb_order.size()); j++) {
						_ptm_order.push_back(0);
						_ele_num.push_back((int) ptm_position_all.at(comb_order.at(j)).second->size());
					}
					_ptm_order.front() -= 1;
					// next ptm
					for (j = 0; j < ((int) _ptm_order.size()); j++) {
						for (l = j; l >= 0; l--) {
							while (_ptm_order.at(l) < (_ele_num.at(l) - 1)) {
								_ptm_order.at(l) += 1;
								for (k = 0; k < l; k++) {
									_ptm_order.at(k) = 0;
								}
								// apply ptm
								if (_iPtmPeptidesSize >= _vPtmPeptides.size()) {
									_vPtmPeptides.push_back(MutatedPeptide());
								}
								mutatedPeptide = &(_vPtmPeptides.at(_iPtmPeptidesSize));
								(*mutatedPeptide) = _mutatedPeptide;
								for (k = (_comb_order.size() - 1); k >= 0; k--) {
									iPtmPos = _ptm_position_all.at(_comb_order.at(k)).first;
									sPtm = _ptm_position_all.at(_comb_order.at(k)).second->at(_ptm_order.at(k)).first;
									mutatedPeptide->dcurrentMass += _ptm_position_all.at(_comb_order.at(k)).second->at(
											_ptm_order.at(k)).second;
									mutatedPeptide->sSeq.insert(iPtmPos + 1, sPtm);
								}
								// cout << mutatedPeptide->sSeq << endl;
								++_iPtmPeptidesSize;
							}
						}
					}
					break;
				}
			}
		}
	}

}

/**
 * mutation: 	2, 3, 4, 5, 6, 7
 * 3, 4 in the same position
 * 5, 6 in the same position
 * two masks:
 * 				0, 1, 1, 0, 0, 0
 * 				0, 0, 0, 1, 1, 0
 */
void ProteinDbParser::generatePositionMasks(vector<Mutation> & _vAllMutations, vector<size_t> & _vIndexDynamicMutations,
		vector<int64_t> & _vMasks) {
	int64_t iMask = 0;
	_vMasks.clear();
	size_t i = 0;
	if (_vIndexDynamicMutations.size() < 2) {
		return;
	}
	for (i = 0; i < _vIndexDynamicMutations.size() - 1; i++) {
		if (_vAllMutations.at(_vIndexDynamicMutations.at(i)).iPos
				== _vAllMutations.at(_vIndexDynamicMutations.at(i + 1)).iPos) {
			iMask += 1;
			iMask = iMask << 1;
		} else {
			if (iMask > 0) {
				iMask += 1;
				iMask = iMask << (_vIndexDynamicMutations.size() - i - 1);
				_vMasks.push_back(iMask);
				iMask = 0;
			}
		}
	}
	if (iMask > 0) {
		if (_vAllMutations.at(_vIndexDynamicMutations.at(i - 1)).iPos
				== _vAllMutations.at(_vIndexDynamicMutations.at(i)).iPos) {
			iMask += 1;
		} else {
			iMask = iMask << 1;
		}
		_vMasks.push_back(iMask);
	}
}

bool ProteinDbParser::isCleavageResidueAfter(char c2) {
	return (sResiduesAfterCleavageSite.find(c2, 0) != string::npos);
}

bool ProteinDbParser::isCleavageResidueBefore(char c1) {
	return (sResiduesBeforeCleavageSite.find(c1, 0) != string::npos);
}

bool ProteinDbParser::isCleavageSite(char c1, char c2) {
	if (c1 == '[') {
		return isCleavageResidueAfter(c2);
	}
	if (c2 == ']') {
		return isCleavageResidueBefore(c1);
	}
	return ((sResiduesBeforeCleavageSite.find(c1, 0) != string::npos)
			&& (sResiduesAfterCleavageSite.find(c2, 0) != string::npos));
}

bool ProteinDbParser::isCleavageSite(string & _seq, Mutation & _mutation, size_t _iCleavagePos) {
	if (_mutation.cType == SUBSITUE) {
		if (_mutation.iPos == _iCleavagePos - 1) {
			return isCleavageSite(_mutation.sSeq.at(0), _seq.at(_iCleavagePos));
		} else if (_mutation.iPos == _iCleavagePos) {
			return isCleavageSite(_seq.at(_iCleavagePos - 1), _mutation.sSeq.at(0));
		} else {
			cout << "Error: bool isCleavageSite(string & _seq, Mutation & _mutation, size_t _iCleavagePos)" << endl;
			exit(1);
		}
	} else {
		cout << "Error: bool isCleavageSite(string & _seq, Mutation & _mutation, size_t _iCleavagePos)" << endl;
		cout << "Indel are under construction." << endl;
		exit(1);
	}
}

bool ProteinDbParser::isCleavageSite(string & _seq, Mutation & _mutationBefore, Mutation & _mutationAfter,
		size_t _iCleavagePos) {
	if (_mutationBefore.cType == SUBSITUE && _mutationAfter.cType == SUBSITUE) {
		if (_mutationBefore.iPos == _iCleavagePos - 1 && _mutationAfter.iPos == _iCleavagePos) {
			return isCleavageSite(_mutationBefore.sSeq.at(0), _mutationAfter.sSeq.at(0));
		} else {
			cout
					<< "Error: bool ProteinDbParser::isCleavageSite(string & _seq, Mutation & _mutationBefore, Mutation & _mutationAfter, size_t _iCleavagePos)"
					<< endl;
			exit(1);
		}
	} else {
		cout
				<< "Error: bool ProteinDbParser::isCleavageSite(string & _seq, Mutation & _mutationBefore, Mutation & _mutationAfter, size_t _iCleavagePos)"
				<< endl;
		cout << "Indel are under construction." << endl;
		exit(1);
	}
}

// organize ptm information
void ProteinDbParser::Initial_PTM_Map() {
	vector<pair<string, double> > ptm_elm;
	pair<string, double> cur_elm;
	int i, residue_id;

	ptm_elm.clear();
	//cout<< orderstring.length() <<endl;
	for (i = 0; i < (int) orderstring.length(); i++) { //initialize the vector
		ptm_map.push_back(ptm_elm);
	}
	for (i = 0; i < ptmlist.size(); i++) {
		// cout<<orderstring<<"aa"<<endl;
		// cout<<ptmlist.residue(i)<<"bb"<<endl;
		residue_id = orderstring.find(ptmlist.residue(i));
		cur_elm = make_pair<string, double>(ptmlist.symbol(i), ptmlist.mass_shift(i));
		ptm_map[residue_id].push_back(cur_elm);
	}
}

void ProteinDbParser::ParseCleavageSites(string & _seq, vector<Mutation> & _vAllMutations, size_t & _iAllMutationsSize,
		multimap<size_t, size_t> & _miiMutations, vector<CleavageSite> & _vAllPosibleCleavageSites,
		size_t & _iAllPosibleCleavageSitesSize, vector<size_t> & _vStaticCleavageSites,
		vector<size_t> & _vDynamicCleavageSites) {
	// each position will only have cleavage site
	size_t i = 0;
	CleavageSite * _pCleavageSite;
	multimap<size_t, size_t>::iterator itBefore, itAfter, itLowBefore, itUpBefore, itLowAfter, itUpAfter;
	_iAllPosibleCleavageSitesSize = 0;
	if (_iAllPosibleCleavageSitesSize >= _vAllPosibleCleavageSites.size()) {
		_vAllPosibleCleavageSites.push_back(CleavageSite());
	}
	_pCleavageSite = &(_vAllPosibleCleavageSites.at(_iAllPosibleCleavageSitesSize));
	_pCleavageSite->iPos = 0;
	_pCleavageSite->iMutationPositionType = 0;
	_pCleavageSite->vMutations.clear();
	_pCleavageSite->bIsDynamicCleavageSite = false;
	_vStaticCleavageSites.clear();
	_vStaticCleavageSites.push_back(0);
	_vDynamicCleavageSites.clear();
	bool bDynamicCleavageSite = false;
	for (; i < _seq.length() - 1; i++) {
		if (isCleavageSite(_seq.at(i), _seq.at(i + 1))) {
			++_iAllPosibleCleavageSitesSize;
			if (_iAllPosibleCleavageSitesSize >= _vAllPosibleCleavageSites.size()) {
				_vAllPosibleCleavageSites.push_back(CleavageSite());
			}
			_pCleavageSite = &(_vAllPosibleCleavageSites.at(_iAllPosibleCleavageSitesSize));
			_pCleavageSite->iPos = i + 1;
			_pCleavageSite->iMutationPositionType = 0;
			_pCleavageSite->vMutations.clear();
			// if there are mutations at position i or i + 1
			itLowBefore = _miiMutations.lower_bound(i);
			itUpBefore = _miiMutations.upper_bound(i);
			itLowAfter = _miiMutations.lower_bound(i + 1);
			itUpAfter = _miiMutations.upper_bound(i + 1);
			bDynamicCleavageSite = false;
			// try the combination of mutations before and after cleavage site
			for (itBefore = itLowBefore; itBefore != itUpBefore; ++itBefore) {
				for (itAfter = itLowAfter; itAfter != itUpAfter; ++itAfter) {
					if (!isCleavageSite(_seq, _vAllMutations.at((*itBefore).second),
							_vAllMutations.at((*itAfter).second), i + 1)) {
						_vDynamicCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
						_pCleavageSite->bIsDynamicCleavageSite = true;
						bDynamicCleavageSite = true;
						break;
					}
				}
			}
			if (bDynamicCleavageSite) {
				continue;
			}
			for (itBefore = itLowBefore; itBefore != itUpBefore; ++itBefore) {
				if (!isCleavageSite(_seq, _vAllMutations.at((*itBefore).second), i + 1)) {
					_vDynamicCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
					_pCleavageSite->bIsDynamicCleavageSite = true;
					bDynamicCleavageSite = true;
					break;
				}
			}
			// as long as one mutation make it not a cleavage site, then it is a dynamic cleavage site
			if (bDynamicCleavageSite) {
				continue;
			}
			for (itAfter = itLowAfter; itAfter != itUpAfter; ++itAfter) {
				if (!isCleavageSite(_seq, _vAllMutations.at((*itAfter).second), i + 1)) {
					_vDynamicCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
					_pCleavageSite->bIsDynamicCleavageSite = true;
					bDynamicCleavageSite = true;
					break;
				}
			}
			if (!bDynamicCleavageSite) {
				_vStaticCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
				_pCleavageSite->bIsDynamicCleavageSite = false;
			}
		} else {
			// if it is not a cleavage site on the original sequence
			// if there are mutations at position i or i + 1
			itLowBefore = _miiMutations.lower_bound(i);
			itUpBefore = _miiMutations.upper_bound(i);
			itLowAfter = _miiMutations.lower_bound(i + 1);
			itUpAfter = _miiMutations.upper_bound(i + 1);
			bDynamicCleavageSite = false;
			for (itBefore = itLowBefore; itBefore != itUpBefore; ++itBefore) {
				for (itAfter = itLowAfter; itAfter != itUpAfter; ++itAfter) {
					if (isCleavageSite(_seq, _vAllMutations.at((*itBefore).second),
							_vAllMutations.at((*itAfter).second), i + 1)) {
						++_iAllPosibleCleavageSitesSize;
						if (_iAllPosibleCleavageSitesSize >= _vAllPosibleCleavageSites.size()) {
							_vAllPosibleCleavageSites.push_back(CleavageSite());
						}
						_pCleavageSite = &(_vAllPosibleCleavageSites.at(_iAllPosibleCleavageSitesSize));
						_pCleavageSite->iPos = i + 1;
						_pCleavageSite->iMutationPositionType = 3;
						_pCleavageSite->vMutations.clear();
						_pCleavageSite->vMutations.push_back(&(_vAllMutations.at((*itBefore).second)));
						_pCleavageSite->vMutations.push_back(&(_vAllMutations.at((*itAfter).second)));
						_vDynamicCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
						_pCleavageSite->bIsDynamicCleavageSite = true;
						bDynamicCleavageSite = true;
						break;
					}
				}
			}
			// as long as one mutation make it a cleavage site, then it is a dynamic cleavage site
			if (bDynamicCleavageSite) {
				continue;
			}
			for (itBefore = itLowBefore; itBefore != itUpBefore; ++itBefore) {
				if (isCleavageSite(_seq, _vAllMutations.at((*itBefore).second), i + 1)) {
					++_iAllPosibleCleavageSitesSize;
					if (_iAllPosibleCleavageSitesSize >= _vAllPosibleCleavageSites.size()) {
						_vAllPosibleCleavageSites.push_back(CleavageSite());
					}
					_pCleavageSite = &(_vAllPosibleCleavageSites.at(_iAllPosibleCleavageSitesSize));
					_pCleavageSite->iPos = i + 1;
					_pCleavageSite->iMutationPositionType = 1;
					_pCleavageSite->vMutations.clear();
					_pCleavageSite->vMutations.push_back(&(_vAllMutations.at((*itBefore).second)));
					_vDynamicCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
					_pCleavageSite->bIsDynamicCleavageSite = true;
					bDynamicCleavageSite = true;
					break;
				}
			}
			// as long as one mutation make it a cleavage site, then it is a dynamic cleavage site
			if (bDynamicCleavageSite) {
				continue;
			}
			for (itAfter = itLowAfter; itAfter != itUpAfter; ++itAfter) {
				if (isCleavageSite(_seq, _vAllMutations.at((*itAfter).second), i + 1)) {
					++_iAllPosibleCleavageSitesSize;
					if (_iAllPosibleCleavageSitesSize >= _vAllPosibleCleavageSites.size()) {
						_vAllPosibleCleavageSites.push_back(CleavageSite());
					}
					_pCleavageSite = &(_vAllPosibleCleavageSites.at(_iAllPosibleCleavageSitesSize));
					_pCleavageSite->iPos = i + 1;
					_pCleavageSite->iMutationPositionType = 2;
					_pCleavageSite->vMutations.clear();
					_pCleavageSite->vMutations.push_back(&(_vAllMutations.at((*itAfter).second)));
					_vDynamicCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
					_pCleavageSite->bIsDynamicCleavageSite = true;
					bDynamicCleavageSite = true;
					break;
				}
			}
		}
	}
	++_iAllPosibleCleavageSitesSize;
	if (_iAllPosibleCleavageSitesSize >= _vAllPosibleCleavageSites.size()) {
		_vAllPosibleCleavageSites.push_back(CleavageSite());
	}
	_pCleavageSite = &(_vAllPosibleCleavageSites.at(_iAllPosibleCleavageSitesSize));
	_pCleavageSite->iPos = _seq.length();
	_pCleavageSite->iMutationPositionType = 0;
	_pCleavageSite->vMutations.clear();
	_pCleavageSite->bIsDynamicCleavageSite = false;
	_vStaticCleavageSites.push_back(_iAllPosibleCleavageSitesSize);
	mscCleavageSites.clear();
	for (size_t i = 0; i < _iAllPosibleCleavageSitesSize; ++i) {
		_pCleavageSite = &(_vAllPosibleCleavageSites.at(i));
		mscCleavageSites[_pCleavageSite->iPos] = _pCleavageSite;
	}
	++_iAllPosibleCleavageSitesSize;
}

void ProteinDbParser::ParseMutationInfo(string & _sInfo, vector<Mutation> & _vAllMutations, size_t & _iAllMutationsSize,
		vector<vector<size_t>> & _vviAllMutationLinks, size_t & _iAllMutationLinksSize,
		multimap<size_t, size_t> & _miiMutations) {
	_iAllMutationsSize = 0;
	int iPos = 0;
	if (bSkipM) {
		--iPos;
	}
	sPositionOfMutation.clear();
	sMutationResidue.clear();
	size_t i;
	Mutation * pMutation;
	for (i = 0; i < _sInfo.size(); i++) {
		if (_sInfo.at(i) == ',') {
			// link information
			break;
		}
		if (isdigit(_sInfo.at(i))) {
			sPositionOfMutation.push_back(_sInfo.at(i));
		} else {
			sMutationResidue.push_back(_sInfo.at(i));
			iPos += atoi(sPositionOfMutation.c_str());
			if (iPos < 0) {
				sPositionOfMutation.clear();
				sMutationResidue.clear();
				continue;
			}
			if (_iAllMutationsSize >= _vAllMutations.size()) {
				_vAllMutations.push_back(Mutation());
			}
			pMutation = &(_vAllMutations.at(_iAllMutationsSize));
			pMutation->iPos = iPos;
			pMutation->iMutationIndex = _iAllMutationsSize;
			pMutation->iMutationLinkIndex = EMPTY;
			if (sMutationResidue.at(0) == INSERT) {
				pMutation->cType = INSERT;
			} else if (sMutationResidue.at(0) == DELETE) {
				pMutation->cType = DELETE;
			} else {
				pMutation->cType = SUBSITUE;
				pMutation->sSeq = sMutationResidue;
			}
			sPositionOfMutation.clear();
			sMutationResidue.clear();
			++_iAllMutationsSize;
		}
	}
	// put mutation into a multimap
	_miiMutations.clear();
	for (size_t j = 0; j < _iAllMutationsSize; j++) {
		_miiMutations.insert(pair<size_t, size_t>(_vAllMutations.at(j).iPos, j));
	}
	// parse link information
	_iAllMutationLinksSize = -1;
	sPositionOfMutation.clear();
	sMutationResidue.clear();
	for (; i < _sInfo.size(); i++) {
		if (_sInfo.at(i) == ',') {
			if (!sPositionOfMutation.empty()) {
				_vviAllMutationLinks.at(_iAllMutationLinksSize).push_back(atoi(sPositionOfMutation.c_str()));
				sPositionOfMutation.clear();
			}
			++_iAllMutationLinksSize;
			if (_iAllMutationLinksSize >= _vviAllMutationLinks.size()) {
				_vviAllMutationLinks.push_back(vector<size_t>());
			}
		}
		if (isdigit(_sInfo.at(i))) {
			sPositionOfMutation.push_back(_sInfo.at(i));
		} else {
			if (!sPositionOfMutation.empty()) {
				_vviAllMutationLinks.at(_iAllMutationLinksSize).push_back(atoi(sPositionOfMutation.c_str()));
				sPositionOfMutation.clear();
			}
		}
	}
	if (!sPositionOfMutation.empty()) {
		_vviAllMutationLinks.at(_iAllMutationLinksSize).push_back(atoi(sPositionOfMutation.c_str()));
		sPositionOfMutation.clear();
	}
	++_iAllMutationLinksSize;
	// put the mutation link to all associated mutations
	for (i = 0; i < _iAllMutationLinksSize; i++) {
		for (size_t j = 0; j < _vviAllMutationLinks.at(i).size(); j++) {
			_vAllMutations.at(_vviAllMutationLinks.at(i).at(j)).iMutationLinkIndex = i;
		}
	}
}

void ProteinDbParser::RemoveIllegalResidue(string& seq) {
	seq.erase(seq.find_last_not_of(" \n\r\t") + 1);
	size_t found = seq.find_first_not_of(slegalChar);
	while (found != string::npos) {
		if (bScreenOutput) {
			cout << "Remove illegal character " << seq.substr(found, 1) << " of " << seq << endl;
		}
		seq.erase(found, 1);
		found = seq.find_first_not_of(slegalChar);
	}
}
