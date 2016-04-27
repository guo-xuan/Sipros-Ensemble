/*
 * TableIsotopeDistribution.cpp
 *
 *  Created on: Oct 1, 2015
 *      Author: xgo
 */

#include "TableIsotopeDistribution.h"

TableIsotopeDistribution::TableIsotopeDistribution() {
	this->iNumberOfUniqueResidue = 0;
	this->iNumberForPrecalculatedResidue = 0;
}

TableIsotopeDistribution::~TableIsotopeDistribution() {
	for (uint i = 0; i < this->vpIsotopeDistributionNterm.size(); i++) {
		for (uint j = 0; j < (*this->vpIsotopeDistributionNterm.at(i)).size(); j++) {
			delete (this->vpIsotopeDistributionNterm.at(i))->at(j);
		}
		delete this->vpIsotopeDistributionNterm.at(i);
	}
	for (uint i = 0; i < this->vpIsotopeDistributionCterm.size(); i++) {
		for (uint j = 0; j < (*this->vpIsotopeDistributionCterm.at(i)).size(); j++) {
			delete (this->vpIsotopeDistributionCterm.at(i))->at(j);
		}
		delete this->vpIsotopeDistributionCterm.at(i);
	}
}

bool compareInt(int _i, int _j) {
	return (_i < _j);
}

bool TableIsotopeDistribution::computeProductIon(string _sSequence,
		vector<vector<double> > & _vvdYionMass, vector<vector<double> > & _vvdYionProb,
		vector<vector<double> > & _vvdBionMass, vector<vector<double> > & _vvdBionProb) {
	int iPeptideLength = 0;
	for (uint i = 1; i < _sSequence.length() - 1; i++) {
		if (!isalpha(_sSequence.at(i))) {
			cout << "only for non-PTM search" << endl;
		} else {
			iPeptideLength = iPeptideLength + 1;
		}
	}
	double dProtonMass = ProNovoConfig::getProtonMass();

	if (iPeptideLength < ProNovoConfig::getMinPeptideLength()) {
		cerr << "ERROR: Peptide sequence is too short " << _sSequence << endl;
		return false;
	}

	_vvdYionMass.clear();
	_vvdYionProb.clear();
	_vvdBionMass.clear();
	_vvdBionProb.clear();

	_vvdYionMass.reserve(iPeptideLength);
	_vvdYionProb.reserve(iPeptideLength);
	_vvdBionMass.reserve(iPeptideLength);
	_vvdBionProb.reserve(iPeptideLength);

	if (_sSequence[0] != '[') {
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}
	vector<int> vi;
	int index = 0;
	IsotopeDistribution pIsotopeDistribution1;
	IsotopeDistribution * pIsotopeDistribution2 = NULL;
	IsotopeDistribution sumDistribution;
	// compute B-ion series
	int j = 0;
	for (j = 0; j < iPeptideLength - 1; j++) {
		if (j < this->iNumberForPrecalculatedResidue) {
			vi.push_back(this->mResidueToIndex.find(_sSequence.at(j + 1))->second);
			sort(vi.begin(), vi.end(), compareInt);
			index = this->getIndexByGivenResidues(vi);
			pIsotopeDistribution1 = *(this->vpIsotopeDistributionNterm.at(j)->at(index));
			_vvdBionMass.push_back(pIsotopeDistribution1.vMass);
			_vvdBionProb.push_back(pIsotopeDistribution1.vProb);
		} else {
			pIsotopeDistribution2 =
					&ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find(
							_sSequence.substr(j + 1, 1))->second;
			sum(&pIsotopeDistribution1, pIsotopeDistribution2, &sumDistribution);
			pIsotopeDistribution1 = sumDistribution;
			sumDistribution.clear();
			_vvdBionMass.push_back(pIsotopeDistribution1.vMass);
			_vvdBionProb.push_back(pIsotopeDistribution1.vProb);
		}
	}

	// compute Y-ion series
	vi.clear();
	for (j = 0; j < iPeptideLength - 1; j++) {
		if (j < this->iNumberForPrecalculatedResidue) {
			vi.push_back(this->mResidueToIndex.find(_sSequence.at(iPeptideLength - j))->second);
			sort(vi.begin(), vi.end(), compareInt);
			index = this->getIndexByGivenResidues(vi);
			pIsotopeDistribution1 = *(this->vpIsotopeDistributionCterm.at(j)->at(index));
			_vvdYionMass.push_back(pIsotopeDistribution1.vMass);
			_vvdYionProb.push_back(pIsotopeDistribution1.vProb);
		} else {
			pIsotopeDistribution2 =
					&ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find(
							_sSequence.substr(iPeptideLength - j, 1))->second;
			sum(&pIsotopeDistribution1, pIsotopeDistribution2, &sumDistribution);
			pIsotopeDistribution1 = sumDistribution;
			sumDistribution.clear();
			_vvdYionMass.push_back(pIsotopeDistribution1.vMass);
			_vvdYionProb.push_back(pIsotopeDistribution1.vProb);
		}
	}

	// change the masses to correct for the proton transfer during peptide bond cleavage
	uint n;
	uint m;
	for (n = 0; n < _vvdYionMass.size(); ++n) {
		for (m = 0; m < _vvdYionMass.at(n).size(); ++m) {
			_vvdYionMass.at(n).at(m) += dProtonMass;
		}
	}

	for (n = 0; n < _vvdBionMass.size(); ++n) {
		for (m = 0; m < _vvdBionMass.at(n).size(); ++m) {
			_vvdBionMass.at(n).at(m) -= dProtonMass;
		}
	}
	return true;
}

void TableIsotopeDistribution::debug() {
	vector<int> vSequence;
	int n = 0;
	for (int i = 0; i < this->iNumberForPrecalculatedResidue; i++) {
		uint j = 0;
		for (j = 0; j < this->vpIsotopeDistributionNterm.at(i)->size(); j++) {
			this->decodeIndexToResidues(j, i + 1, this->iNumberOfUniqueResidue, vSequence);
			n = this->getIndexByGivenResidues(vSequence);
			if (n != (int) j) {
				cout << "error, " << i << ", " << j << ", " << n << endl;
			}
		}
	}
}

void calculateCardinalityOfMultiset(const int & _n, const int & _k, int & _result) {
	double result = 1;
	for (int i = 0; i < _k; i++) {
		result *= ((double) _n + i) / ((double) _k - i);
	}
	_result = round(result);
}

void TableIsotopeDistribution::decodeIndexToResidues(const int & _index, const int & _size,
		const int & _uniqueResidue, vector<int> & _vector) {
	//get the most significant digit first
	int iPreviousDigit = 0;
	int iCurrentDigit = 0;
	int i = 0;
	int iTemp = 0;
	int count = _index;
	_vector.clear();
	for (iPreviousDigit = _uniqueResidue - 1; iPreviousDigit >= 0; iPreviousDigit--) {
		iTemp = this->vviIndexDictionary[_size - 1].at(iPreviousDigit);
		if (count >= iTemp) {
			iPreviousDigit++;
			count -= iTemp;
			break;
		}
	}
	if (iPreviousDigit < 0)
		iPreviousDigit = 0;
	_vector.push_back(iPreviousDigit);
	//get the rest digits
	for (i = _size - 2; i >= 0; i--) {
		for (iCurrentDigit = _uniqueResidue - 1; iCurrentDigit >= iPreviousDigit; iCurrentDigit--) {
			if (iCurrentDigit == 0) {
				iTemp = 0;
			} else {
				if (iPreviousDigit == 0) {
					iTemp = this->vviIndexDictionary[i].at(iCurrentDigit - 1);
				} else {
					iTemp = this->vviIndexDictionary[i].at(iCurrentDigit - 1)
							- this->vviIndexDictionary[i].at(iPreviousDigit - 1);
				}
			}
			if (count >= iTemp) {
				count -= iTemp;
				iPreviousDigit = iCurrentDigit;
				_vector.push_back(iPreviousDigit);
				break;
			}
		}
	}
}
/**
 * _positionEnd excluded
 */
IsotopeDistribution* TableIsotopeDistribution::getDistributionByGivenResidues(
		const string & _sSequence, const int & _positionBegin, const int & _positionEnd) {
	int index = this->getIndexByGivenResidues(_sSequence, _positionBegin, _positionEnd);
	return this->vpIsotopeDistributionNterm.at(_positionEnd - _positionBegin)->at(index);
}

int TableIsotopeDistribution::getIndexByGivenResidues(const string & _sSequence,
		const int & _positionBegin, const int & _positionEnd) {
	vector<int> vi;
	map<char, int>::iterator ite;
	for (int i = _positionBegin; i < _positionEnd; i++) {
		ite = this->mResidueToIndex.find(_sSequence.at(i));
		if (ite != this->mResidueToIndex.end()) {
			vi.push_back(ite->second);
		} else {
			cout << "error in from sequence to index" << endl;
			exit(1);
		}
	}
	sort(vi.begin(), vi.end(), compareInt);
	return this->getIndexByGivenResidues(vi);
}

int TableIsotopeDistribution::getIndexByGivenResidues(const vector<int> & _viSequence) {
	int index = 0;
	//process the first digit
	if (_viSequence.at(0) > 0) {
		index += this->vviIndexDictionary.at(_viSequence.size() - 1).at(_viSequence[0] - 1);
	}
	//process the rest digits
	for (uint i = 1; i < _viSequence.size(); i++) {
		if (_viSequence.at(i) > _viSequence.at(i - 1)) {
			index += this->vviIndexDictionary.at(_viSequence.size() - i - 1).at(
					_viSequence.at(i) - 1);
			if (_viSequence.at(i - 1) > 0) {
				index -= this->vviIndexDictionary.at(_viSequence.size() - i - 1).at(
						_viSequence.at(i - 1) - 1);
			}
		}
	}
	return index;
}

/**
 * _positionEnd excluded
 */
int TableIsotopeDistribution::getIndexByGivenResidues(const vector<int> & _viSequence,
		const int & _positionBegin, const int & _positionEnd) {
	int index = 0;
	//process the first digit
	if (_viSequence[_positionBegin] > 0) {
		index += this->vviIndexDictionary.at(_positionEnd - _positionBegin - 1).at(
				_viSequence.at(_positionBegin) - 1);
	}
	//process the rest digits
	for (int i = _positionBegin + 1; i < _positionEnd; i++) {
		if (_viSequence.at(i) > _viSequence.at(i - 1)) {
			index += this->vviIndexDictionary.at(
					_positionEnd - _positionBegin - 1 - (i - _positionBegin)).at(
					_viSequence.at(i) - 1);
			if (_viSequence.at(i - 1) > 0) {
				index -= this->vviIndexDictionary.at(
						_positionEnd - _positionBegin - 1 - (i - _positionBegin)).at(
						_viSequence.at(i - 1) - 1);
			}
		}
	}
	return index;
}

void TableIsotopeDistribution::initializeEverything(const int & _numberOfPrecalculatedResidues) {

	this->setResidueIndex();
	this->initializeIndexDictionary(iNumberOfUniqueResidue, _numberOfPrecalculatedResidues);
	this->setVectorSize(iNumberOfUniqueResidue, _numberOfPrecalculatedResidues);
}

void TableIsotopeDistribution::initializeIndexDictionary(const int & _nResidues,
		const int & _nPrecalculatedResidues) {
	int n = 0;
	for (int i = 0; i < _nPrecalculatedResidues; i++) {
		vector<int> _vCount(_nResidues);
		fill(_vCount.begin(), _vCount.end(), 0);
		for (int j = 0; j < _nResidues; j++) {
			calculateCardinalityOfMultiset(_nResidues - j, i, n);
			_vCount.at(j) = n;
			if (j > 0) {
				_vCount.at(j) += _vCount.at(j - 1);
			}
		}
		this->vviIndexDictionary.push_back(_vCount);
	}
}

IsotopeDistribution * sum(const IsotopeDistribution * distribution0,
		const IsotopeDistribution * distribution1, IsotopeDistribution * sumDistribution) {
	double currentMass;
	double currentProb;
	size_t iSizeDistribution0 = distribution0->vMass.size();
	size_t iSizeDistribution1 = distribution1->vMass.size();
	for (size_t k = 0; k < iSizeDistribution0 + iSizeDistribution1 - 1; k++) {
		double sumweight = 0, summass = 0;
		size_t start = k < (iSizeDistribution1 - 1) ? 0 : k - iSizeDistribution1 + 1; // max(0, k-f_n+1)
		size_t end = k < (iSizeDistribution0 - 1) ? k : iSizeDistribution0 - 1; // min(g_n - 1, k)
		for (size_t i = start; i <= end; i++) {
			double weight = distribution0->vProb[i] * distribution1->vProb[k - i];
			double mass = distribution0->vMass[i] + distribution1->vMass[k - i];
			sumweight += weight;
			summass += weight * mass;
		}
		currentMass = summass / sumweight;
		currentProb = sumweight;
		sumDistribution->vMass.push_back(currentMass);
		sumDistribution->vProb.push_back(currentProb);
	}
	int iSizeSumDistribution;
	int i;

	//prune small probabilities
	vector<double>::iterator iteProb = sumDistribution->vProb.begin();
	vector<double>::iterator iteMass = sumDistribution->vMass.begin();
	while (iteProb != sumDistribution->vProb.end()) {
		if ((*iteProb) > ProNovoConfig::configIsotopologue.ProbabilityCutoff)
			break;
		iteProb++;
		iteMass++;
	}
	sumDistribution->vProb.erase(sumDistribution->vProb.begin(), iteProb);
	sumDistribution->vMass.erase(sumDistribution->vMass.begin(), iteMass);

	while (1) {
		if (sumDistribution->vProb.size() == 0)
			break;
		if (sumDistribution->vProb.back() > ProNovoConfig::configIsotopologue.ProbabilityCutoff)
			break;
		sumDistribution->vProb.pop_back();
		sumDistribution->vMass.pop_back();
	}

	// normalize the probability space to 1
	double sumProb = 0;
	iSizeSumDistribution = sumDistribution->vMass.size();
	for (i = 0; i < iSizeSumDistribution; ++i)
		sumProb += sumDistribution->vProb[i];

	if (sumProb <= 0)
		return sumDistribution;

	for (i = 0; i < iSizeSumDistribution; ++i)
		sumDistribution->vProb[i] = sumDistribution->vProb[i] / sumProb;
	return sumDistribution;

}

void TableIsotopeDistribution::precalculateIsotopicDistribution(const vector<int> & _viSequence,
		IsotopeDistribution * _sumDistribution) {
	if (_viSequence.size() == 1) {
		(*_sumDistribution) = ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find(
				this->mIndexToResidue.find(_viSequence.at(0))->second)->second;
	} else {
		int index = this->getIndexByGivenResidues(_viSequence, 1, _viSequence.size());
		IsotopeDistribution * distribution1 = (this->vpIsotopeDistributionNterm.at(
				_viSequence.size() - 2))->at(index);
		IsotopeDistribution * distribution2 =
				&ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find(
						this->mIndexToResidue.find(_viSequence.at(0))->second)->second;
		sum(distribution1, distribution2, _sumDistribution);
	}
}

bool compareVector(const vector<int> & _v1, const vector<int> & _v2) {
	if (_v1.size() != _v2.size()) {
		return false;
	}
	for (uint i = 0; i < _v1.size(); i++) {
		if (_v1.at(i) != _v2.at(i)) {
			return false;
		}
	}
	return true;
}

bool TableIsotopeDistribution::setResidueIndex() {

	this->mResidueToIndex.clear();
	this->mIndexToResidue.clear();
	map<string, vector<int> >::iterator ite1;
	map<char, int>::iterator ite2;
	map<string, vector<int> > mResidueAtomicComposition =
			ProNovoConfig::configIsotopologue.mResidueAtomicComposition;
	int index = 0;
	bool bFind = false;
	for (ite1 = mResidueAtomicComposition.begin(); ite1 != mResidueAtomicComposition.end();
			ite1++) {
		if (isalpha(ite1->first[0]) && (ite1->first.size() == 1)) {
			bFind = false;
			for (ite2 = this->mResidueToIndex.begin(); ite2 != this->mResidueToIndex.end();
					ite2++) {
				if (compareVector(ite1->second,
						(mResidueAtomicComposition.find(string(1, ite2->first)))->second)) {
					this->mResidueToIndex[ite1->first.at(0)] = ite2->second;
					bFind = true;
					break;
				}
			}
			if (bFind == true) {
				continue;
			}
			this->mResidueToIndex[ite1->first.at(0)] = index;
			this->mIndexToResidue[index] = ite1->first;
			index++;
		}

	}
	this->iNumberOfUniqueResidue = index;
	return true;
}

void TableIsotopeDistribution::setVectorSize(const int & _nUniqueResidues,
		const int & _nPrecalculatedResidues) {
	this->iNumberForPrecalculatedResidue = _nPrecalculatedResidues;
	this->vpIsotopeDistributionNterm.resize(_nPrecalculatedResidues);
	this->vpIsotopeDistributionCterm.resize(_nPrecalculatedResidues);
	int count = 0;
	for (int i = 0; i < _nPrecalculatedResidues; i++) {
		calculateCardinalityOfMultiset(_nUniqueResidues, i + 1, count);
		vector<IsotopeDistribution*> * vpNterm = new vector<IsotopeDistribution*>(count);
		this->vpIsotopeDistributionNterm.at(i) = vpNterm;
		vector<IsotopeDistribution*> * vpCterm = new vector<IsotopeDistribution*>(count);
		this->vpIsotopeDistributionCterm.at(i) = vpCterm;
	}
}

void TableIsotopeDistribution::startPrecalculation() {
	vector<int> vSequence;
	for (int i = 0; i < this->iNumberForPrecalculatedResidue; i++) {
		int size = this->vpIsotopeDistributionNterm.at(i)->size();
		int j = 0;
#pragma omp parallel for \
	private(vSequence, j) \
		schedule(guided)
		for (j = 0; j < size; j++) {
			this->decodeIndexToResidues(j, i + 1, this->iNumberOfUniqueResidue, vSequence);
			IsotopeDistribution* sumDistribution = new IsotopeDistribution();
			this->precalculateIsotopicDistribution(vSequence, sumDistribution);
			this->vpIsotopeDistributionNterm.at(i)->at(j) = sumDistribution;
		}

	}

	//add Cterm
	IsotopeDistribution * cterm =
			&ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find("Cterm")->second;
	for (int i = 0; i < this->iNumberForPrecalculatedResidue; i++) {
		int size = this->vpIsotopeDistributionNterm.at(i)->size();
		int j = 0;
#pragma omp parallel for \
	private(j) \
		schedule(guided)
		for (j = 0; j < size; j++) {
			IsotopeDistribution* sumDistribution = new IsotopeDistribution();
			IsotopeDistribution * distribution1 = (this->vpIsotopeDistributionNterm.at(i))->at(j);
			sum(distribution1, cterm, sumDistribution);
			sumDistribution->swap();
			this->vpIsotopeDistributionCterm.at(i)->at(j) = sumDistribution;
		}

	}

	//add Nterm
	IsotopeDistribution * nterm =
			&ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find("Nterm")->second;
	IsotopeDistribution sumDistribution;
	for (int i = 0; i < this->iNumberForPrecalculatedResidue; i++) {
		int size = this->vpIsotopeDistributionNterm.at(i)->size();
		int j = 0;
#pragma omp parallel for \
		private(j, sumDistribution) \
			schedule(guided)
		for (j = 0; j < size; j++) {
			IsotopeDistribution * distribution1 = (this->vpIsotopeDistributionNterm.at(i))->at(j);
			sumDistribution.vMass.clear();
			sumDistribution.vProb.clear();
			sum(distribution1, nterm, &sumDistribution);
			sumDistribution.swap();
			*(this->vpIsotopeDistributionNterm.at(i)->at(j)) = sumDistribution;
		}

	}
}
