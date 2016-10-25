/*
 * ProteinDbParser.h
 *
 *  Created on: Sep 29, 2016
 *      Author: xgo
 */

#ifndef PROTEINDBPARSER_H_
#define PROTEINDBPARSER_H_

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>

#include "peptide.h"
#include "proNovoConfig.h"
#include "ptm.h"

#define EMPTY	1000000
#define INSERT	'^'
#define DELETE	'`'
#define SUBSITUE	'S'

struct Mutation {
	size_t iPos;
	char cType;
	string sSeq;
	size_t iMutationLinkIndex;
	size_t iMutationIndex;
};

struct CleavageSite {
	/**
	 * iPos Position of the cleavage sites, cut between KX, the position of X
	 * iVariantIndex	if it is EMPTY, the cleavage site is from original sequence.
	 * 					if it is not EMPTY, the cleavage site is from a mutation.
	 * vMutations 	if size is EMPTY, the cleavage site is a static cleavage site from original sequence.
	 * 				if it is not EMPTY, the cleavage site is from some mutations.
	 */
	size_t iPos;
	size_t iMutationIndex;
	bool bIsDynamicCleavageSite;
	vector<Mutation *> vMutations;
	// 0: no mutation; 1: before cleavage site; 2 after cleavage site; 3: both
	size_t iMutationPositionType;
	bool operator<(const CleavageSite &cleavage) const {
		return iPos < cleavage.iPos;
	}
};

struct ProteinSegment {
	// the original sequence without any variants
	string sSeq;
	// Position on the original sequence, first included, second not included
	pair<size_t, size_t> pPostion;
	vector<size_t> vPositiveDynamicCleavageSites;
	char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
};

struct MutatedPeptide {
	string sSeq;
	ProteinSegment * proteinSegment;
	pair<size_t, size_t> pPostion;
	char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
	double dcurrentMass;
};

class LookupTable {

private:
	vector<bool> vb;
	vector<size_t> vSet;
public:
	void clean();
	bool exist(size_t _index);
	inline const vector<size_t> & get() {
		return vSet;
	}
	bool push(size_t _index);
};

class ProteinDbParser {
public:
	ProteinDbParser(bool _bScreenOutput);
	~ProteinDbParser();

	bool getNextProtein();
	size_t getNumMutatedPeptide();
	size_t getNumPtmPeptide();
	const vector<MutatedPeptide> & peptideModify(size_t _iMutatedPeptideIndex);
	const vector<MutatedPeptide> & proteinDigest();
	void setPeptide(size_t _index, Peptide * _pPeptide);

private:
	// all protein characters
	string orderstring;
	// from configure file
	PTM_List ptmlist;
	// re-organized based on ptmlist
	vector<vector<pair<string, double> > > ptm_map;
	// position, and vector of PTMs with symbol and mass shift
	// all potential ptm position for the current peptide
	vector<pair<size_t, vector<pair<string, double> > *> > ptm_position_all;
	// based on given comb_order
	vector<int> comb_order, ptm_order, ele_num;
	// the maximum number of PTMs allowed in a peptide, set this variable from ProNovoConfig::getMaxPTMcount()
	size_t imaxPTM;
	// if true, allows standard output
	bool bScreenOutput;
	// the filename of protein DB, set this variable from ProNovoConfig::getFASTAfilename();
	string sDBFilename;
	// the file operator of protein DB
	ifstream db_stream;
	// read a line from protein DB
	string sLine;
	//current protein id
	int iProteinId;
	// the current protein name
	string scurrentProteinName;
	// next protein name; if no more protein, then "" empty
	string snextProteinName;
	// original protein sequence from the protein DB
	string sProteinSeq;
	// if mutation info exists, it is surrounded by sOpenCharMutation and sCloseCharMutation
	static string sOpenCharMutation;
	static string sEndCharMutation;
	// all mutations
	vector<Mutation> vAllMutations;
	// the size of vAllMutations;
	size_t iAllMutationsSize;
	// positions of mutations, index of mutations, only those mutations won't affect cleavage sits
	multimap<size_t, size_t> miiMutations;
	// all mutation links
	vector<vector<size_t>> vviAllMutationLinks;
	// the size of mutation links
	size_t iAllMutationLinksSize;
	// all cleavage sites on the protein sequence, each position will only have cleavage site
	vector<CleavageSite> vAllPosibleCleavageSites;
	// poistion of cleavage sites
	map<size_t, CleavageSite*> mscCleavageSites;
	// the size of vAllPosibleCleavageSites
	size_t iAllPosibleCleavageSitesSize;
	// indexes of cleavage sites not affected by mutations
	vector<size_t> vStaticCleavageSites;
	// indexes cleavage sites affected by mutations
	vector<size_t> vDynamicCleavageSites;
	// Residues identified before the cleavage site
	string sResiduesBeforeCleavageSite;
	// Residues identified after the cleavage site
	string sResiduesAfterCleavageSite;
	// all the mutated peptides in this proteins
	vector<MutatedPeptide> vMutatedPeptides;
	// the size of vMutatedPeptides
	size_t iMutatedPeptidesSize;
	// cleavage sites to do the protein digestions, not all cleavage site will be used
	vector<CleavageSite*> vCombinationCleavageSites;
	// dynamic cleavage sites must be covered by a peptide
	vector<size_t> vSelectedDynamicCleavageSites;
	// the protein segments after digestion
	vector<ProteinSegment> vProteinSegments;
	// the size of vProteinSegments
	size_t iProteinSegmentsSize;
	// all the modified peptide for a mutated peptide
	vector<MutatedPeptide> vPtmPeptides;
	// size of vPtmPeptides
	size_t iPtmPeptidesSize;
	// the link mask used to make sure there is no link conflict
	vector<int64_t> vLinkMasks;
	// the position mask used to make sure there is no position conflict
	vector<int64_t> vPositionMasks;
	// mutations in the range of a peptide
	vector<size_t> vMutationsInPeptide;
	// all possible residues
	string slegalChar;
	// temporary storing the mutation position info
	string sPositionOfMutation;
	// temporary storing the mutation string
	string sMutationResidue;
	// temporary storing the mutation info from the protein ID
	string sMutationInfo;
	// temporary storing the begin position of mutation info string in the protein ID sequence
	size_t iMutationInfoOpenPos;
	// temporary storing the end position of mutation info string in the protein ID sequence
	size_t iMutationInfoClosePos;
	// maximum allowed cleavage sites, get it from configure file
	int iMaxMissedCleavageSites;
	// quick check if two mutations are in the same position
	LookupTable lookupTablePosition;
	// N term mass
	double dTerminusMassN;
	// C term mass
	double dTerminusMassC;
	// mass for each possible residue
	map<char, double> mResiduleMass;
	// minimum peptide length
	size_t iMinPeptideLen;
	// maximum peptide length
	size_t iMaxPeptideLen;
	// bool skip the first M or not
	bool bSkipM;

	void applyMutations(ProteinSegment & _proteinsegment, vector<Mutation> & _vAllVariants,
			vector<size_t> & _vIndexDynamicMutations, int64_t _iMasks, vector<MutatedPeptide> & _vMutatedPeptides,
			size_t & _iMutatedPeptidesSize, LookupTable & _lookupTablePosition);
	void digest(string & _seq, vector<CleavageSite> & _vAllPosibleCleavageSites, size_t _iAllPosibleCleavageSitesSize, vector<size_t> & _vStaticCleavageSites,
			vector<size_t> & vDynamicCleavageSites, vector<CleavageSite*> & vCombinationCleavageSites,
			vector<size_t> & vSelectedDynamicCleavageSites, vector<ProteinSegment> & _vProteinSegments,
			size_t & _iProteinSegmentsSize);
	void generateLinkMasks(vector<vector<size_t>> & _vviAllMutationLinks, const vector<size_t> & _vLinkIndexes,
			vector<Mutation> & _vAllMutations, vector<size_t> & _vIndexDynamicMutations, vector<int64_t> & _vMasks);
	void generateMutatedPeptides(ProteinSegment & _proteinsegment, vector<MutatedPeptide> & _vMutatedPeptides,
			size_t & iMutatedPeptidesSize, multimap<size_t, size_t> & _miiMutations, vector<Mutation> & _vAllMutations,
			vector<vector<size_t>> & _vviAllMutationLinks, LookupTable & _lookupTablePosition,
			vector<int64_t> & _vLinkMasks, vector<int64_t> & _vPositionMasks, vector<size_t> & _vDynamicMutations);
	void generateProteinSegment(string & _seq, vector<CleavageSite*> & _vCleavageSites,
			vector<size_t> & _vCleavageSitesMustHave, vector<ProteinSegment> & _vProteinSegments,
			size_t & _iProteinSegmentsSize);
	void generatePtmPeptide(MutatedPeptide & _mutatedPeptide, vector<MutatedPeptide> & _vPtmPeptides,
			size_t & _iPtmPeptidesSize, vector<pair<size_t, vector<pair<string, double> > *> > & _ptm_position_all,
			vector<int> & _comb_order, vector<int> & _ptm_order, vector<int> & _ele_num);
	void generatePositionMasks(vector<Mutation> & _vAllMutations, vector<size_t> & _vIndexDynamicMutations,
			vector<int64_t> & _vMasks);
	bool isCleavageResidueAfter(char c2);
	bool isCleavageResidueBefore(char c1);
	bool isCleavageSite(char c1, char c2);
	bool isCleavageSite(string & _seq, Mutation & _mutation, size_t _iCleavagePos);
	bool isCleavageSite(string & _seq, Mutation & _mutationBefore, Mutation & _mutationAfter, size_t _iCleavagePos);
	void Initial_PTM_Map();
	void ParseCleavageSites(string & _seq, vector<Mutation> & _vAllMutations, size_t & _iAllMutationsSize,
			multimap<size_t, size_t> & _miiMutations, vector<CleavageSite> & _vAllPosibleCleavageSites,
			size_t & _iAllPosibleCleavageSitesSize, vector<size_t> & _vStaticCleavageSites,
			vector<size_t> & _vDynamicCleavageSites);
	void ParseMutationInfo(string & _sInfo, vector<Mutation> & _vAllMutations, size_t & _iAllMutationsSize,
			vector<vector<size_t>> & _vviAllMutationLinks, size_t & _iAllMutationLinksSize,
			multimap<size_t, size_t> & _miiMutations);
	void RemoveIllegalResidue(string & seq);
};

#endif /* PROTEINDBPARSER_H_ */
