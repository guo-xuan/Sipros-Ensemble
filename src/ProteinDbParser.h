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
	vector<size_t> * vLinkedVariants;
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
	vector<size_t> vMutations;
	bool operator<(const CleavageSite &cleavage) const {
		return iPos < cleavage.iPos;
	}
};

struct ProteinSegment {
	// the original sequence without any variants
	string sSeq;
	// Position on the original sequence, first included, second not included
	pair<int, int> pPostion;
	vector<size_t> vStaticVariants;
};

struct MutatedPeptide {
	string sSeq;
	pair<int, int> pPostion;
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
	const vector<MutatedPeptide> & peptideDigest(size_t _iMutatedPeptideIndex);
	const vector<MutatedPeptide> & proteinDigest();

private:
	//
	string orderstring;
	// from configure file
	PTM_List ptmlist;
	// re-organized based on ptmlist
	vector<vector<pair<string, double> > > ptm_map;
	// position, and vector of PTMs with symbol and mass shift
	// all potential ptm position for the current peptide
	vector<pair<size_t, vector<pair<string, double> > > > ptm_position_all;
	// based on given comb_order
	vector<size_t> comb_order, ptm_order, ele_num;
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

	string scurrentProteinName;

	string snextProteinName;

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
	//
	vector<CleavageSite*> vCombinationCleavageSites;
	//
	vector<size_t> vSelectedDynamicCleavageSites;
	// the protein segments after digestion
	vector<ProteinSegment> vProteinSegments;
	// the size of vProteinSegments
	size_t iProteinSegmentsSize;
	// all the modified peptide for a mutated peptide
	vector<MutatedPeptide> vPtmPeptides;
	// size of vPtmPeptides
	size_t iPtmPeptidesSize;

	string slegalChar;
	string sPositionOfMutation;
	string sMutationResidue;
	string sMutationInfo;
	size_t iMutationInfoOpenPos;
	size_t iMutationInfoClosePos;

	int iMaxMissedCleavageSites;

	LookupTable lookupTableIndex;
	LookupTable lookupTablePosition;

	void applyMutations(ProteinSegment & _proteinsegment, vector<Mutation> & _vAllVariants,
			vector<size_t> & _vIndexStaticMutations, vector<size_t> & _vIndexDynamicMutations, int64_t _iMasks,
			vector<MutatedPeptide> & _vMutatedPeptides);
	void digest(string & _seq, vector<CleavageSite> & _vAllPosibleCleavageSites, vector<size_t> & _vStaticCleavageSites,
			vector<size_t> & vDynamicCleavageSites, vector<CleavageSite*> & vCombinationCleavageSites,
			vector<size_t> & vSelectedDynamicCleavageSites, vector<ProteinSegment> & _vProteinSegments,
			size_t & _iProteinSegmentsSize);
	void generateLinkMasks(vector<vector<size_t>> & _vviAllMutationLinks, vector<size_t> & _vLinkIndexes,
			vector<Mutation> & _vAllVariants, vector<size_t> & _vIndexDynamicMutations, vector<int64_t> & _vMasks);
	void generateMutatedPeptides(ProteinSegment & _proteinsegment, vector<MutatedPeptide> & _vMutatedPeptides,
			size_t & iMutatedPeptidesSize, multimap<size_t, size_t> & _miiMutations, vector<Mutation> & _vAllVariants,
			vector<vector<size_t>> & _vviAllMutationLinks, LookupTable & _lookupTableIndex,
			LookupTable & _lookupTablePosition);
	void generateProteinSegment(string & _seq, vector<CleavageSite*> & _vCleavageSites,
			vector<size_t> & _vCleavageSitesMustHave, vector<ProteinSegment> & _vProteinSegments, size_t & _iProteinSegmentsSize);
	void generatePtmPeptide(MutatedPeptide & _mutatedPeptide, vector<MutatedPeptide> & _vPtmPeptides,
			size_t & _iPtmPeptidesSize, vector<pair<size_t, vector<pair<string, double> > > > & _ptm_position_all,
			vector<size_t> & _comb_order, vector<size_t> & _ptm_order, vector<size_t> & _ele_num);
	void generatePositionMasks(vector<Mutation> & _vAllVariants, vector<size_t> & _vIndexDynamicMutations,
			vector<int64_t> & _vMasks);
	bool isCleavageSite(char c1, char c2);
	bool isCleavageSite(string & _seq, Mutation & _mutation, size_t _iCleavagePos);
	bool isCleavageSite(string & _seq, Mutation & _mutationBefore, Mutation & _mutationAfter, size_t _iCleavagePos);
	bool isConflicting(vector<size_t> & _vStaticVariants);
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
