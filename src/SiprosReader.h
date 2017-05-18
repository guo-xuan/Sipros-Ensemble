/*
 * SiprosReader.h
 *
 *  Created on: May 17, 2017
 *      Author: xgo
 */

#ifndef SIPROSREADER_H_
#define SIPROSREADER_H_



class SiprosReader{
public:
	bool static MzmlReader(std::string & filename_str, std::vector<Spectrum> *  _vSpectra);
};



#endif /* SIPROSREADER_H_ */
