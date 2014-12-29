/*
 * AHCluster.h
 *
 *  Created on: Dec 28, 2014
 *      Author: mingming
 */

#ifndef AHCLUSTER_H_
#define AHCLUSTER_H_

#include <vector>
#include <string>

#include "MultiAlign.h"

using namespace std;

class AHCluster {
public:
	AHCluster();
	virtual ~AHCluster();

	vector<vector<string> > rep;
	map<string,string> seqs;
	MultiAlign* align;

	void init();

};

#endif /* AHCLUSTER_H_ */
