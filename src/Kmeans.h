/*
 * Kmeans.h
 *
 *  Created on: Jan 28, 2014
 *      Author: mingmingliu
 */

#ifndef KMEANS_H_
#define KMEANS_H_

#include <vector>
#include <map>
#include <string>
#include "MultiAlign.h"

using namespace std;

class Kmeans{
public:
	int k;
	int targetG;
	int len;
	string targetS;

	map<string,string> seqs;
	vector<string> centers;
	vector<string> heads;
	vector<int> labels;
	MultiAlign* align;
	vector<vector<int> > groups;

	Kmeans(int k, MultiAlign * align);
	int editDist(string seq1, string seq2);
	void init();
	void updateCenter(int gnum);
	void printCluter(int gnum,string filename);
	void runKmeans();

};



#endif /* KMEANS_H_ */
