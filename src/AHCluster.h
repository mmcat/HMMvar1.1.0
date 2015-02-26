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
#include <map>

#include "MultiAlign.h"

using namespace std;

struct node{
	int left;
	int right;
	double dist;
};



class AHCluster {
public:

	virtual ~AHCluster();

	AHCluster(MultiAlign * align);

	vector<node*> result;
	int npoints;
	int alignlen;
	string targetS;
	int targetG;
	int groups;
	static const double A = 0.75;

	map<string,vector<string> > seqs;
	map<string,vector<string> > seqs_bk;
	vector<int> clusterid;
	vector<string> heads;
	vector<string> heads_bk;
	vector<vector<double> > distmatrix;
	MultiAlign* align;

	void init();
	void cleanData();
	int selectRange(int pos, int* start,int* end);
	double alignSeqSim(string seq1,string seq2,int* n1,int *n2);
	void runAHC();
	void getOptimized();
	void getDistMat();
	double getEntropy(int i);
	double getEntropyDiff(int s1,int s2);
	double find_closest_pair(int* lp, int* rp);
	void cuttree();
	bool validCol(int pos);
	int printCluster(int gnum,string filename);
	map<char,int> getBKcount(int col);
	map<char,int> getMergeCount(int s1,int s2,int col);
//	int getUniqPermutate(string s);


};

#endif /* AHCLUSTER_H_ */
