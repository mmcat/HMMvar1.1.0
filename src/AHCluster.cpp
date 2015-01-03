/*
 * AHCluster.cpp
 *
 *  Created on: Dec 28, 2014
 *      Author: mingming
 */

#include "AHCluster.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <limits.h>

using namespace std;

AHCluster::~AHCluster() {
	// TODO Auto-generated destructor stub
}

AHCluster::AHCluster(MultiAlign *align){
	this->align = align;
	npoints = 0;
	alignlen=0;
}

void AHCluster::init(){
	ifstream fin(align->align_output_file_stockholm_.c_str());
	if(fin.is_open()){
		string line;
		getline(fin,line);
		getline(fin,line);
		while(getline(fin,line)){
			string::size_type p = line.find(" ");
			if (p!=std::string::npos){
				string def = line.substr(0,p);
				string seqseg = line.substr(p+1);
				if(seqs.find(def)==seqs.end()){
					vector<string> tmp(1,seqseg);
					seqs[def]=tmp;
					heads.push_back(def);
				}
				else{
					seqs[def][0]+=seqseg;
				}
			}
		}
		npoints=heads.size();
		alignlen=seqs[heads[0]][0].length();

	}
	else{
		fprintf(stderr,"Cannot open the alignment file!\n");
		exit(-1);
	}
	fin.close();

}

void AHCluster::runAHC(){


	getDistMat();
	if(distmatrix.size()==0) return;
	int nnodes = npoints-1;
	for(int i=0;i<nnodes;i++){
		node* newnode = new node;
		result.push_back(newnode);
	}

	vector<int> distid(npoints,0);
	map<string,vector<string> > newseqs = seqs;
	vector<string> newheads = heads;

	for (int i = 0; i < npoints; i++) distid[i] = i;
	for(int i=0;i<nnodes;i++){
		int is = 1;
		int js = 0;

		result[i]->dist = find_closest_pair(&is, &js);
		result[i]->left = distid[js];
		result[i]->right = distid[is];

		cout<<result[i]->left<<","<<result[i]->right<<":"<<result[i]->dist<<endl;
		cout<<heads[js]<<" & "<<heads[is]<<endl;

		/* Make node js the new node */
		for(int j=0;j<seqs[heads[is]].size();j++)
			seqs[heads[js]].push_back(seqs[heads[is]][j]);

		seqs.erase(heads[is]);
		heads.erase(heads.begin()+is);
		distmatrix.erase(distmatrix.begin()+is);

		for(int j=is;j<distmatrix.size();j++)
			distmatrix[j].erase(distmatrix[j].begin()+is);

		 /* Fix the distances */
		distid.erase(distid.begin()+is);
		distid[js] = -i-1;

		for(int j=0;j<js;j++){
			distmatrix[js][j]=getEntropyDiff(js,j);
		}
		for(int j=js+1;j<distmatrix.size();j++){
			distmatrix[j][js]=getEntropyDiff(j,js);
		}
	}

//	cout<<heads.size()<<" "<<seqs.size()<<" "<<seqs[heads[0]].size()<<endl;
	for(int i=0;i<result.size();i++){
		cout<<-(i+1)<<"\t"<<result[i]->left<<"\t"<<result[i]->right<<"\t"<<result[i]->dist<<endl;
	}


	seqs = newseqs;
	heads = newheads;

}

double AHCluster::find_closest_pair(int* lp, int* rp){
	double d = INT_MAX;
	for(int i=0;i<distmatrix.size();i++){
		for(int j=0;j<i;j++){
			if(distmatrix[i][j]<d){
				d = distmatrix[i][j];
				*lp = i;
				*rp = j;
			}
		}
	}

	return d;

}



double AHCluster::getEntropyDiff(int s1,int s2){
	double s=0;

	int n1 = seqs[heads[s1]].size();
	int n2 = seqs[heads[s2]].size();
//	cout<<heads[s1]<<" "<<heads[s2]<<endl;
	for(int i=0;i<alignlen;i++){
		map<char,int> aacount_bk = getBKcount(i);
		map<char,int> aacount_m = getMergeCount(s1,s2,i);
		/*calculate observed entropy*/
		double obj = 0;
		for(map<char,int>::iterator it = aacount_m.begin();it!=aacount_m.end();it++){
			cout<<it->first<<" "<<it->second<<endl;
			double c = 0;
			for(int j=1;j<=it->second;j++) c+=log(j);
			obj+=c;
		}
		/*calculate expected entropy*/
		double exp = 0;
		for(map<char,int>::iterator it = aacount_bk.begin();it!=aacount_bk.end();it++){
			cout<<it->first<<" "<<it->second<<" "<<(double)it->second*(n1+n2)/npoints<<endl;

			double f=lgamma((double)it->second*(n1+n2)/npoints+1);
			exp+=f;
		}
		s+=exp-obj;
	}

	s=s/alignlen;
	double sp=0;
	for(int i=1;i<=n1+n2;i++)
		sp+=log(i);

	return A*s+(1-A)*sp;

}



void AHCluster::cuttree(){
	vector<int> clusterid_tmp(npoints,-1);
	int cnum=0;
	int i, j, k;
	double opt = INT_MAX;
	for(int nclusters=2;nclusters<=2;nclusters++){
		int icluster = 0;
		const int n = npoints-nclusters; /* number of nodes to join */
		vector<int> nodeid(n,-1);

		for (i = npoints-2; i >= n; i--)
		{ k = result[i]->left;
	    	if (k>=0)
	    	{ clusterid_tmp[k] = icluster;
	    	icluster++;
	    	}
	    	k = result[i]->right;
	    	if (k>=0)
	    	{ clusterid_tmp[k] = icluster;
	    	icluster++;
	    	}
		}

		double s=0;

		for (i = n-1; i >= 0; i--)
		{ if(nodeid[i]<0)
	    	{ j = icluster;
	    	nodeid[i] = j;
	    	s+=result[i]->dist;
	    	icluster++;
	    	}
	    	else j = nodeid[i];
	    	k = result[i]->left;
	    	if (k<0) nodeid[-k-1] = j;else clusterid_tmp[k] = j;
	    	k = result[i]->right;
	    	if (k<0) nodeid[-k-1] = j;else clusterid_tmp[k] = j;
		}

		if(s<opt){
			opt = s;
			cnum = nclusters;
			clusterid = clusterid_tmp;
		}

	}
	for(int i=0;i<heads.size();i++){
		cout<<i<<"\t"<<heads[i]<<"\t"<<clusterid[i]<<endl;
	}
	return;

}

void AHCluster::printCluster(int gnum,string filename){

}

void AHCluster::getDistMat(){
	for(int i=0;i<heads.size();i++){
		vector<double> dist;
		for(int j=0;j<i;j++){
			double d = getEntropyDiff(i,j);
			dist.push_back(d);
		}
		distmatrix.push_back(dist);
	}

}



map<char,int> AHCluster::getBKcount(int col){
	map<char,int> aacount;
	for(auto it=seqs.begin();it!=seqs.end();it++){
		vector<string> cluster = it->second;
		for(int i=0;i<cluster.size();i++){
			aacount[cluster[i][col]]++;
		}
	}
	return aacount;
}

map<char,int> AHCluster::getMergeCount(int s1,int s2,int col){
	map<char,int> aacount;
	for(int i=0;i<seqs[heads[s1]].size();i++)
		aacount[seqs[heads[s1]][i][col]]++;
	for(int i=0;i<seqs[heads[s2]].size();i++)
			aacount[seqs[heads[s2]][i][col]]++;
	return aacount;
}

