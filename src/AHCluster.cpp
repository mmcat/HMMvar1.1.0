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
					if(targetS.find(def)!=std::string::npos) targetS = def;

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
	seqs_bk = seqs;
	heads_bk = heads;

}

int AHCluster::selectRange(int pos,int* start,int* end){

	seqs = seqs_bk;
	heads = heads_bk;
	npoints = heads.size();
	alignlen = seqs[heads[0]][0].length();

	int align_pos = 0,charn=0;

	while(charn<pos){
	    if(seqs[targetS][0][align_pos++]!='-') charn++;
	}


	*start=align_pos-1;
	*end=align_pos-1;
	while(*start>=0 && seqs[targetS][0][*start]!='-') (*start)--;
	(*start)++;
	while(*end<alignlen && seqs[targetS][0][*end]!='-') (*end)++;
	(*end)--;

	alignlen = *end-*start+1;
	for(int i=0;i<heads.size();){
		int gapn=0;
		for(int j=*start;j<=*end;j++)
			if(seqs[heads[i]][0][j]=='-') gapn++;
		if(gapn==alignlen) {
			seqs.erase(heads[i]);
			heads.erase(heads.begin()+i);
		}
		else{
			seqs[heads[i]][0] = seqs[heads[i]][0].substr(*start,alignlen);
			i++;
		}
	}

	npoints = heads.size();

	ofstream fout("/Users/mingmingliu/Documents/study/2014spring/class_variants/data/tmp1.sto.aln");
		fout<<"# STOCKHOLM 1.0\n\n";
		for(int i=0;i<heads.size();i++)
			fout<<heads[i]<<" "<<seqs[heads[i]][0]<<endl;
		fout<<"//";
		fout.close();


	return align_pos;

}

void AHCluster::runAHC(){


	getDistMat();
	result.clear();
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

//		cout<<result[i]->left<<","<<result[i]->right<<":"<<result[i]->dist<<endl;
//		cout<<heads[js]<<" & "<<heads[is]<<endl;

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
//			cout<<it->first<<" "<<it->second<<endl;
			double c = 0;
			for(int j=1;j<=it->second;j++) c+=log(j);
			obj+=c;
		}
		/*calculate expected entropy*/
		double exp = 0;
		for(map<char,int>::iterator it = aacount_bk.begin();it!=aacount_bk.end();it++){
//			cout<<it->first<<" "<<it->second<<" "<<(double)it->second*(n1+n2)/npoints<<endl;

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
	for(int nclusters=1;nclusters<=npoints;nclusters++){
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
	    	if (k<0) nodeid[-k-1] = j; else clusterid_tmp[k] = j;
          	k = result[i]->right;
	    	if (k<0) nodeid[-k-1] = j;else clusterid_tmp[k] = j;

		}

		if(s<opt){
			opt = s;
			cnum = nclusters;
			clusterid = clusterid_tmp;
		}

	}
	groups = cnum;
	for(int i=0;i<heads.size();i++){
		if(targetS==heads[i]) targetG = clusterid[i];
		cout<<i<<"\t"<<heads[i]<<"\t"<<clusterid[i]<<endl;
	}
	return;

}

int AHCluster::printCluster(int gnum,string filename){
	ofstream fout(filename.c_str());
	fout<<"# STOCKHOLM 1.0\n\n";
	if(gnum==-1){
		for(int i=0;i<clusterid.size();i++){
			fout<<clusterid[i]<<"_"<<heads[i]<<" "<<seqs[heads[i]][0]<<endl;
		}
		fout<<"//";
		fout.close();

	}
	else{
		int c=0;
		for(int i=0;i<clusterid.size();i++){

			if(clusterid[i]==gnum){
				c++;
				fout<<heads[i]<<" "<<seqs[heads[i]][0]<<endl;
			}
		}
		fout<<"//";
		fout.close();
		if(c<2) return 0;
	}
	return 1;


}

void AHCluster::cleanData(){
	int j,i=1;
	while(i<heads.size())
	{

		for(j=0;j<i;j++){
			cout<<heads[i]<<"\t"<<heads[j]<<endl;
			int n1=0,n2=0;
			double sim = alignSeqSim(seqs[heads[j]][0],seqs[heads[i]][0],&n1,&n2);

			if(sim>0.95){
				if(n1>=n2 && heads[i]!=targetS){
					seqs.erase(heads[i]);
					heads.erase(heads.begin()+i);
				}
				else{
					seqs.erase(heads[j]);
					heads[j]=heads[i];
					heads.erase(heads.begin()+i);
				}
				break;
			}
		}
		if(j==i)i++;

	}
	for(int icol=0;icol<alignlen;){
		int irow;
		int gapn=0;
		for(irow=0;irow<heads.size();irow++){
			if(seqs[heads[irow]][0][icol]=='-') gapn++;
		}
		if(gapn>=irow*0.99){
			for(int i=0;i<heads.size();i++)
				seqs[heads[i]][0].erase(seqs[heads[i]][0].begin()+icol);
		}
		else icol++;
	}
	npoints = heads.size();
	alignlen = seqs[heads[0]][0].length();

	seqs_bk = seqs;
	heads_bk = heads;


	ofstream fout("/Users/mingmingliu/Documents/study/2014spring/class_variants/data/tmp.sto.aln");
	fout<<"# STOCKHOLM 1.0\n\n";
	for(int i=0;i<heads.size();i++)
		fout<<heads[i]<<" "<<seqs[heads[i]][0]<<endl;
	fout<<"//";
	fout.close();

}

double AHCluster::alignSeqSim(string seq1,string seq2,int* n1,int* n2){

	double dist = 0;
	int n=0;
	for(int i=0;i<seq1.length();i++){
		if(seq1[i]!='-') (*n1)++;
		if(seq2[i]!='-') (*n2)++;
		if(seq1[i]!='-' && seq2[i]!='-')
		{
			n++;
			if(seq1[i]==seq2[i]) dist++;

		}
	}
	if(n==0) return 1;

	return dist/n;
}


void AHCluster::getDistMat(){
	distmatrix.clear();
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

