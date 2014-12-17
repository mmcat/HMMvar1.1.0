/*
 * Kmeans.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: mingmingliu
 */
#include "Kmeans.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <limits.h>

using namespace std;

Kmeans::Kmeans(int k, MultiAlign* align){
	this->k = k;
	this->align = align;
	len = 0;
	targetG = 0;

}

void Kmeans::init(){
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
					seqs[def] = seqseg;
					heads.push_back(def);
				}
				else{
					seqs[def]+=seqseg;
				}
			}
		}
		int n=heads.size();
		map<int,bool> mp;
		srand (time(NULL));
		for(int i=0;i<k;i++){
			int m=rand()%n;
			while(mp.find(m)!=mp.end())
				m=rand()%n;
			mp[m]=true;
			centers.push_back(seqs[heads[m]]);
		}
	}
	else{
		fprintf(stderr,"Cannot open the alignment file!\n");
		exit(-1);
	}
	fin.close();
	labels=vector<int>(seqs.size(),-1);
	len = seqs.begin()->second.length();
	groups=vector<vector<int> >(k,vector<int>(1));

}

void Kmeans::updateCenter(int gnum){
	string newcenter="";

	for(int j=0;j<len;j++){
		map<char,int> cnt;

		for(int i=0;i<groups[gnum].size();i++){

			char aa = seqs[heads[groups[gnum][i]]][j];
			if(cnt.find(aa)==cnt.end()) {
				cnt[aa] = 1;
			}
			else {
				cnt[aa]++;
			}
		}
		int max = INT_MIN;
		char thisaa;
		for(map<char,int>::iterator it = cnt.begin();it!=cnt.end();it++){
			if(it->second>max){
				max = it->second;
				thisaa = it->first;
			}
		}
		newcenter+=thisaa;

	}
	centers[gnum] = newcenter;
}

void Kmeans::runKmeans(){
	int changes = seqs.size();
	int stop = 3;
    cout<<heads.back()<<endl;
	while(changes>stop){
		for(int i=0;i<groups.size();i++)
			groups[i].clear();
		changes = 0;
		int g;
		for(int i=0;i<heads.size();i++){
			int min = INT_MAX;
            if(i==heads.size()-1){
            	int mm = 0;
            	cout<<heads[i]<<"\n";
            }
			for(int j=0;j<k;j++){
				int d = editDist(seqs[heads[i]],centers[j]);
				if(d<min){
					min = d;
					g = j;
				}
			}
			groups[g].push_back(i);
			if(g!=labels[i]) {
				changes++;
				labels[i] = g;
			}

			if(this->targetS.find(heads[i])!=std::string::npos){
				targetG=g;
			}

		}

		for(int i=0;i<groups.size();i++){
			updateCenter(i);
		}


	}

}

int Kmeans::editDist(string seq1, string seq2){


	int dist = 0;

	for(int i=0;i<len;i++){
		if(seq1[i]!=seq2[i]) dist++;
	}

	return dist;
}

void Kmeans::printCluter(int gnum,string filename){
	ofstream fout(filename.c_str());
	int block_seq_len = 60;
	int MIN_DEF_LEN = 15;
	int blocks = ceil((double)len/block_seq_len);

	if(fout.is_open()){
		fout<<"# STOCKHOLM 1.0\n\n";
		string def;
		string l_def;
		for(int i=0;i<blocks-1;i++){
			for(int j=0;j<groups[gnum].size();j++){
				def = heads[groups[gnum][j]];
				l_def = def;
				if(def.length()<MIN_DEF_LEN)
					l_def=def+string(MIN_DEF_LEN-def.length(),'X');
				fout<<l_def<<" "<<seqs[def].substr(i*block_seq_len,block_seq_len)<<"\n";
			}
			fout<<"\n\n";
		}

		for(int j=0;j<groups[gnum].size();j++){
			def = heads[groups[gnum][j]]; //TODO: add to main branch
			l_def = def;
			if(def.length()<MIN_DEF_LEN)
				l_def=def+string(MIN_DEF_LEN-def.length(),'X');
			fout<<l_def<<" "<<seqs[def].substr((blocks-1)*block_seq_len,block_seq_len)<<"\n";
		}

		fout<<"\n\n//";

	}
	fout.close();

}



