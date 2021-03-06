/*
 * =====================================================================================
 *
 *       Filename:  compensatoryIndels.cpp
 *
 *    Description:  find all compensatory set of variants for a given gene
 *
 *        Version:  1.0
 *        Created:  05/01/2013 09:35:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingming Liu 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "omp.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <time.h>
#include "math.h"

using namespace std;
const int BIGSUM=4;
const int MAXMUT=5;
const int MAXSETSIZE=10;
const int MAXSETNUM=20;
const int RANDNUM=200;
const int RANGE = 20; //indels occur within 20 bp
int base=0;
struct var{
	int diff;
	string allele;
	string mutid;
	int pos;
};

void dynFindSet(vector<var>& arrayVar);
vector< vector<int> > traceBack(vector< vector<bool> > & T,vector<var> &arrayVar);
void readVar(ifstream &fin, vector<vector<var> > &arrayVar);
vector<string> split(string ss,string seperator);
void outPut(vector< vector<bool> > &T,vector<var> &arrayVar,vector< vector<int> >& cur_sets,vector<int> s,int target,int srow,int scol);
void findComp (vector<var> &allVar,vector<var> &comp,int index);
bool mysort(var var1,var var2);

int main(int argc,const char* argv[]){
	ifstream fin;
	if(!argv[1]){
		cerr<<"Please input the variant file!"<<endl;
		exit(0);
	}
	fin.open(argv[1]);
	vector< vector<var> > arrayVarSet;
	readVar(fin,arrayVarSet);
	for(int i=0;i<arrayVarSet.size();i++){
		vector<var> arrayVar = arrayVarSet[i];
		cout<<"#find compensary for var:"<<arrayVar[0].mutid<<'\t'<<arrayVar[0].allele<<'\t'<<arrayVar[0].pos<<endl;
		dynFindSet(arrayVar);
	}
	return 0;

}

void dynFindSet(vector<var>& arrayVar){
	int sump=0;
	int sumn=0;
	vector< vector<bool> > table;
	for(vector<var>::iterator it=arrayVar.begin();it!=arrayVar.end();++it){
		if(it->diff>0)
			sump+=it->diff;
		else
			sumn+=it->diff;
	}
	base = abs(sumn);
	if(sump>BIGSUM) sump = BIGSUM;
	if(base>BIGSUM) base = BIGSUM;
	int sum = sump+base;
	for(int i=0;i<arrayVar.size();i++){
		vector<bool> t(sum+1,false);
		table.push_back(t);
	}
	for(int i=0;i<arrayVar.size();i++){
		if(i==0){
			for(int j=0;j<table[i].size();j++)
				if(j-base==arrayVar[i].diff){
					table[i][j]=true;
					break;
				}
		}
		else{
			for(int j=0;j<table[i].size();j++){
				int diff = arrayVar[i].diff;
				if(table[i-1][j] || diff==j-base||(j-diff>=0&&j-diff<=sum&&table[i-1][j-diff])) 
					table[i][j] = true;
			}
		}
	}
	vector< vector<int> > sets = traceBack(table,arrayVar);
	ofstream fout("table.txt");
	if(fout.is_open()){
		for(int i=0;i<table.size();i++){
			for(int j=0;j<table[i].size();j++)
				fout << table[i][j]<<"\t";
			fout<<"\n";
		}
		
  //		fout<<"\n\n\n";
  		for(int i=0;i<sets.size();i++){
			for(int j=0;j<sets[i].size();j++)
				cout<<arrayVar[0].mutid<<"_"<<i<<"\t"<<arrayVar[sets[i][j]].allele<<"\t"<<arrayVar[sets[i][j]].pos<<"\t"<<arrayVar[sets[i][j]].mutid<<"\n";
			
		}
	}


}

vector< vector<int> > traceBack(vector< vector<bool> > & T,vector<var> &arrayVar){
	int col=T[0].size();
	int row = T.size();
	vector< vector<int> > sets;
#pragma omp parallel
	{
		#pragma omp for
		for(int ele=-base-(-base%3);ele<col-base;ele+=3){
			int i = row-1;
			int j = ele+base;
			if(!T[i][j]){cout<<"#no subset for target:"<<ele<<endl;ele+=3;continue;}
			vector< vector<int> > cur_sets;
			vector<int> s;
			outPut(T,arrayVar,cur_sets,s,ele,i,j);
			#pragma omp critical
			{
				for(int i=0;i<cur_sets.size();i++){
					sets.push_back(cur_sets[i]);
				}
			}
		}
	}	
	return sets;

}

void  outPut(vector< vector<bool> > &T,vector<var> &arrayVar, vector< vector<int> > &cur_sets,vector<int> s,int target,int srow,int scol){
	int row = T.size();
	int col = T[0].size();
	if(srow<0||scol<0||!T[srow][scol]||cur_sets.size()>=MAXSETNUM) return;

	vector<int> ss = s;
	if(srow==0&&scol-base==arrayVar[srow].diff&&s.size()<MAXSETSIZE){
		s = ss;
		s.push_back(srow);
		cur_sets.push_back(s);
/*		#pragma omp critical
		{
			
			for(int i=0;i<s.size();i++){
				cout<<i<<"\t"<<arrayVar[s[i]].allele<<"\t"<<arrayVar[s[i]].pos<<"\t"<<arrayVar[s[i]].mutid<<endl;
			}
//			cout<<"target:"<<target<<endl;
			
		}*/
	}
	if(scol-arrayVar[srow].diff<col && scol-arrayVar[srow].diff>=0 && srow>0 && T[srow-1][scol-arrayVar[srow].diff]){
		s = ss;
		s.push_back(srow);
		outPut(T,arrayVar,cur_sets,s,target,srow-1,scol-arrayVar[srow].diff);
	}
        s = ss;	
	outPut(T,arrayVar,cur_sets,s,target,srow-1,scol);
	
}

bool mysort(var var1,var var2){
	return var1.pos<var2.pos;
}

void findComp (vector<var> &allVar,vector<var> &comp,int index){
	int curPos = allVar[index].pos;
	comp.push_back(allVar[index]);
	
	for(int i=index-1;i>=0;--i){
		if(curPos - allVar[i].pos<=RANGE && allVar[i].pos!=comp.back().pos)
			comp.push_back(allVar[i]);
		else if(allVar[i].pos==comp.back().pos) continue;
		else break;
	}
	int big = comp[0].pos;
	for(int j=index+1;j<allVar.size();j++){
		if(allVar[j].pos - curPos<=RANGE && allVar[j].pos!=big){
			comp.push_back(allVar[j]);
			big = comp.back().pos;
		}
		else if(allVar[j].pos==big) continue;
		else break;
	}
}

void readVar(ifstream &fin, vector< vector<var> > &arrayVarSet){
	string line;
	vector<var> allVar;
	if(fin.is_open()){
		while(getline(fin,line)){
//			cout<<line<<endl;
			var v;
			string wt,mt;
			vector<string> tokens = split(line,"\t");
			v.mutid = tokens[0];
			v.allele = tokens[1];
			v.pos = atof(tokens[2].c_str());
			string::size_type p = v.allele.find("/");
			if(p==string::npos){
				cerr<<"variant format error:"<<v.allele<<" MUT_ID:"<<v.mutid<<endl;
			}
			else{
				wt = v.allele.substr(0,p);
				if(wt=="-") wt="";
				mt = v.allele.substr(p+1);
				if(mt=="-") mt="";
				int diff=wt.length()-mt.length();
				v.diff = diff;
			}
			//small insertions and deletions
			if(v.diff%3!=0&&wt.length()<MAXMUT&&mt.length()<MAXMUT)
				allVar.push_back(v);

		}
		fin.close();

	}
	else{
		cerr<<"The variant file open failed:NO such file or directory!"<<endl;
		return;
	}

	//random sampling
/*	vector<var> rnd_allVar;
	int tot = allVar.size();
	int* visited = new int[tot];
	for(int i=0;i<tot;i++) visited[i] = 0;
	srand(time(NULL));
	if(tot>RANDNUM){
		int n=0;
		while(n<RANDNUM){
			int index = rand()%tot;
			if(visited[index]==1) continue;
			else{
				visited[index] = 1;
				rnd_allVar.push_back(allVar[index]);
				n++;
			}
		}
	}*/
//	int mm=0;
   	sort(allVar.begin(),allVar.end(),mysort);
	for(int i=0;i<allVar.size();i++){
		vector<var> comp;
		findComp(allVar,comp,i);
		if(comp.size()>0)
			arrayVarSet.push_back(comp);
	}

}

vector<string> split(string ss,string seperator){
        vector<string> tokens;
        int s=0;
        string str = ss;
        while(s<str.length()){
                string::size_type p=str.find(seperator,s);
                if(p!=string::npos) {
                        tokens.push_back(str.substr(s,p-s));
                        s=p+seperator.length();
                }
                else {
                        p=str.find("\n",s);
                        if(p!=string::npos)
                                tokens.push_back(str.substr(s,p-s));
                        else tokens.push_back(str.substr(s));
                        break;
                }
        }
        return tokens;
}


