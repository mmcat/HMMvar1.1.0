/*
 * Copyright (C) 2013 Virginia Tech

 * This file is part of HMMvar. HMMvar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Name: ScoreVar.cpp
 * Author: Mingming Liu
 */

#include <string>
#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "ScoreVar.h"
#include "CodeonAlphabeta.h"


using namespace std;

ScoreVar::ScoreVar(){

}

ScoreVar::~ScoreVar(){
	if (!tmp_dir_.empty()) {
			char rm_temp_dir_command[BUF_SIZE_MED];
			sprintf(rm_temp_dir_command, "rm -rf %s", tmp_dir_.c_str());
			int ret;
			ret = system(rm_temp_dir_command);
			if (ret != 0) {
				fprintf(stderr, "removing temporary directory failed (%d)\n", ret);
			}
		}

}



bool ScoreVar::SetQuerySequenceFromFastaFile(const char* id, string fasta_file){

	if(this->CreateTempDir()!=0) exit(-1);

	Log("loading query sequence from a FASTA file...\n", true);
	return query_seq_.SetSequenceFromFastaFile(id, fasta_file);
}

void ScoreVar::translate(const Sequence& nt, Sequence* aa,int type){
	aa->id_ = nt.id_;
	aa->def_ = nt.def_;
	string seqaa;
	string seqnt = nt.seq_;
        if(seqnt.empty()) {fprintf(stderr,"the sequence is empty\n");return;}
	transform(seqnt.begin(),seqnt.end(),seqnt.begin(), ::tolower);
	replace(seqnt.begin(),seqnt.end(),'u','t');
	set_Codeon();
	for(int i=0; i<seqnt.length()-codeonsize+1;i+=codeonsize){
		string triplet = seqnt.substr(i,codeonsize);
		if(triplet.compare(codeongap)==0){
			seqaa+=gap;
		}
		else if(Codeon.find(triplet) != Codeon.end()){
			seqaa+=Tables[type-1].substr(Codeon[triplet],1);
		}
		else seqaa+='X';

	}
	unsigned stop = seqaa.find('*');
	aa->seq_ = seqaa.substr(0,stop);


}


int ScoreVar::getScore(string hmm_path,string wtaa_file_name,string multi_align_file){


	if(hmm_output_file_name_.empty()){
		hmm_output_file_name_ = "/hmm_out_tmp";
	}

	buildHMM(hmm_path,multi_align_file);
//	string query_seq_aa;
//	this->nt2aa(query_seq_.seq_, &query_seq_aa);
	double bs_wild = searchHMM(hmm_path,wtaa_file_name);
	int len = query_seq_.seq_.length()/3;
	double null_pro = getNullPro(len);
	double pro_wild = null_pro*exp(bs_wild*CONST_LOG2/INTSCALE);

	if (variants_.size()==0) {
		fprintf(stderr,"No variants found to be predicted!\n");
		exit(-1);
	}

	if(score_out_file_name_.empty()) score_out_file_name_ = variant_file_name_+".out";
	FILE* fp=fopen(score_out_file_name_.c_str(),"w");

	for(map<string,vector<variant> >::iterator it = variants_.begin(); it!=variants_.end();it++){
		Sequence mut_seq_nt,mut_seq_aa;
		this->getMutantSeqFromVariantsSet(it->second,&mut_seq_nt);
<<<<<<< HEAD
        if(mut_seq_nt.seq_.empty()){fprintf(fp,"%s\t%.3f\n",it->first.c_str(),999.0);continue;}
		if(seq_type_=="prot"){mut_seq_aa = mut_seq_nt;}
		else{translate(mut_seq_nt,&mut_seq_aa,1);};
=======
                if(mut_seq_nt.seq_.empty()){fprintf(fp,"%s\t%.3f\n",it->first.c_str(),999.0);continue;}
		translate(mut_seq_nt,&mut_seq_aa,1);
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
		string mut_aa_file_name=tmp_dir_+"/query_aa_file_mt";
		mut_seq_aa.Print(mut_aa_file_name);
		double bs_mut = searchHMM(hmm_path,mut_aa_file_name);
		double mt_proba = null_pro*exp(bs_mut*CONST_LOG2/INTSCALE);
//		it->wt_proba_ = pro_wild;
<<<<<<< HEAD
		double odds = getOdds(pro_wild,mt_proba);
		double diffs = getDiffs(bs_wild,bs_mut);
		for(vector<variant>::iterator itt = it->second.begin();itt!=it->second.end();++itt){
			if(itt->type=="snp")
				fprintf(fp,"%s\t%.15f\t%s\n",it->first.c_str(),diffs,itt->commends_.c_str());
			else
				fprintf(fp,"%s\t%.15f\t%s\n",it->first.c_str(),odds,itt->commends_.c_str());
		}
=======
//		double odds = getOdds(pro_wild,mt_proba);
		double diffs = getDiffs(bs_wild,bs_mut);
		fprintf(fp,"%s\t%.3f\n",it->first.c_str(),diffs);
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d

	}



	fclose(fp);


	return 0;
}

int ScoreVar::getVariants(string filename){
<<<<<<< HEAD
    variant_file_name_ = filename;
=======
        variant_file_name_ = filename;
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
	FILE* fp = fopen(filename.c_str(),"r");
	if(fp==NULL){printf("Cannot find variants file:No such file!\n");return -1;}
	char buffer[BUF_SIZE_MED];
	char *a;
	int rn = 0;
	while(fgets(buffer,BUF_SIZE_MED,fp)!=NULL){
<<<<<<< HEAD
		printf("%s\n",buffer);
		if(buffer[0]=='#')
			continue;
=======
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
		rn++;
		variant var;
		a = strtok(buffer,"\t\n ");
		var.set_id_ = string(a);
		a = strtok(NULL,"\t\n");
		var.varinat_str_ = string(a);
<<<<<<< HEAD
		string::size_type pos = var.varinat_str_.find('/');
		if(pos==string::npos){printf("Error: Variant format is not correct:%s\n",var.varinat_str_.c_str());exit(0);}
		var.wild = var.varinat_str_.substr(0,pos);
		var.mutant = var.varinat_str_.substr(pos+1);
		if(var.wild!="-"&&var.wild.length()==1&&var.mutant!="-"&&var.mutant.length()==1)
			var.type = "snp";
		else var.type = "indel";

		a = strtok(NULL,"\t\n ");
		if(!a) {printf("Error: File format is not correct at row %d (skipped): %s\n",rn,filename.c_str());continue;}
		var.pos = atof(a);
		var.variant_id_ = this->query_seq_.id_;
		a = strtok(NULL,"\n");
		if(a) var.commends_=string(a);
=======
		a = strtok(NULL,"\t\n ");
		if(!a) {printf("Error: File format is not correct at row %d: %s\n",rn,filename.c_str());exit(0);}
		var.pos = atof(a);
		var.variant_id_ = this->query_seq_.id_;
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
		if(variants_.find(var.set_id_)==variants_.end()){
			vector<variant> var_set;
			var_set.push_back(var);
			variants_.insert(pair<string,vector<variant> >(var.set_id_,var_set));
		}
		else{
			variants_[var.set_id_].push_back(var);

		}

		}

	return 0;
}

double ScoreVar::getNullPro(int len){
	return exp(len*log(P1)+log(1.0-P1));
}

void ScoreVar::getMutantSeqFromVariants(variant var, Sequence* mut_seq){
	mut_seq->id_=query_seq_.id_;
	mut_seq->def_ = query_seq_.def_;
	mut_seq->seq_ = query_seq_.seq_;
	unsigned pos = var.varinat_str_.find('/');
	string w = var.varinat_str_.substr(0,pos);
	string m = var.varinat_str_.substr(pos+1);
	string s = query_seq_.seq_.substr(var.pos-1,w.size());
	if(w.compare(s)!=0){
		printf("Warning:The variant %s %d is not match with the sequence! Please check the variant.\n", var.varinat_str_.c_str(),var.pos);
	}
	if(m[0]=='-') mut_seq->seq_.erase(var.pos-1,w.size());
	else
		mut_seq->seq_.replace(var.pos-1,w.size(),m);

}

bool myComp(variant var1,variant var2){
	return var1.pos<var2.pos;
}

void ScoreVar::getMutantSeqFromVariantsSet(vector<variant> &var, Sequence* mut_seq){
	mut_seq->id_ = query_seq_.id_;
	mut_seq->def_ = query_seq_.def_;
	string str= query_seq_.seq_;
	sort(var.begin(),var.end(),myComp);
	int start = 0;
	for(vector<variant>::iterator it=var.begin();it!=var.end();++it){
<<<<<<< HEAD
//		string::size_type pos = it->varinat_str_.find('/');
//		if(pos==string::npos){printf("Error: Variant format is not correct:%s\n",it->varinat_str_.c_str());exit(0);}
		string w = it->wild;/*it->varinat_str_.substr(0,pos);*/
		string m = it->mutant;/*it->varinat_str_.substr(pos+1);*/
=======
		string::size_type pos = it->varinat_str_.find('/');
		if(pos==string::npos){printf("Error: Variant format is not correct:%s\n",it->varinat_str_.c_str());exit(0);}
		string w = it->varinat_str_.substr(0,pos);
		string m = it->varinat_str_.substr(pos+1);
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
		string s = query_seq_.seq_.substr(it->pos-1,w.size());
		if(w.compare(s)!=0&&w!="-"){
				printf("Warning:The variant %s %d is not match with the sequence! Please check the variant.\n", it->varinat_str_.c_str(),it->pos);
			}
<<<<<<< HEAD
//		else {
=======
		else {
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
			string seg = str.substr(start,it->pos-start-1);
			if(w=="-" && m!="-"){
				seg = str.substr(start,it->pos-start);
				mut_seq->seq_.append(seg+m);
			}
			else if(w!="-" && m=="-") mut_seq->seq_.append(seg);
			else mut_seq->seq_.append(seg+m);
			start = it->pos+w.length()-1;

<<<<<<< HEAD
//		}
=======
		}
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d

	}
	mut_seq->seq_.append(str.substr(start));
}



double ScoreVar::getOdds(double p1,double p2){
	if(p2==0 || p1==1) return 999;
	else return p1*(1-p2)/p2*(1-p1);

}

double ScoreVar::getDiffs(double b1,double b2){
	return b1-b2;
}
int ScoreVar::buildHMM(string hmm_path,string multi_align_file){
	char hmmbuild_cmd[BUF_SIZE_MED];
	string hmm_out_file = tmp_dir_+hmm_output_file_name_;
	if(hmm_path.empty()) hmm_path = HMMPATH;
	sprintf(hmmbuild_cmd,"%shmmbuild --amino --cpu 10 %s %s ",hmm_path.c_str(),hmm_out_file.c_str(),multi_align_file.c_str());

	Log("build hidden Markov model...\n",true);

	int return_code = system(hmmbuild_cmd);
		if(return_code!=0){
			fprintf(stderr, "hmmbuild failed (exit:%d)\n", return_code);
					exit(-1);

		}
		return 0;


}

double ScoreVar::searchHMM(string hmm_path,string aa_file_name){
	char hmmsearch_cmd[BUF_SIZE_MED];
	string hmm_out_file = tmp_dir_+hmm_output_file_name_;
	if(hmm_path.empty()) hmm_path = HMMPATH;
	sprintf(hmmsearch_cmd,"%shmmsearch %s %s|grep -A4 \"Scores for complete sequences\"|tail -n 1",hmm_path.c_str(), hmm_out_file.c_str(),aa_file_name.c_str());
	Log("maching hidden Markov model...\n",true);

	FILE *PIPE = popen(hmmsearch_cmd,"r");

	if(PIPE==NULL){
		fprintf(stderr, "hmmsearch failed (exit:%s)\n", hmmsearch_cmd);
		exit(-1);
	}

	char buf[BUF_SIZE_LARGE];
	char* str;
	fgets(buf,BUF_SIZE_MED,PIPE);
	str = strtok(buf, " \r\n");
	str = strtok(NULL," \r\n");
        pclose(PIPE);
	if(str==NULL) return 0;
	return atof(str);
}

int ScoreVar::CreateTempDir(){
	char temp_dir[BUF_SIZE_MED];
		sprintf(temp_dir, "%s/hmmVarXXXXXX", P_tmpdir);

		if (mkdtemp(temp_dir) == NULL) {
			fprintf(stderr, "Error in creating temporary directory\n");
			return -1;
		}

		tmp_dir_ = temp_dir;

		return 0;
}

int ScoreVar::parseBlastResults(string blastdb_cmd,string blastdb,string subject_sequence_file_name,string blast_out_file_name,double cutoff){
	if (subject_sequence_file_name.empty()) return -1;


	char tmp_file_id[BUF_SIZE_MED];
	FILE* fp = fopen(blast_out_file_name.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << blast_out_file_name << endl;
		return -1;
	}

	char buf[BUF_SIZE_MED];

	char* db_id;
	char* gi;
	char* pident;
	int num_seqs = 0;

	sprintf(tmp_file_id, "%s/ids", this->tmp_dir_.c_str());
	FILE* fp_id_out = fopen(tmp_file_id, "w");

	char last_id[BUF_SIZE_MED] = "";

	// get subject seqs
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_MED, fp) == NULL)
			break;

		db_id = strtok(buf, "\t\n");
		gi = strtok(NULL, "\t\n");
		pident = strtok(NULL,"\t\n");

		char id[BUF_SIZE_MED]="";

		if (db_id == NULL || gi == NULL) {
<<<<<<< HEAD
			Log("Parsing Error,Skipped!\n", true);
//			return -1;
			continue;
=======
			Log("Parsing Error\n", true);
			return -1;
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
		}

		strcpy(id, db_id);

		if (db_id != NULL && gi != NULL) {
			if (atof(pident)<cutoff*100 || strcmp(last_id, id) == 0) {
				continue;
			}

			fprintf(fp_id_out, "%s\n", id);
			strcpy(last_id, id);
			num_seqs++;
		} else {
			Log("Wrong entry in blast out\n", true);
			break;
		}
	}

	fclose(fp);
	fclose(fp_id_out);

	this->CreateFastaFileUsingBlastdbcmd(blastdb_cmd,blastdb,subject_sequence_file_name,tmp_file_id);

	return num_seqs;

}

int ScoreVar::CreateFastaFileUsingBlastdbcmd(string blastdb_cmd,string blastdb,string subject_seq_file_name,string file_id){
	char blastdbcmd[BUF_SIZE_MED];
	if(blastdb_cmd.empty()) blastdb_cmd = BLASTDBCMD;
	if(blastdb.empty()) blastdb = BLASTDB;
		sprintf(blastdbcmd, "%s -db %s -dbtype 'prot' -entry_batch %s > %s",
				blastdb_cmd.c_str(),
				blastdb.c_str(),
				file_id.c_str(),
				subject_seq_file_name.c_str()
				);
		system(blastdbcmd);

	return 0;
}

int ScoreVar::setHomoseq(string blastdb_cmd,string blastdb,string subject_sequence_file_name,string blast_out_file_name,double cutoff){


	Log("filtering related sequences...\n", true);

	return this->parseBlastResults(blastdb_cmd,blastdb,subject_sequence_file_name,blast_out_file_name,cutoff);

}

