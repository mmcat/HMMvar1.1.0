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
 * Name: ScoreVar.h
 * Author: Mingming Liu
 */

#ifndef SCOREVAR_H_
#define SCOREVAR_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Sequence.h"
#include "variant.h"


using namespace std;

class ScoreVar{
public:
	ScoreVar();
	virtual ~ScoreVar();

	Sequence query_seq_;
	Sequence query_seq_bk;
	vector < Sequence > subject_seqs_;
	vector<double> scores;
	vector<double> wtscores;
	vector<double> mtscores;
	vector<string> ids;
	map<string,vector <variant> > variants_;
	map<string,vector <variant> > variants_bk;
//	string hmm_input_file_name_;
	string hmm_output_file_name_;
	string save_hmm_file_name_;
        string variant_file_name_;
	string tmp_dir_;
	string score_out_file_name_;
	string seq_type_;

//	string hmm_command_;

//	void SetHMMParameters(string input,string output);
	bool SetQuerySequenceFromFastaFile(const char* id, string fasta_file);

	int CreateTempDir();
	int parseBlastResults(string blastdb_cmd,string blastdb,string subject_sequence_file_name,string blast_out_file_name,double cutoff);
	int setHomoseq(string blastdb_cmd,string blastdb,string subject_sequence_file_name,string blast_out_file_name,double cutoff);
	int getScore(string hmmbuild_file,string hmm_path,string wtaa_file_name,string multi_align_file);
	int getVariants(string filename);
	double getNullPro(int len);
	double getOdds(double p1,double p2);
        double getDiffs(double b1,double b2);
	void getMutantSeqFromVariants(variant var,Sequence* mut_seq);
	void getMutantSeqFromVariantsSet(vector<variant> &var,Sequence* mut_seq);
	void translate(const Sequence& aa, Sequence * nt,int type);
	void getVarInRange(int start,int end);


private:

	int SetSequences();
	int CreateFastaFileUsingBlastdbcmd(string blastdb_cmd,string blastdb,string subject_seq_file_name,string file_id);
	int buildHMM(string hmm_path,string multi_align_file);
	double searchHMM(string hmmbuild_file,string hmm_path,string aa_file_name);


};




#endif /* SCOREVAR_H_ */
