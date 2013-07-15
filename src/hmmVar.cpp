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
 * Name: hmmVar.cpp
 * Author: Mingming Liu
 */


#include <iostream>
#include <string>

#include "BlastSearch.h"
#include "MultiAlign.h"
#include "ScoreVar.h"
#include "Option.h"

using namespace std;

int main(int argc, char** argv){

	Option opt;
	opt.SetOptions(argc,argv);

	ScoreVar score;
    score.SetQuerySequenceFromFastaFile(opt.query_file_name_.c_str(),opt.query_file_name_);

    Sequence wtaa;
    int type=1;        //different type of codeon table
<<<<<<< HEAD
    string wtaa_query_file_name;
// step 1: find seed protein (query protein sequence)
    if(opt.seq_type_ != "" && opt.seq_type_=="prot"){
    	wtaa_query_file_name = opt.query_file_name_;
    	score.seq_type_ = "prot";
    }
    else{
    	score.translate(score.query_seq_, &wtaa,type);
    	wtaa_query_file_name=score.tmp_dir_+"/query_aa_file_wt";
    	wtaa.Print(wtaa_query_file_name);
    }
// step 2: search for homologous sequences by psiblast
	blastSearch blast(wtaa_query_file_name, opt.blast_output_file_name_,opt.subject_sequences_file_name_,opt.save_blastout_,score.tmp_dir_);

	if(opt.blast_output_file_name_.empty())
		blast.searchHomoSeq(opt.blastdb_file_name_,opt.psiblast_command_);
=======
// step 1: find seed protein (query protein sequence)
    score.translate(score.query_seq_, &wtaa,type);
    string wtaa_query_file_name=score.tmp_dir_+"/query_aa_file_wt";
    wtaa.Print(wtaa_query_file_name);
// step 2: search for homologous sequences by psiblast
	blastSearch blast(wtaa_query_file_name, opt.blast_output_file_name_,opt.subject_sequences_file_name_,opt.save_blastout_,score.tmp_dir_);
	blast.searchHomoSeq(opt.blastdb_file_name_,opt.psiblast_command_);
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d

// step 3: make multiple sequence alignment (MSA)
	MultiAlign align(score.tmp_dir_,opt.save_muscle_);
    if (score.setHomoseq(opt.blastdbcmd_command_,opt.blastdb_file_name_,blast.subject_sequences_file_name_,blast.blast_output_file_name_,CUTOFF)!=-1 \
    		&& align.make_multi_align(opt.muscle_command_,blast.subject_sequences_file_name_,opt.muscle_output_file_name_)!=-1)

    	align.fasta2stockholm(align.align_output_file_fasta_);


    score.getVariants(opt.variants_file_name_);  // get the mutation (alleles and positions)
// step 4&5: make hidden markov model and scoring
	score.getScore(opt.hmmer_command_,wtaa_query_file_name,align.align_output_file_stockholm_);

	cout<<"done!"<<endl;


	return 0;
}

