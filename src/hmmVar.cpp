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
	opt.SetOptions(argc,argv); // processing the arguments from command  line.

	ScoreVar score;
    score.SetQuerySequenceFromFastaFile(opt.query_file_name_.c_str(),opt.query_file_name_);

    Sequence wtaa;     // wild type amino acid
    int type=1;        //different type of codeon table
    string wtaa_query_file_name;
// step 1: target seed protein (query protein sequence). The query sequence could be protein sequence or cDNA sequence.
// Translate the query sequence to protein sequence if it is a cDNA sequence.
    Log("======Step 1: target seed protein (query protein sequence) ======\n",true);
    if(opt.seq_type_ != "" && opt.seq_type_=="prot"){
    	wtaa_query_file_name = opt.query_file_name_;
    	score.seq_type_ = "prot";
    }
    else{
    	score.translate(score.query_seq_, &wtaa,type);
    	wtaa_query_file_name=score.tmp_dir_+"/query_aa_file_wt";
    	wtaa.Print(wtaa_query_file_name);
    }
// step 2: search for homologous sequences by psiblast. The parameters of psiblast are set by default (Please see the paper for details).
// If there is a blast result file exists in advance, this step is skiiped automatically.
    Log("======Step 2: search for homologous sequences by psiblast ======\n",true);
	blastSearch blast(wtaa_query_file_name, opt.blast_output_file_name_,opt.subject_sequences_file_name_,opt.save_blastout_,score.tmp_dir_);

	if(opt.blast_output_file_name_.empty() && opt.hmmbuild_output_file_name_.empty()){
        
		blast.searchHomoSeq(opt.blastdb_file_name_,opt.psiblast_command_);
    }

// step 3: parse blast results and make multiple sequence alignment (MSA) by muscle.
	Log("======Step 3: make multiple sequence alignment (MSA) ======\n",true);
	MultiAlign align(score.tmp_dir_,opt.save_muscle_);
    if(opt.hmmbuild_output_file_name_.empty()){

        if (score.setHomoseq(opt.blastdbcmd_command_,opt.blastdb_file_name_,blast.subject_sequences_file_name_,blast.blast_output_file_name_,CUTOFF)!=-1 \
    		  && align.make_multi_align(opt.muscle_command_,blast.subject_sequences_file_name_,opt.muscle_output_file_name_)!=-1)

    	   align.fasta2stockholm(align.align_output_file_fasta_);
    }


    score.getVariants(opt.variants_file_name_);  // get the mutation (alleles and positions)
    
// step 4&5: build hidden Markov model and scoring
    Log("======Step 4 & 5: make hidden Markov model and scoring ======\n",true);
	score.getScore(opt.hmmbuild_output_file_name_,opt.hmmer_command_,wtaa_query_file_name,align.align_output_file_stockholm_);

	cout<<"done!"<<endl;


	return 0;
}

