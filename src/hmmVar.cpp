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
#include <fstream>

#include "BlastSearch.h"
#include "MultiAlign.h"
#include "ScoreVar.h"
#include "Option.h"
#include "Kmeans.h"
#include "AHCluster.h"

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
/*    if(opt.hmmbuild_output_file_name_.empty()){

        if (score.setHomoseq(opt.blastdbcmd_command_,opt.blastdb_file_name_,blast.subject_sequences_file_name_,blast.blast_output_file_name_,CUTOFF)!=-1 \
    		  && align.make_multi_align(opt.muscle_command_,blast.subject_sequences_file_name_,opt.muscle_output_file_name_)!=-1)

    	   align.fasta2stockholm(align.align_output_file_fasta_);


    }*/
    Log("======Step 4: clustering the MSA======\n",true);
    align.align_output_file_stockholm_="/Users/mingmingliu/Documents/study/2014spring/class_variants/data/PTPRD.uniprot90.aln.stk";
    AHCluster ahc(&align);
    ahc.targetS = score.query_seq_.def_;
    ahc.init();
    ahc.cleanData();
//	ahc.runAHC();
//	ahc.cuttree();


/*
    Kmeans clustering(4,&align);
    clustering.targetS = score.query_seq_.def_;
    clustering.init();
    clustering.runKmeans();
*/



    score.getVariants(opt.variants_file_name_);  // get the mutation (alleles and positions)
    score.query_seq_bk = score.query_seq_;
    
    Log("======Step 5: make hidden Markov model and scoring ======\n",true);
//    while(score.variants_bk.size()>0){
    ofstream result_out("/Users/mingmingliu/Documents/study/2014spring/class_variants/data/pdprd_cluster_result.txt");
    for(map<string,vector<variant> >::iterator it=score.variants_bk.begin();it!=score.variants_bk.end();++it){
    	if(it->second.empty()) continue;

    	variant seed = *(it->second.begin());
    	int pos = seed.pos;
    	int a_start,a_end;
    	/*update alignment*/

    	int align_pos = ahc.selectRange(pos,&a_start,&a_end);

    	align.align_output_file_stockholm_ = "/Users/mingmingliu/Documents/study/2014spring/class_variants/data/tmp1.sto.aln";

    	/*update query seq*/
    	score.query_seq_ = score.query_seq_bk;
    	score.query_seq_.seq_ = score.query_seq_.seq_.substr(pos-(align_pos-a_start),ahc.alignlen);
    	int s_start = pos-(align_pos-a_start); // base 0
    	int s_end = pos+(a_end-align_pos); // base 0

    	wtaa_query_file_name = score.tmp_dir_+"/query_aa_file_wt";
    	ofstream fout(wtaa_query_file_name.c_str());
    	fout<<score.query_seq_.def_<<endl;
    	fout<<score.query_seq_.seq_;
    	fout.close();



    	/*update variants*/
    	score.getVarInRange(s_start,s_end);

    	ahc.runAHC();
    	ahc.cuttree();

       string group_align_file_name = score.tmp_dir_+"/group_align_output_file_stockholm";
       string filename = group_align_file_name;
       vector<vector<double> > grouped_wtscores;
       vector<vector<double> > grouped_mtscores;
       vector<vector<string> > grouped_ids;

       score.wtscores.clear();
       score.mtscores.clear();
       score.ids.clear();

       score.getScore(opt.hmmbuild_output_file_name_,opt.hmmer_command_,wtaa_query_file_name,align.align_output_file_stockholm_);
       grouped_wtscores.push_back(score.wtscores);
       grouped_mtscores.push_back(score.mtscores);
       grouped_ids.push_back(score.ids);

       score.wtscores.clear();
       score.mtscores.clear();
       score.ids.clear();

       string labeled_all_seqs_filename = filename+"_labeled";
       ahc.printCluster(-1,labeled_all_seqs_filename);

       group_align_file_name=filename+"_target"+to_string(ahc.targetG);
       ahc.printCluster(ahc.targetG,group_align_file_name);
       score.getScore(opt.hmmbuild_output_file_name_,opt.hmmer_command_,wtaa_query_file_name,group_align_file_name);
       grouped_wtscores.push_back(score.wtscores);
       grouped_mtscores.push_back(score.mtscores);
       grouped_ids.push_back(score.ids);

       for(int i=0;i<ahc.groups;i++){
    	   if(i==ahc.targetG) continue;
    	   group_align_file_name=filename+to_string(i);
    	   if(ahc.printCluster(i,group_align_file_name)==0) continue;
//    	   ahc.printCluster(i,group_align_file_name);
    	   score.wtscores.clear();
           score.mtscores.clear();
           score.ids.clear();
           score.getScore(opt.hmmbuild_output_file_name_,opt.hmmer_command_,wtaa_query_file_name,group_align_file_name);
           grouped_wtscores.push_back(score.wtscores);
           grouped_mtscores.push_back(score.mtscores);
           grouped_ids.push_back(score.ids);
       }

       score.query_seq_ = score.query_seq_bk;

  /*     for(int i=0;i<grouped_wtscores.size();i++){
       	for(int j=0;j<grouped_wtscores[i].size();j++){
       		cout<<"("<<grouped_ids[i][j]<<":"<<grouped_wtscores[i][j]<<","<<grouped_mtscores[i][j]<<")"<<",";
       	}
       	cout<<"\n";
       }*/
//       cout<<ahc.groups<<"\t"<<grouped_wtscores[0].size()<<endl;

       for(int j=0;j<grouped_wtscores[0].size();j++){
    	   result_out<<grouped_ids[0][j]<<"\t";
    	   int i;
    	   for(i=0;i<grouped_wtscores.size()-1;i++)

    		   result_out<<grouped_wtscores[i][j]-grouped_mtscores[i][j]<<"\t";

    	   result_out<<grouped_wtscores[i][j]-grouped_mtscores[i][j]<<"\n";

       }
    }
    result_out.close();
   	cout<<"All steps are done!"<<endl;


   	return 0;
}

