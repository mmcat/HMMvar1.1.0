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
 * Name: MultiAlign.cpp
 * Author: Mingming Liu
 */

#include <string>
#include <stdio.h>
#include <map>
#include <vector>

#include "MultiAlign.h"


using namespace std;

MultiAlign::MultiAlign(string path,string alignout_save_name){
	output_path_ = path;
        alignout_save_file_name_ = alignout_save_name;

}

MultiAlign::~MultiAlign(){

}



int MultiAlign::make_multi_align(string muscle_cmd, string align_input_file,string align_output_file){
	align_output_file_fasta_ = align_output_file;
	if(align_output_file.empty()) align_output_file_fasta_ = output_path_+"/align_out_tmp";
	return this->runMuscle(muscle_cmd,align_input_file);

}

int MultiAlign::runMuscle(string muscle_cmd,string align_input_file){

	char muscle_command[BUF_SIZE_MED];
	if(muscle_cmd.empty()) muscle_cmd = MUSCLE;
	sprintf(muscle_command,"%s -in %s -out %s -maxiters 1 -diags -sv -distance1 kbit20_3", muscle_cmd.c_str(),
			align_input_file.c_str(),align_output_file_fasta_.c_str());

	Log("making multiple sequence alignment...\n", true);
	int return_code = system(muscle_command);
	if(return_code!=0){
		fprintf(stderr, "muscle failed (exit:%d)\n", return_code);
				exit(-1);

	}
	if (!alignout_save_file_name_.empty()) {
					char cp_command[BUF_SIZE_MED];
					sprintf(cp_command, "cp %s %s", align_output_file_fasta_.c_str(), alignout_save_file_name_.c_str());
					// will replace system
					int ret;
					ret = system(cp_command);
					if (ret != 0) {
						fprintf(stderr, "saving blastout file failed (%d)\n", ret);
					}
				}
	return 0;
}

int MultiAlign::fasta2stockholm(string oldfile){
	if(oldfile.empty())
		return -1;
	if(align_output_file_stockholm_.empty()){
		align_output_file_stockholm_ = output_path_+"/align_out_stockholm";
	}
	FILE* fpout=fopen(align_output_file_stockholm_.c_str(),"w");
	FILE* fpin = fopen(oldfile.c_str(),"r");

	fprintf(fpout,"# STOCKHOLM 1.0\n\n");

	map <string, vector <string> > seq_map;
	char buf[BUF_SIZE_LARGE];
	string name;

	int numblock;

	while(fgets(buf,BUF_SIZE_LARGE,fpin) != NULL){
		if(buf[0]=='>') {
			name.assign(buf,1,15);
			if(seq_map.find(name)!=seq_map.end()) seq_map[name].clear();
			seq_map.insert(make_pair(name,vector<string>()));

		}
		else{

			try{
				string seq = string(buf);
				seq_map.at(name).push_back(seq);
//				throw name.c_str();
			}
			catch(char* name){
				printf("Cannot find %s\n",name);

			}
		}
	}
	map<string,vector<string> >::iterator it = seq_map.begin();
	numblock = seq_map[it->first].size();
	for(int i=0; i<numblock;i++){
		for (map<string,vector<string> >::iterator it = seq_map.begin(); it!=seq_map.end();it++ ){

			fprintf(fpout,"%s %s",it->first.c_str(),it->second[i].c_str());
		}
		fprintf(fpout,"\n\n");
	}

	fprintf(fpout,"//\n");


	fclose(fpin);
	fclose(fpout);
	return 0;
}
