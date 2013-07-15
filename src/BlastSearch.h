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
 * Name: BlastSearch.h
 * Author: Mingming Liu
 */

#ifndef BLASTSEARCH_H_
#define BLASTSEARCH_H_

#include <string>

#include "Common.h"

using namespace std;

class blastSearch{
public:
	blastSearch(string blast_query,string blast_output,string subject_seqs,string blastout_save_name,string path);
	virtual ~blastSearch();
    string blastout_save_file_name_ ;
	string blast_output_file_name_;
//	string blast_db_name_;
	string blast_query_file_name_;
//	string psiblast_command_;
	string subject_sequences_file_name_;
	string temp_dir_;

	void searchHomoSeq(string blastdb_file_name,string psiblast_cmd);

private:

	int runBlast(string blastdb_file_name,string psiblast_cmd);




};


#endif /* BLASTSEARCH_H_ */
