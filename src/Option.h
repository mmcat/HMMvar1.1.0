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
 * Name: Option.h
 * Author: Mingming Liu
 */

#ifndef OPTION_H_
#define OPTION_H_

#include<string>
using namespace std;

class Option{
public:
	Option();
	virtual ~Option();

	string save_blastout_;
	string save_muscle_;
	string save_hmmer_;
    string seq_type_;


	string blast_output_file_name_;
	string hmmbuild_output_file_name_;
	string psiblast_command_;
	string blastdbcmd_command_;
	string query_file_name_;
	string blastdb_file_name_;
	string variants_file_name_;

	string muscle_output_file_name_;
	string muscle_command_;

	string hmmer_output_file_name_;
	string hmmer_command_;

	string subject_sequences_file_name_;


	int SetOptions(int argc, char** argv);

	void PrintOptions(FILE* out);
	void PrintUsage();


};


#endif /* OPTION_H_ */
