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
 * Name: MultiAlign.h
 * Author: Mingming Liu
 */

#ifndef MULTIALIGN_H_
#define MULTIALIGN_H_

#include <string>

#include "Common.h"

using namespace std;

class MultiAlign{
public:
	MultiAlign(string path,string alignout_save_name);
	virtual ~MultiAlign();

//	string align_input_file_name_;
//	string align_output_file_name_;
//	string muscle_command_;
	string align_output_file_fasta_;
	string align_output_file_stockholm_;
	string output_path_;
        string alignout_save_file_name_;

//	void setParameters(string input,string output);
	int make_multi_align(string muscle_cmd,string align_input_file, string align_output_file);
	int fasta2stockholm(string oldfile);

private:
	int runMuscle(string muscle_cmd,string align_input_file);


};


#endif /* MULTIALIGN_H_ */
