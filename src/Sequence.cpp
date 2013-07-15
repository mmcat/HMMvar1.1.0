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
 * Name: Sequence.cpp
 * Author: Mingming Liu
 */

#include "Sequence.h"

#include <iostream>
#include <stdio.h>
using namespace std;

Sequence::Sequence(){

}

Sequence::~Sequence(){

}

void Sequence::Print(string& filename){
	if(filename.empty()) {fprintf(stderr,"Print sequence: No output file name provided!\n");exit(0);}


	FILE* fp = fopen(filename.c_str(),"w");
	fprintf(fp,"%s\n%s\n",this->def_.c_str(),this->seq_.c_str());
	fclose(fp);
}

bool Sequence::SetSequenceFromFastaFile(const char* id, string filename){

	FILE* fp = fopen(filename.c_str(), "r");

		if (fp == NULL) {
			Log("File open error\n", true);
			cerr << filename << endl;
			exit(1);
			return false;
		}

		id_ = id;

		char* strptr;
		char buf[BUF_SIZE_LARGE];
		while (!feof(fp)) {
			if (fgets(buf, BUF_SIZE_LARGE, fp) == NULL)
				break;

			if (buf[0] == '>') {
				// def line
				strptr = strtok(buf,"\n");
				def_.assign(strptr, 100);
			} else {
				strptr = strtok(buf, " \r\n");
				if (strptr == NULL)
					break;
				seq_.append(strptr);
			}
		}

		fclose(fp);

	return true;
}

