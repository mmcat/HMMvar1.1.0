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
 * Name: Sequence.h
 * Author: Mingming Liu
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include "Common.h"

using namespace std;

class Sequence {
public:
	Sequence();
	virtual ~Sequence();

	bool SetSequence();
	bool SetSequence(char* id, string blast_db_file, string blastdbcmd);
	bool SetSequenceFromFastaFile(const char* id, string filename);
	void Print(string& filename);

	string id_;
	string seq_;
	string def_;
//	int cluster_id_;
	double e_value_;
	double bit_score_;
	double weight_;

	static char buf_[BUF_SIZE_LARGE];
};


#endif /* SEQUENCE_H_ */
