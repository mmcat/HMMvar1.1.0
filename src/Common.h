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
 * Name: Common.h
 * Author: Mingming Liu
 */

#ifndef COMMON_H_
#define COMMON_H_

#define BUF_SIZE_LARGE 0x100000
#define BUF_SIZE_MED 0x10000
#define BUF_SIZE_SMALL 0x1000


#define P1 350./351
#define CONST_LOG2 0.69314718055994529
#define INTSCALE 1000.0
#define CUTOFF 0.5

#define BLASTDB "/research/mingming/ncbi-blast-2.2.26+/db/nr"
#define BLASTDBCMD "/home/mingming/bin/blastdbcmd"
#define PSIBLAST "/home/mingming/bin/psiblast"
#define MUSCLE "/home/mingming/bin/muscle"
#define HMMPATH "/home/mingming/bin/"

#define VERSION "1.1"

#include <string>
#include <cstring>
#include <cstdlib>


using namespace std;



void Log(FILE* out, string message, bool with_time);
void Log(string message, bool with_time);



#endif /* COMMON_H_ */
