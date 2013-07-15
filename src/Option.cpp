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
 * Name: Option.cpp
 * Author: Mingming Liu
 */

#include "Option.h"
#include <cstdio>
#include <getopt.h>
#include <stdlib.h>

char usage[] =
		"USAGE:\
<<<<<<< HEAD
		\n  hmmvar [Options]\
		\n\nOptions:\
		\n  -q <query sequence filename (required)>\
		\n  -d <database search>\
		\n  -v <variant filename (required)>\
		\n  --psiblastcmd <path to PSIBLAST>\
		\n  --musclecmd <path to MUSCLE>\
		\n  --hmmercmd <path to HMMER>\
		\n  --blastdbcmd <path to BLASTDBCMD>\
		\n  --seqtype <sequence type, prot or nucl, default nucl>\
		\n  --save_blastout <save psiblast output to filename>\
		\n  --blastout <psiblast output filename>\
		\n  Note: Please set all the paths correctly if you are not usinig the default paths.\
		\n\nExample: hmmvar -q <query filename> -v <variants filename>\n";
=======
		\n  HMMvar [Options]\
		\n\nOptions:\
		\n  -q <string>\
		\n  -d <string>\
		\n  -v <string>\
		\n  --psiblastcmd <string>\
		\n  --musclecmd <string>\
		\n  --hmmercmd <string>\
		\n  --blastdbcmd <string>\
		\n  --subject_sequence <string>\
		\n";
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d


Option::Option(){
}

Option::~Option(){

}

int Option::SetOptions(int argc, char** argv){
	int c;
	while (1)
		{
			static struct option long_options[] =
			{
					{"blastout", required_argument, 0, 1},
					{"save_blastout", required_argument, 0, 2},
					{"psiblastcmd", required_argument, 0, 3},
					{"blastdbcmd", required_argument, 0, 4},

					{"save_muscleout", required_argument,0,5},
					{"musclecmd", required_argument,0, 6},

					{"save_hmmerout",required_argument,0,7},
					{"hmmercmd",required_argument, 0, 8},

					{"subject_sequence",required_argument,0,9},

<<<<<<< HEAD
					{"seqtype",required_argument,0,10},

=======
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
			};
			/* getopt_long stores the option index here. */
			int option_index = 0;
			c = getopt_long(argc, argv, "q:d:v:",
							long_options, &option_index);
			/* Detect the end of the options. */
					if (c == -1)
						break;

					switch (c)
					{
					case 0:
						break;
					case 1:
						blast_output_file_name_ = optarg;
						break;
					case 2:
						save_blastout_ = optarg;
						break;
					case 3:
						psiblast_command_ = optarg;
						break;
					case 4:
						blastdbcmd_command_ = optarg;
						break;
					case 5:
						save_muscle_ = optarg;
						break;
					case 6:
						muscle_command_ = optarg;
						break;
					case 7:
						save_hmmer_=optarg;
						break;
					case 8:
						hmmer_command_ = optarg;
						break;
					case 9:
						subject_sequences_file_name_ = optarg;
						break;
<<<<<<< HEAD
					case 10:
					    seq_type_ = optarg;
						break;
=======
>>>>>>> c0904697a27e9226687ff8cad71154c30eece81d
					case 'q':
						query_file_name_ = optarg;
						break;
					case 'd':
						blastdb_file_name_ = optarg;
						break;
					case 'v':
						variants_file_name_ = optarg;
						break;
					case 'b':
						blast_output_file_name_ = optarg;
						break;

					default: /* '?' */
						this->PrintUsage();
						exit(EXIT_FAILURE);
					}

		}
	if (query_file_name_.empty()
				|| variants_file_name_.empty()
				) {
			this->PrintUsage();
			exit(EXIT_FAILURE);
	}

	return 1;
}

void Option::PrintUsage()
{
	fprintf(stderr, "%s", usage);
}

