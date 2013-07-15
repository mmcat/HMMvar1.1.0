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
 * Name: CodeonAlphabeta.h
 * Author: Mingming Liu
 */

#ifndef CODEONALPHABETA_H_
#define CODEONALPHABETA_H_

#include <string>
#include <map>

using namespace std;



const string Tables[6] = {
		      "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
		      "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
		      "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
		      "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
		      "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
		      "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

};

const string Starts[6] = {
		      "---M---------------M---------------M----------------------------",
		      "--------------------------------MMMM---------------M------------",
		      "----------------------------------MM----------------------------",
		      "--MM---------------M------------MMMM---------------M------------",
		      "---M----------------------------MMMM---------------M------------",
		      "-----------------------------------M----------------------------",

};

const char nucleios[4]={'t','c','a','g'};
const int codeonsize=3;
const int ntsize=4;
const char gap = '-';
const string codeongap="---";

map <string,int> Codeon;
//const map <char, vector<char> > IUPAC_DNA={{'A',{'A'}},{'M',{'A','C'}},{'T',{'T'}}};
//const map <char, vector<char> > IUPAC_AA;
void set_Codeon(){
	int id=0;
	for(int i=0; i<ntsize; i++)
		for(int j=0; j<ntsize;j++)
			for(int k=0;k<ntsize;k++){
				char buf[4];
				sprintf(buf,"%c%c%c",nucleios[i],nucleios[j],nucleios[k]);
				string triple=string(buf);
				Codeon.insert(make_pair(triple,id));
				id++;
			}

}



#endif /* CODEONALPHABETA_H_ */
