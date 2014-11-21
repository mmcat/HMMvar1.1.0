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
 * Name: variant.h
 * Author: Mingming Liu
 */

#ifndef VARIANT_H_
#define VARIANT_H_

#include <string>

using namespace std;

class variant{
public:
	variant();
	virtual ~variant();
	string varinat_str_;
	string variant_id_;
	string set_id_;
	string commends_;
	string type;
	string wild;
	string mutant;
//	double wt_proba_;
//	double mt_proba_;
//	double odds_;

	int pos;


};


#endif /* VARIANT_H_ */
