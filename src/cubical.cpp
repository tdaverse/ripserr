/*
 This file is an altered form of the Cubical Ripser software created by
 Takeki Sudo and Kazushi Ahara. Details of the original software are below the
 dashed line.
 -Raoul Wadhwa
 -------------------------------------------------------------------------------
 Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
 This file is part of CubicalRipser_2dim.
 CubicalRipser: C++ system for computation of Cubical persistence pairs
 Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
 CubicalRipser is free software: you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the
 Free Software Foundation, either version 3 of the License, or (at your option)
 any later version.
 CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips
 persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
 We rearrange his codes of Ripser and add some new ideas for optimization on it
 and modify it for calculation of a Cubical filtration.
 This part of CubicalRiper is a calculator of cubical persistence pairs for
 2 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.
 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
 You should have received a copy of the GNU Lesser General Public License along
 with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <Rcpp.h>

using namespace std;

/*****birthday_index*****/
class BirthdayIndex
{
  // member vars
public:
  double birthday;
  int index;
  int dim;
  
  // constructors
  BirthdayIndex(double b, int i, int d) : birthday(b), index(i), dim(d) {}
  BirthdayIndex() : BirthdayIndex(0, -1, 1) {}
  BirthdayIndex(const BirthdayIndex& b) : BirthdayIndex(b.birthday, b.index, b.dim) {}
  
  // copy
  void copyBirthdayIndex(BirthdayIndex v)
  {
    birthday = v.birthday;
    index = v.index;
    dim = v.dim;
  }
  
  // getters
  double getBirthday() { return birthday; }
  int getIndex() { return index; }
  int getDimension() { return dim; }
};

// utility method
bool cmpBday(const BirthdayIndex& o1, const BirthdayIndex& o2)
{
  return (o1.birthday == o2.birthday ? o1.index < o2.index : o1.birthday < o2.birthday);
}

struct BirthdayIndexComparator
{
  bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const
  {
    return cmpBday(o1, o2);
  }
};

struct BirthdayIndexInverseComparator
{
  bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const
  {
    return !cmpBday(o1, o2);
  }
};

/*****write_pairs*****/
class WritePairs
{
  // member vars
public:
  int64_t dim;
  double birth;
  double death;
  
  // constructor
  WritePairs(int64_t setDim, double setBirth, double setDeath) : dim(setDim), birth(setBirth), death(setDeath) {}
  
  // getters
  int64_t getDimension() { return dim; }
  double getBirth() { return birth; }
  double getDeath() { return death; }
};
