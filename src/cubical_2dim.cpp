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
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdint>

using namespace std;

/*****birthday_index*****/
class BirthdayIndex
{
  //member vars
public:
  double birthday;
  int index;
  int dim;

  // default constructor
  BirthdayIndex()
  {
    birthday = 0;
    index = -1;
    dim = 1;
  }

  // individual params constructor
  BirthdayIndex(double _b, int _i, int _d)
  {
    birthday = _b;
    index = _i;
    dim = _d;
  }

  // copy/clone constructor
  BirthdayIndex(const BirthdayIndex& b)
  {
    birthday = b.birthday;
    index = b.index;
    dim = b.dim;
  }

  // copy method
  void copyBirthdayIndex(BirthdayIndex v)
  {
    birthday = v.birthday;
    index = v.index;
    dim = v.dim;
  }

  // getters
  double getBirthday()
  {
    return birthday;
  }
  long getIndex()
  {
    return index;
  }
  int getDimension()
  {
    return dim;
  }

  // print member vals to console
  void print()
  {
    cout << "(dob:" << birthday << "," << index << ")" << endl;
  }
};

struct BirthdayIndexComparator
{
  bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const
  {
    if (o1.birthday == o2.birthday)
    {
      if (o1.index < o2.index)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      if (o1.birthday > o2.birthday)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  }
};

struct BirthdayIndexInverseComparator
{
  bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const
  {
    if (o1.birthday == o2.birthday)
    {
      if (o1.index < o2.index)
      {
        return false;
      }
      else
      {
        return true;
      }
    }
    else
    {
      if (o1.birthday > o2.birthday)
      {
        return false;
      }
      else
      {
        return true;
      }
    }
  }
};

/*****dense_cubical_grids*****/
enum file_format
{
  DIPHA,
  PERSEUS
};

class DenseCubicalGrids // file_read
{
public:
  double threshold;
  int dim;
  int ax, ay;
  double dense2[2048][1024];
  file_format format;

  // constructor (w/ file read)
  DenseCubicalGrids(const string& filename, double _threshold, file_format _format)
  {
    threshold = _threshold;
    format = _format;

    if (format == DIPHA) // ???.complex, DIPHA format
    {
      ifstream reading_file;

      ifstream fin(filename, ios::in | ios::binary);
      cout << filename << endl;

      int64_t d;
      fin.read((char*) &d, sizeof(int64_t)); // magic number
      assert(d == 8067171840);
      fin.read((char*) &d, sizeof(int64_t)); // type number
      assert(d == 1);
      fin.read((char*) &d, sizeof(int64_t)); //data num
      fin.read((char*) &d, sizeof(int64_t)); // dim
      dim = d;
      assert(dim == 2);
      fin.read((char*) &d, sizeof(int64_t));
      ax = d;
      fin.read((char*) &d, sizeof(int64_t));
      ay = d;
      assert(0 < ax && ax < 2000 && 0 < ay && ay < 1000);
      cout << "ax : ay = " << ax << " : " << ay << endl;

      double dou;
      for (int y = 0; y < ay + 2; ++y)
      {
        for (int x = 0; x < ax + 2; ++x)
        {
          if (0 < x && x <= ax && 0 < y && y <= ay)
          {
            if (!fin.eof())
            {
              fin.read((char*) &dou, sizeof(double));
              dense2[x][y] = dou;
            }
            else
            {
              cout << "file endof error " << endl;
            }
          }
          else
          {
            dense2[x][y] = threshold;
          }
        }
      }
      fin.close();
    }
    else if(format == PERSEUS)// PERSEUS format
    {
      ifstream reading_file;
      reading_file.open(filename.c_str(), ios::in);
      cout << filename << endl;

      string reading_line_buffer;
      getline(reading_file, reading_line_buffer);
      dim = atoi(reading_line_buffer.c_str());
      getline(reading_file, reading_line_buffer);
      ax = atoi(reading_line_buffer.c_str());
      getline(reading_file, reading_line_buffer);
      ay = atoi(reading_line_buffer.c_str());
      assert(0 < ax && ax < 2000 && 0 < ay && ay < 1000);
      cout << "ax : ay = " << ax << " : " << ay << endl;

      for (int y = 0; y <ay + 2; ++y)
      {
        for (int x = 0; x < ax + 2; ++x)
        {
          if (0 < x && x <= ax && 0 < y && y <= ay)
          {
            if (!reading_file.eof())
            {
              getline(reading_file, reading_line_buffer);
              dense2[x][y] = atoi(reading_line_buffer.c_str());
              if (dense2[x][y] == -1)
              {
                dense2[x][y] = threshold;
              }
            }
          }
          else
          {
            dense2[x][y] = threshold;
          }
        }
      }
    }
  }

  // getter
  double getBirthday(int index, int dim)
  {
    int cx = index & 0x07ff;
    int cy = (index >> 11) & 0x03ff;
    int cm = (index >> 21) & 0xff;

    switch (dim)
    {
      case 0:
        return dense2[cx][cy];
      case 1:
        switch (cm)
        {
          case 0:
            return max(dense2[cx][cy], dense2[cx + 1][cy]);
          default:
            return max(dense2[cx][cy], dense2[cx][cy + 1]);
        }
      case 2:
        return max(max(dense2[cx][cy], dense2[cx + 1][cy]), max(dense2[cx][cy + 1], dense2[cx + 1][cy + 1]));
    }
    return threshold;
  }
};

/*****write_pairs*****/
class WritePairs
{
public:
  int64_t dim;
  double birth;
  double death;

  // constructor
  WritePairs(int64_t _dim, double _birth, double _death)
  {
    dim = _dim;
    birth = _birth;
    death = _death;
  }

  // getters
  int64_t getDimension()
  {
    return dim;
  }
  double getBirth()
  {
    return birth;
  }
  double getDeath()
  {
    return death;
  }
};

/*****columns_to_reduce*****/
class ColumnsToReduce
{
  // member vars
public:
  vector<BirthdayIndex> columns_to_reduce;
  int dim;
  int max_of_index;

  // constructor
  ColumnsToReduce(DenseCubicalGrids* _dcg)
  {
    dim = 0;
    int ax = _dcg->ax;
    int ay = _dcg->ay;
    max_of_index = 2048 * (ay + 2);
    int index;
    double birthday;
    for (int y = ay; y > 0; --y)
    {
      for (int x = ax; x > 0; --x)
      {
        birthday = _dcg->dense2[x][y];
        index = x | (y << 11);
        if (birthday != _dcg -> threshold)
        {
          columns_to_reduce.push_back(BirthdayIndex(birthday, index, 0));
        }
      }
    }
    sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndexComparator());
  }

  // getter (length of member vector)
  int size()
  {
    return columns_to_reduce.size();
  }
};

/*****simplex_coboundary_enumerator*****/
class SimplexCoboundaryEnumerator
{
  // member vars
public:
  BirthdayIndex simplex;
  DenseCubicalGrids* dcg;
  int dim;
  double birthtime;
  int ax, ay;
  int cx, cy, cm;
  int count;
  BirthdayIndex nextCoface;
  double threshold;

  // constructor
  SimplexCoboundaryEnumerator()
  {
    nextCoface = BirthdayIndex(0, -1, 1);
  }

  // member methods
  void setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg)
  {
    simplex = _s;
    dcg = _dcg;
    dim = simplex.dim;
    birthtime = simplex.birthday;
    ax = _dcg->ax;
    ay = _dcg->ay;

    cx = (simplex.index) & 0x07ff;
    cy = (simplex.index >> 11) & 0x03ff;
    cm = (simplex.index >> 21) & 0xff;

    threshold = _dcg->threshold;
    count = 0;
  }
  bool hasNextCoface()
  {
    int index = 0;
    double birthday = 0;
    switch (dim)
    {
      case 0:
        for (int i = count; i < 4; i++)
        {
          switch (i)
          {
            case 0: // y+
              index = (1 << 21) | ((cy) << 11) | (cx);
              birthday = max(birthtime, dcg->dense2[cx  ][cy+1]);
              break;
            case 1: // y-
              index = (1 << 21) | ((cy-1) << 11) | (cx);
              birthday = max(birthtime, dcg->dense2[cx  ][cy-1]);
              break;
            case 2: // x+
              index = (0 << 21) | ((cy) << 11) | (cx);
              birthday = max(birthtime, dcg->dense2[cx+1][cy  ]);
              break;
            case 3: // x-
              index = (0 << 21) | ((cy) << 11) | (cx-1);
              birthday = max(birthtime, dcg->dense2[cx-1][cy  ]);
              break;
          }

          if (birthday != threshold)
          {
            count = i + 1;
            nextCoface = BirthdayIndex(birthday, index, 1);
            return true;
          }
        }
        return false;
      case 1:
        switch (cm)
        {
          case 0:
            if (count == 0) // upper
            {
              count++;
              index = ((cy) << 11) | cx;
              birthday = max(max(birthtime, dcg->dense2[cx][cy + 1]), dcg->dense2[cx + 1][cy + 1]);
              if (birthday != threshold)
              {
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            if (count == 1) // lower
            {
              count++;
              index = ((cy - 1) << 11) | cx;
              birthday = max(max(birthtime, dcg->dense2[cx][cy - 1]), dcg->dense2[cx + 1][cy - 1]);
              if (birthday != threshold)
              {
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            return false;
          case 1:
            if (count == 0) // right
            {
              count ++;
              index = ((cy) << 11) | cx;
              birthday = max(max(birthtime, dcg->dense2[cx + 1][cy]), dcg->dense2[cx + 1][cy + 1]);
              if (birthday != threshold)
              {
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            if (count == 1) //left
            {
              count++;
              index = ((cy) << 11) | (cx - 1);
              birthday = max(max(birthtime, dcg->dense2[cx - 1][cy]), dcg->dense2[cx - 1][cy + 1]);
              if (birthday != threshold)
              {
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            return false;
        }
    }
    return false;
  }

  // getter
  BirthdayIndex getNextCoface()
  {
    return nextCoface;
  }
};

// placeholder for compilation to work - REMOVE once cubical_2dim is fully ported
int main()
{
  return 0;
}
