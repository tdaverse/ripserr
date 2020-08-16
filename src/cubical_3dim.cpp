/*
 This file is an altered form of the Cubical Ripser software created by
 Takeki Sudo and Kazushi Ahara. Details of the original software are below the
 dashed line.
 -Raoul Wadhwa
 -------------------------------------------------------------------------------
 Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
 This file is part of CubicalRipser_3dim.
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
#include <cstdint>

using namespace std;

/*****birthday_index*****/
class BirthdayIndex
{
  // member vars
public:
  double birthday;
  int index;
  int dim;

  // default empty constructor
  BirthdayIndex()
  {
    birthday = 0;
    index = -1;
    dim = 1;
  }

  // member setter constructor
  BirthdayIndex(double _b, int _i, int _d)
  {
    birthday = _b;
    index = _i;
    dim = _d;
  }

  // clone/copy constructor
  BirthdayIndex(const BirthdayIndex& b)
  {
    birthday = b.birthday;
    index = b.index;
    dim = b.dim;
  }

  // clone/copy method
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

  // output methods
  void print()
  {
    cout << "(dob:" << birthday << "," << index << ")" << endl;
  }
  void VertexPrint()
  {
    int px = index & 0x01ff;
    int py = (index >> 9) & 0x01ff;
    int pz = (index >> 18) & 0x01ff;
    int pm = (index >> 27) & 0xff;

    cout << "birthday : (m, z, y, x) = " << birthday << " : (" << pm << ", " << pz << ", " << py << ", " << px << ")" << endl;
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

/*****coeff*****/
class Coeff
{
  // member vars
public:
  int cx, cy, cz, cm;

  // constructor
  Coeff()
  {
    cx = 0;
    cy = 0;
    cz = 0;
    cm = 0;
  }

  // setters
  void setXYZ(int _cx, int _cy, int _cz)
  {
    cx = _cx;
    cy = _cy;
    cz = _cz;
    cm = 0;
  }
  void setXYZM(int _cx, int _cy, int _cz, int _cm)
  {
    cx = _cx;
    cy = _cy;
    cz = _cz;
    cm = _cm;
  }
  void setIndex(int index)
  {
    cx = index & 0x01ff;
    cy = (index >> 9) & 0x01ff;
    cz = (index >> 18) & 0x01ff;
    cm = (index >> 27) & 0xff;
  }

  // getter
  int getIndex()
  {
    return cx | cy << 9 | cz << 18 | cm << 27;
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

/*****vertices*****/
class Vertices
{
  // member vars
public:
  Coeff** vertex;
  int dim;
  int ox, oy, oz;
  int type;

  // constructor
  Vertices()
  {
    dim = 0;
    vertex = new Coeff*[8];
    for(int d = 0; d < 8; ++d)
    {
      vertex[d] = new Coeff();
    }
  }

  // setter
  void setVertices(int _dim, int _ox, int _oy, int _oz, int _om) // 0 cell
  {
    dim = _dim;
    ox = _ox;
    oy = _oy;
    oz = _oz;
    type = _om;

    if (dim == 0)
    {
      vertex[0]->setXYZ(_ox, _oy, _oz);
    }
    else if (dim == 1)
    {
      switch(_om)
      {
        case 0:
          vertex[0] -> setXYZ(_ox, _oy, _oz);
          vertex[1] -> setXYZ(_ox + 1, _oy, _oz);
          break;

        case 1:
          vertex[0] -> setXYZ(_ox, _oy, _oz);
          vertex[1] -> setXYZ(_ox, _oy + 1, _oz);
          break;

        default:
          vertex[0] -> setXYZ(_ox, _oy, _oz);
          vertex[1] -> setXYZ(_ox, _oy, _oz + 1);
          break;
      }
    }
    else if (dim == 2)
    {
      switch (_om)
      {
        case 0: // x - y
          vertex[0] -> setXYZ(_ox, _oy, _oz);
          vertex[1] -> setXYZ(_ox + 1, _oy, _oz);
          vertex[2] -> setXYZ(_ox + 1, _oy + 1, _oz);
          vertex[3] -> setXYZ(_ox, _oy + 1, _oz);
          break;

        case 1: // z - x
          vertex[0] -> setXYZ(_ox, _oy, _oz);
          vertex[1] -> setXYZ(_ox, _oy, _oz + 1);
          vertex[2] -> setXYZ(_ox + 1, _oy, _oz + 1);
          vertex[3] -> setXYZ(_ox + 1, _oy, _oz);
          break;

        default: // y - z
          vertex[0] -> setXYZ(_ox, _oy, _oz);
          vertex[1] -> setXYZ(_ox, _oy + 1, _oz);
          vertex[2] -> setXYZ(_ox, _oy + 1, _oz + 1);
          vertex[3] -> setXYZ(_ox, _oy, _oz + 1);
          break;
      }
    }
    else if (dim == 3) // cube
    {
      vertex[0] -> setXYZ(_ox, _oy, _oz);
      vertex[1] -> setXYZ(_ox + 1, _oy, _oz);
      vertex[2] -> setXYZ(_ox + 1, _oy + 1, _oz);
      vertex[3] -> setXYZ(_ox, _oy + 1, _oz);
      vertex[4] -> setXYZ(_ox, _oy, _oz + 1);
      vertex[5] -> setXYZ(_ox + 1, _oy, _oz + 1);
      vertex[6] -> setXYZ(_ox + 1, _oy + 1, _oz + 1);
      vertex[7] -> setXYZ(_ox, _oy + 1, _oz + 1);
    }
  }
};

/*****dense_cubical_grids*****/

/*****columns_to_reduce*****/

/*****union_find*****/

/*****simplex_coboundary_enumerator*****/

/*****joint_pairs*****/

/*****compute_pairs*****/

/*****cubicalripser_3dim*****/
//PLACEHOLDER - WILL REPLACE
int main()
{
  return 0;
}
