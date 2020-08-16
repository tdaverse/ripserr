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
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cassert>

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
enum file_format
{
  DIPHA,
  PERSEUS
};

class DenseCubicalGrids // file_read
{
  // member vars
public:
  double threshold;
  int dim;
  int ax, ay, az;
  double dense3[512][512][512];
  file_format format;

  // constructor
  DenseCubicalGrids(const string& filename, double _threshold, file_format _format)
  {
    threshold = _threshold;
    format = _format;

    switch (format)
    {
      case DIPHA:
      {
        ifstream reading_file;

        ifstream fin( filename, ios::in | ios::binary );
        cout << filename << endl;

        int64_t d;
        fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
        assert(d == 8067171840);
        fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
        assert(d == 1);
        fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
        fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim
        dim = d;
        assert(dim == 3);
        fin.read( ( char * ) &d, sizeof( int64_t ) );
        ax = d;
        fin.read( ( char * ) &d, sizeof( int64_t ) );
        ay = d;
        fin.read( ( char * ) &d, sizeof( int64_t ) );
        az = d;
        assert(0 < ax && ax < 510 && 0 < ay && ay < 510 && 0 < az && az < 510);
        cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

        double dou;
        for(int z = 0; z < az + 2; ++z)
        {
          for (int y = 0; y < ay + 2; ++y)
          {
            for (int x = 0; x < ax + 2; ++x)
            {
              if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az)
              {
                if (!fin.eof())
                {
                  fin.read( ( char * ) &dou, sizeof( double ) );
                  dense3[x][y][z] = dou;
                }
                else
                {
                  cerr << "file endof error " << endl;
                }
              }
              else
              {
                dense3[x][y][z] = threshold;
              }
            }
          }
        }
        fin.close();
        break;
      }

      case PERSEUS:
      {
        ifstream reading_file;
        reading_file.open(filename.c_str(), ios::in);

        string reading_line_buffer;
        getline(reading_file, reading_line_buffer);
        dim = atoi(reading_line_buffer.c_str());
        getline(reading_file, reading_line_buffer);
        ax = atoi(reading_line_buffer.c_str());
        getline(reading_file, reading_line_buffer);
        ay = atoi(reading_line_buffer.c_str());
        getline(reading_file, reading_line_buffer);
        az = atoi(reading_line_buffer.c_str());
        assert(0 < ax && ax < 510 && 0 < ay && ay < 510 && 0 < az && az < 510);
        cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

        for(int z = 0; z < az + 2; ++z)
        {
          for (int y = 0; y <ay + 2; ++y)
          {
            for (int x = 0; x < ax + 2; ++x)
            {
              if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az)
              {
                if (!reading_file.eof())
                {
                  getline(reading_file, reading_line_buffer);
                  dense3[x][y][z] = atoi(reading_line_buffer.c_str());
                  if (dense3[x][y][z] == -1)
                  {
                    dense3[x][y][z] = threshold;
                  }
                }
              }
              else
              {
                dense3[x][y][z] = threshold;
              }
            }
          }
        }
        break;
      }
    }
  }

  // getters
  double getBirthday(int index, int dim)
  {
    int cx = index & 0x01ff;
    int cy = (index >> 9) & 0x01ff;
    int cz = (index >> 18) & 0x01ff;
    int cm = (index >> 27) & 0xff;

    switch (dim)
    {
      case 0:
        return dense3[cx][cy][cz];
      case 1:
        switch (cm)
        {
          case 0:
            return max(dense3[cx][cy][cz], dense3[cx + 1][cy][cz]);
          case 1:
            return max(dense3[cx][cy][cz], dense3[cx][cy + 1][cz]);
          case 2:
            return max(dense3[cx][cy][cz], dense3[cx][cy][cz + 1]);
        }
      case 2:
        switch (cm)
        {
          case 0: // x - y (fix z)
            return max(max(dense3[cx][cy][cz], dense3[cx + 1][cy][cz]), max(dense3[cx + 1][cy + 1][cz], dense3[cx][cy + 1][cz]));
          case 1: // z - x (fix y)
            return max(max(dense3[cx][cy][cz], dense3[cx][cy][cz + 1]), max(dense3[cx + 1][cy][cz + 1], dense3[cx + 1][cy][cz]));
          case 2: // y - z (fix x)
            return max(max(dense3[cx][cy][cz], dense3[cx][cy + 1][cz]), max(dense3[cx][cy + 1][cz + 1], dense3[cx][cy][cz + 1]));
        }
      case 3:
        return max(max(max(dense3[cx][cy][cz], dense3[cx + 1][cy][cz]), max(dense3[cx + 1][cy + 1][cz], dense3[cx][cy + 1][cz])), max(max(dense3[cx][cy][cz + 1], dense3[cx + 1][cy][cz + 1]), max(dense3[cx + 1][cy + 1][cz + 1], dense3[cx][cy + 1][cz + 1])));
    }
    return threshold;
  }
  void GetSimplexVertices(int index, int dim, Vertices* v)
  {
    int cx = index & 0x01ff;
    int cy = (index >> 9) & 0x01ff;
    int cz = (index >> 18) & 0x01ff;
    int cm = (index >> 27) & 0xff;

    v -> setVertices(dim ,cx, cy, cz , cm);
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
    int az = _dcg->az;
    max_of_index = 512 * 512 * (az + 2);
    int index;
    double birthday;

    for (int z = az; z > 0; --z)
    {
      for (int y = ay; y > 0; --y)
      {
        for (int x = ax; x > 0; --x)
        {
          birthday = _dcg->dense3[x][y][z];
          index = x | (y << 9) | (z << 18);
          if (birthday != _dcg->threshold)
          {
            columns_to_reduce.push_back(BirthdayIndex(birthday, index, 0));
          }
        }
      }
    }
    sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndexComparator());
  }

  // getter (length of vector)
  int size()
  {
    return columns_to_reduce.size();
  }
};

/*****union_find*****/
class UnionFind
{
  // member vars
public:
  int max_of_index;
  vector<int> parent;
  vector<double> birthtime;
  vector<double> time_max;
  DenseCubicalGrids* dcg;

  // constructor
  UnionFind(int moi, DenseCubicalGrids* _dcg) : parent(moi), birthtime(moi), time_max(moi) // Thie "n" is the number of cubes.
  {
    dcg = _dcg;
    max_of_index = moi;

    for (int i = 0; i < moi; ++i)
    {
      parent[i] = i;
      birthtime[i] = dcg->getBirthday(i, 0);
      time_max[i] = dcg->getBirthday(i, 0);
    }
  }

  // member methods
  int find(int x) // Thie "x" is Index.
  {
    int y = x, z = parent[y];
    while (z != y)
    {
      y = z;
      z = parent[y];
    }
    y = parent[x];
    while (z != y)
    {
      parent[x] = z;
      x = y;
      y = parent[x];
    }
    return z;
  }

  void link(int x, int y)
  {
    x = find(x);
    y = find(y);
    if (x == y) return;
    if (birthtime[x] > birthtime[y])
    {
      parent[x] = y;
      birthtime[y] = min(birthtime[x], birthtime[y]);
      time_max[y] = max(time_max[x], time_max[y]);
    }
    else if (birthtime[x] < birthtime[y])
    {
      parent[y] = x;
      birthtime[x] = min(birthtime[x], birthtime[y]);
      time_max[x] = max(time_max[x], time_max[y]);
    }
    else //birthtime[x] == birthtime[y]
    {
      parent[x] = y;
      time_max[y] = max(time_max[x], time_max[y]);
    }
  }
};

/*****simplex_coboundary_enumerator*****/
class SimplexCoboundaryEnumerator
{
  // member vars
public:
  BirthdayIndex simplex;
  DenseCubicalGrids* dcg;
  Vertices* vtx;
  double birthtime;
  int ax, ay, az;
  int cx, cy, cz;
  int count;
  BirthdayIndex nextCoface;
  double threshold;

  // constructor
  SimplexCoboundaryEnumerator()
  {
    vtx = new Vertices();
    nextCoface = BirthdayIndex(0, -1, 1);
  }

  // setter
  void setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg)
  {
    simplex = _s;
    dcg = _dcg;
    _dcg->GetSimplexVertices(simplex.index, simplex.dim, vtx);
    birthtime = simplex.birthday;
    ax = _dcg->ax;
    ay = _dcg->ay;
    az = _dcg->az;

    threshold = _dcg->threshold;
    count = 0;
  }

  // member methods
  bool hasNextCoface()
  {
    int index = 0;
    double birthday = 0;
    cx = vtx->ox;
    cy = vtx->oy;
    cz = vtx->oz;
    switch (vtx->dim)
    {
      case 0: // dim0
        for (int i = count; i < 6; ++i)
        {
          switch (i)
          {
            case 0:
              index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
              birthday = max(birthtime, dcg->dense3[cx][cy][cz + 1]);
              break;

            case 1:
              index = (2 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
              birthday = max(birthtime, dcg->dense3[cx][cy][cz - 1]);
              break;

            case 2:
              index = (1 << 27) | (cz << 18) | (cy << 9) | cx;
              birthday = max(birthtime, dcg -> dense3[cx][cy + 1][cz]);
              break;

            case 3:
              index = (1 << 27) | (cz << 18) | ((cy - 1) << 9) | cx;
              birthday = max(birthtime, dcg -> dense3[cx][cy - 1][cz]);
              break;

            case 4:
              index = (0 << 27) | (cz << 18) | (cy << 9) | cx;
              birthday = max(birthtime, dcg -> dense3[cx + 1][cy][cz]);
              break;

            case 5:
              index = (0 << 27) | (cz << 18) | (cy << 9) | (cx - 1);
              birthday = max(birthtime, dcg -> dense3[cx - 1][cy][cz]);
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

      case 1: // dim1
        switch (vtx->type)
        {
          case 0: // dim1 type0 (x-axis -> )
            for (int i = count; i < 4; ++i)
            {
              switch (i)
              {
                case 0:
                  index = (1 << 27) | (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx + 1][cy][cz + 1]});
                  break;

                case 1:
                  index = (1 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx + 1][cy][cz - 1]});
                  break;

                case 2:
                  index = (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx + 1][cy + 1][cz]});
                  break;

                case 3:
                  index = (cz << 18) | ((cy - 1) << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx + 1][cy - 1][cz]});
                  break;
              }

              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            return false;

          case 1: // dim1 type1 (y-axis -> )
            for (int i = count; i < 4; ++i)
            {
              switch(i)
              {
                case 0:
                  index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx][cy + 1][cz + 1]});
                  break;

                case 1:
                  index = (2 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx][cy + 1][cz - 1]});
                  break;

                case 2:
                  index = (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz]});
                  break;

                case 3:
                  index = (cz << 18) | (cy << 9) | (cx - 1);
                  birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz]});
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            return false;

          case 2: // dim1 type2 (z-axis -> )
            for(int i = count; i < 4; ++i)
            {
              switch (i)
              {
                case 0:
                  index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx][cy + 1][cz + 1]});
                  break;

                case 1:
                  index = (2 << 27) | (cz << 18) | ((cy - 1) << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx][cy - 1][cz + 1]});
                  break;

                case 2:
                  index = (1 << 27) | (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy][cz + 1]});
                  break;

                case 3:
                  index = (1 << 27) | (cz << 18) | (cy << 9) | (cx - 1);
                  birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy][cz + 1]});
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex(birthday, index, 2);
                return true;
              }
            }
            return false;
        }
        return false;

      default: // dim2
        switch (vtx->type)
        {
          case 0: // dim2 type0 (fix z)
            for (int i = count; i < 2; ++i)
            {
              switch (i)
              {
                case 0: // upper
                  index = (cz << 18) | (cy << 9) | cx;
                  birthday = max(birthtime, max(max(dcg->dense3[cx][cy][cz + 1], dcg->dense3[cx + 1][cy][cz + 1]), max(dcg->dense3[cx][cy + 1][cz + 1], dcg->dense3[cx + 1][cy + 1][cz + 1])));
                  break;

                case 1: // lower
                  index = ((cz - 1) << 18) | (cy << 9) | cx;
                  birthday = max(birthtime, max(max(dcg->dense3[cx][cy][cz - 1], dcg->dense3[cx + 1][cy][cz - 1]), max(dcg->dense3[cx][cy + 1][cz - 1], dcg->dense3[cx + 1][cy + 1][cz - 1])));
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex(birthday, index, 3);
                return true;
              }
            }
            return false;

          case 1: // dim2 type1 (fix y)
            for (int i = count; i < 2; ++i)
            {
              switch (i)
              {
                case 0: // left
                  index = (cz << 18) | (cy << 9) | cx;
                  birthday = max(birthtime, max(max(dcg->dense3[cx][cy + 1][cz], dcg->dense3[cx + 1][cy + 1][cz]), max(dcg->dense3[cx][cy + 1][cz + 1], dcg->dense3[cx + 1][cy + 1][cz + 1])));
                  break;

                case 1: //right
                  index = (cz << 18) | ((cy - 1) << 9) | cx;
                  birthday = max(birthtime, max(max(dcg->dense3[cx][cy - 1][cz], dcg->dense3[cx + 1][cy - 1][cz]), max(dcg->dense3[cx][cy - 1][cz + 1], dcg->dense3[cx + 1][cy - 1][cz + 1])));
                  break;

              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex(birthday, index, 3);
                return true;
              }
            }
            return false;

          case 2: // dim2 type2 (fix x)
            for (int i = count; i < 2; ++i)
            {
              switch (i)
              {
                case 0: // left
                  index = (cz << 18) | (cy << 9) | cx;
                  birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz], dcg -> dense3[cx + 1][cy][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
                  break;

                case 1: //right
                  index = (cz << 18) | (cy << 9) | (cx - 1);
                  birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz], dcg -> dense3[cx - 1][cy][cz + 1],dcg -> dense3[cx - 1][cy + 1][cz + 1]});
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex(birthday, index, 3);
                return true;
              }
            }
            return false;
        }
        return false;
    }
  }

  BirthdayIndex getNextCoface()
  {
    return nextCoface;
  }
};

/*****joint_pairs*****/

/*****compute_pairs*****/

/*****cubicalripser_3dim*****/
//PLACEHOLDER - WILL REPLACE
int main()
{
  return 0;
}
