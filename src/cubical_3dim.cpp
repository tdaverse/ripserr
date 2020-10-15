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
#include <cassert>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>
#include <queue>
#include <Rcpp.h>

using namespace std;

/*****birthday_index*****/
class BirthdayIndex3
{
public:
  double birthday;
  int index;
  int dim;
  
  BirthdayIndex3() : birthday(0), index(-1), dim(1) {};
  BirthdayIndex3(double _b, int _i, int _d) : birthday(_b), index(_i), dim(_d) {};
  BirthdayIndex3(const BirthdayIndex3& b) : birthday(b.birthday), index(b.index), dim(b.dim) {};
  
  void copyBirthdayIndex3(BirthdayIndex3 v)
  {
    birthday = v.birthday;
    index = v.index;
    dim = v.dim;
  }
  
  double getBirthday() { return birthday; }
  long getIndex() { return index; }
  int getDimension() { return dim; }
};

bool bdayCmp(const BirthdayIndex3& o1, const BirthdayIndex3& o2)
{
  return (o1.birthday == o2.birthday) ? (o1.index < o2.index) : (o1.birthday > o2.birthday);
}

struct BirthdayIndex3Comparator
{
  bool operator()(const BirthdayIndex3& o1, const BirthdayIndex3& o2) const
  {
    return bdayCmp(o1, o2);
  }
};

struct BirthdayIndex3InverseComparator
{
  bool operator()(const BirthdayIndex3& o1, const BirthdayIndex3& o2) const
  {
    return !bdayCmp(o1, o2);
  }
};

/*****coeff*****/
class Coeff
{
public:
  int cx, cy, cz, cm;
  
  Coeff() : cx(0), cy(0), cz(0), cm(0) {};
  
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
  
  int getIndex()
  {
    return cx | cy << 9 | cz << 18 | cm << 27;
  }
};

/*****write_pairs*****/
class WritePairs3
{
public:
  int64_t dim;
  double birth;
  double death;
  
  WritePairs3(int64_t _dim, double _birth, double _death) : dim(_dim), birth(_birth), death(_death) {};
  
  int64_t getDimension() { return dim; }
  double getBirth() { return birth; }
  double getDeath() { return death; }
};

/*****vertices*****/
class Vertices
{
public:	
  Coeff* vertex[8];
  int dim; 
  int ox, oy, oz;
  int type;
  
  Vertices() : dim(0)
  {
    for (int d = 0; d < 8; ++d)
      vertex[d] = new Coeff();
  }
  
  // free pointers
  ~Vertices()
  {
    for (int d = 0; d < 8; d++)
      delete vertex[d];
  }
  
  void setVertices(int _dim, int _ox, int _oy, int _oz, int _om) // 0 cell
  {
    dim = _dim;
    ox = _ox;
    oy = _oy;
    oz = _oz;
    type = _om;
    
    if (dim == 0)
    {
      vertex[0] -> setXYZ(_ox, _oy, _oz);
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
class DenseCubicalGrids3
{
public:
  double threshold;
  int dim;
  int ax, ay, az;
  double dense3[512][512][512];
  
  DenseCubicalGrids3(const Rcpp::NumericVector& image, double _threshold, int nx, int ny, int nz) : threshold(_threshold), ax(nx), ay(ny), az(nz)
  {
    dim = 3;
    
    // set everything to threshold
    for (int i = 0; i < 512; i++)
      for (int j = 0; j < 512; j++)
        for (int k = 0; k < 512; k++)
          dense3[i][j][k] = threshold;
    
    // set values based on image
    for (int i = 0; i < ax * ay * az; i++)
      dense3[i % ax + 1][i / ax % ay + 1][i / (ax * ay) % az + 1] = image(i);
  }
  
  double getBirthday(int index, int dim)
  {
    int cx = index & 0x01ff;
    int cy = (index >> 9) & 0x01ff;
    int cz = (index >> 18) & 0x01ff;
    int cm = (index >> 27) & 0xff;
    
    switch(dim)
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
            return max({dense3[cx][cy][cz], dense3[cx + 1][cy][cz], 
                        dense3[cx + 1][cy + 1][cz], dense3[cx][cy + 1][cz]});
          case 1: // z - x (fix y)
            return max({dense3[cx][cy][cz], dense3[cx][cy][cz + 1], 
                        dense3[cx + 1][cy][cz + 1], dense3[cx + 1][cy][cz]});
          case 2: // y - z (fix x)
            return max({dense3[cx][cy][cz], dense3[cx][cy + 1][cz], 
                        dense3[cx][cy + 1][cz + 1], dense3[cx][cy][cz + 1]});
        }
      case 3:
        return max({dense3[cx][cy][cz], dense3[cx + 1][cy][cz], 
                    dense3[cx + 1][cy + 1][cz], dense3[cx][cy + 1][cz],
                    dense3[cx][cy][cz + 1], dense3[cx + 1][cy][cz + 1],
                    dense3[cx + 1][cy + 1][cz + 1], dense3[cx][cy + 1][cz + 1]});
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

/*****union_find*****/
class UnionFind3
{
public:
  int max_of_index;
  vector<int> parent;
  vector<double> birthtime;
  vector<double> time_max;
  DenseCubicalGrids3* dcg;
  
  UnionFind3(int moi, DenseCubicalGrids3* _dcg) : parent(moi), birthtime(moi), time_max(moi)
  {
    dcg = _dcg;
    max_of_index = moi;
    
    for(int i = 0; i < moi; ++i){
      parent[i] = i;
      birthtime[i] = dcg -> getBirthday(i, 0);
      time_max[i] = dcg -> getBirthday(i, 0);
    }
  }
  
  int find(int x) // Thie "x" is Index.
  {
    int y = x,
        z = parent[y];
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
    else if(birthtime[x] < birthtime[y])
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

/*****columns_to_reduce*****/
class ColumnsToReduce3
{
public:
  vector<BirthdayIndex3> columns_to_reduce;
  int dim;
  int max_of_index;
  
  ColumnsToReduce3(DenseCubicalGrids3* _dcg)
  { 
    dim = 0;
    int ax = _dcg -> ax;
    int ay = _dcg -> ay;
    int az = _dcg -> az;
    max_of_index = 512 * 512 * (az + 2);
    int index;
    double birthday;
    
    for(int z = az; z > 0; --z)
      for (int y = ay; y > 0; --y)
        for (int x = ax; x > 0; --x)
        {
          birthday = _dcg -> dense3[x][y][z];
          index = x | (y << 9) | (z << 18);
          if (birthday != _dcg -> threshold) columns_to_reduce.push_back(BirthdayIndex3(birthday, index, 0));
        }
        sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndex3Comparator());
  }
  
  int size() { return columns_to_reduce.size(); }
};

/*****simplex_coboundary_estimator*****/
class SimplexCoboundaryEnumerator3
{
public:
  BirthdayIndex3 simplex;
  DenseCubicalGrids3* dcg;
  Vertices* vtx;
  double birthtime;
  int ax, ay, az;
  int cx, cy, cz;
  int count;
  BirthdayIndex3 nextCoface;
  double threshold;
  
  SimplexCoboundaryEnumerator3()
  {
    vtx = new Vertices();
    nextCoface = BirthdayIndex3(0, -1, 1);
  }
  
  // free pointers
  ~SimplexCoboundaryEnumerator3()
  {
    delete vtx;
  }
  
  void setSimplexCoboundaryEnumerator3(BirthdayIndex3 _s, DenseCubicalGrids3* _dcg)
  {
    simplex = _s;
    dcg = _dcg;
    _dcg -> GetSimplexVertices(simplex.index, simplex.dim, vtx);
    birthtime = simplex.birthday;
    ax = _dcg -> ax;
    ay = _dcg -> ay;
    az = _dcg -> az;
    
    threshold = _dcg -> threshold;
    count = 0;
  }
  
  bool hasNextCoface()
  {
    int index = 0;
    double birthday = 0;
    cx = vtx -> ox;
    cy = vtx -> oy;
    cz = vtx -> oz;
    switch (vtx->dim) {
    case 0: // dim0
      for (int i = count; i < 6; ++i) {
        switch (i){
        case 0:
          index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
          birthday = max(birthtime, dcg -> dense3[cx][cy][cz + 1]);
          break;
          
        case 1:
          index = (2 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
          birthday = max(birthtime, dcg -> dense3[cx][cy][cz - 1]);
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
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex3(birthday, index, 1);
          return true;
        }
      }
      return false;
      
    case 1: // dim1
      switch (vtx->type) {
      case 0: // dim1 type0 (x-axis -> )
        for(int i = count; i < 4; ++i){
          switch(i){
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
          
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex3(birthday, index, 2);
            return true;
          }
        }
        return false;
        
      case 1: // dim1 type1 (y-axis -> )
        for(int i = count; i < 4; ++i){
          switch(i){
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
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex3(birthday, index, 2);
            return true;
          }
        }
        return false;
        
      case 2: // dim1 type2 (z-axis -> )
        for(int i = count; i < 4; ++i){
          switch(i){
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
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex3(birthday, index, 2);
            return true;
          }
        }
        return false;
      }
      return false;
      
    default: // dim2
      switch (vtx->type) {
      case 0: // dim2 type0 (fix z)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // upper
            index = (cz << 18) | (cy << 9) | cx;
            birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx + 1][cy][cz + 1], 
                           dcg -> dense3[cx][cy + 1][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
            break;
            
          case 1: // lower
            index = ((cz - 1) << 18) | (cy << 9) | cx;
            birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx + 1][cy][cz - 1], 
                           dcg -> dense3[cx][cy + 1][cz - 1],dcg -> dense3[cx + 1][cy + 1][cz - 1]});
            break;
            
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex3(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 1: // dim2 type1 (fix y)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // left
            index = (cz << 18) | (cy << 9) | cx;
            birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx + 1][cy + 1][cz], 
                           dcg -> dense3[cx][cy + 1][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
            break;
            
          case 1: //right
            index = (cz << 18) | ((cy - 1) << 9) | cx;
            birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx + 1][cy - 1][cz], 
                           dcg -> dense3[cx][cy - 1][cz + 1],dcg -> dense3[cx + 1][cy - 1][cz + 1]});
            break;
            
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex3(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 2: // dim2 type2 (fix x)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // left
            index = (cz << 18) | (cy << 9) | cx;
            birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz], 
                           dcg -> dense3[cx + 1][cy][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
            break;
            
          case 1: //right
            index = (cz << 18) | (cy << 9) | (cx - 1);
            birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz], 
                           dcg -> dense3[cx - 1][cy][cz + 1],dcg -> dense3[cx - 1][cy + 1][cz + 1]});
            break;
            
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex3(birthday, index, 3);
            return true;
          }
        }
        return false;
      }
    return false;
    }
  }
  
  BirthdayIndex3 getNextCoface() { return nextCoface; }
};

/*****joint_pairs*****/
class JointPairs3
{
  int n; // the number of cubes
  int ctr_moi;
  int ax, ay, az;
  DenseCubicalGrids3* dcg;
  ColumnsToReduce3* ctr;
  vector<WritePairs3> *wp;
  Vertices* vtx;
  double u, v;
  vector<int64_t> cubes_edges;
  vector<BirthdayIndex3> dim1_simplex_list;
  
public:
  JointPairs3(DenseCubicalGrids3* _dcg, ColumnsToReduce3* _ctr, vector<WritePairs3> &_wp)
  {
    dcg = _dcg;
    ax = dcg -> ax;
    ay = dcg -> ay;
    az = dcg -> az;
    ctr = _ctr; // ctr is "0-dim"simplex list.
    ctr_moi = ctr -> max_of_index;
    n = ctr -> columns_to_reduce.size();
    
    wp = &_wp;
    vtx = new Vertices();
    
    for(int x = 1; x <= ax; ++x)
      for(int y = 1; y <= ay; ++y)
        for(int z = 1; z <= az; ++z)
          for(int type = 0; type < 3; ++type)
          {
            int index = x | (y << 9) | (z << 18) | (type << 27);
            double birthday = dcg -> getBirthday(index, 1);
            
            if(birthday < dcg -> threshold)
              dim1_simplex_list.push_back(BirthdayIndex3(birthday, index, 1));
          }
    
    sort(dim1_simplex_list.rbegin(), dim1_simplex_list.rend(), BirthdayIndex3Comparator());
  }
  
  // free pointers
  ~JointPairs3()
  {
    delete vtx;
  }
  
  void joint_pairs_main()
  {
    cubes_edges.resize(2);
    UnionFind3 dset(ctr_moi, dcg);
    ctr -> columns_to_reduce.clear();
    ctr -> dim = 1;
    double min_birth = dcg -> threshold;
    
    for(auto e : dim1_simplex_list)
    {
      cubes_edges.clear();
      dcg -> GetSimplexVertices(e.getIndex(), 1, vtx);
      
      cubes_edges.push_back(vtx -> vertex[0] -> getIndex());
      cubes_edges.push_back(vtx -> vertex[1] -> getIndex());
      
      u = dset.find(cubes_edges[0]);
      v = dset.find(cubes_edges[1]);
      
      if(min_birth >= min(dset.birthtime[u], dset.birthtime[v]))
        min_birth = min(dset.birthtime[u], dset.birthtime[v]);
      
      if(u != v)
      {
        double birth = max(dset.birthtime[u], dset.birthtime[v]);
        double death = max(dset.time_max[u], dset.time_max[v]);
        
        if (birth == death)
          dset.link(u, v);
        else
        {
          wp -> push_back(WritePairs3(0, birth, death));
          dset.link(u, v);
        }
      }
      else // If two values have same "parent", these are potential edges which make a 2-simplex.
        ctr -> columns_to_reduce.push_back(e);
    }
    
    wp -> push_back(WritePairs3(-1, min_birth, dcg -> threshold));
    sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndex3Comparator());
  }
};

/*****compute_pairs*****/
template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};

class ComputePairs3
{
public:
  DenseCubicalGrids3* dcg;
  ColumnsToReduce3* ctr;
  hash_map<int, int> pivot_column_index;
  int ax, ay, az;
  int dim;
  vector<WritePairs3> *wp;
  
  ComputePairs3(DenseCubicalGrids3* _dcg, ColumnsToReduce3* _ctr, vector<WritePairs3> &_wp)
  {
    dcg = _dcg;
    ctr = _ctr;
    dim = _ctr -> dim;
    wp = &_wp;
    
    ax = _dcg -> ax;
    ay = _dcg -> ay;
    az = _dcg -> az;
  }
  
  void compute_pairs_main()
  {
    pivot_column_index = hash_map<int, int>();
    vector<BirthdayIndex3> coface_entries;
    auto ctl_size = ctr -> columns_to_reduce.size();
    SimplexCoboundaryEnumerator3 cofaces;
    unordered_map<int, priority_queue<BirthdayIndex3, vector<BirthdayIndex3>, BirthdayIndex3Comparator>> recorded_wc;
    
    pivot_column_index.reserve(ctl_size);
    recorded_wc.reserve(ctl_size);
    
    for(int i = 0; i < ctl_size; ++i) {
      if (i % 5000 == 0) {
        Rcpp::checkUserInterrupt();
      }
      
      auto column_to_reduce = ctr -> columns_to_reduce[i]; 
      priority_queue<BirthdayIndex3, vector<BirthdayIndex3>, BirthdayIndex3Comparator> 
        working_coboundary;
      double birth = column_to_reduce.getBirthday();
      
      int j = i;
      BirthdayIndex3 pivot(0, -1, 0);
      bool might_be_apparent_pair = true;
      bool goto_found_persistence_pair = false;
      
      do {
        auto simplex = ctr -> columns_to_reduce[j]; // get CTR[i] 
        coface_entries.clear();
        cofaces.setSimplexCoboundaryEnumerator3(simplex, dcg);// make coface data
        
        while (cofaces.hasNextCoface() && !goto_found_persistence_pair) { // repeat there remains a coface
          BirthdayIndex3 coface = cofaces.getNextCoface();
          coface_entries.push_back(coface);
          if (might_be_apparent_pair && (simplex.getBirthday() == coface.getBirthday())) { // If bt is the same, go thru
            if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) { // If coface is not in pivot list
              pivot.copyBirthdayIndex3(coface); // I have a new pivot
              goto_found_persistence_pair = true; // goto (B)
            } else { // If pivot list contains this coface,
              might_be_apparent_pair = false; // goto (A)
            }
          }
        }
        
        if (!goto_found_persistence_pair) { // (A) If pivot list contains this coface,
          auto findWc = recorded_wc.find(j); // we seek wc list by 'j'
          
          if(findWc != recorded_wc.end()){ // If the pivot is old,
            auto wc = findWc -> second;
            while(!wc.empty()){ // we push the data of the old pivot's wc
              auto e = wc.top();
              working_coboundary.push(e);
              wc.pop();
            }
          } else { // If the pivot is new,
            for(auto e : coface_entries){ // making wc here
              working_coboundary.push(e);
            }
          }
          pivot = get_pivot(working_coboundary); // getting a pivot from wc
          
          if (pivot.getIndex() != -1) { // When I have a pivot, ...
            auto pair = pivot_column_index.find(pivot.getIndex());
            if (pair != pivot_column_index.end()) {	// If the pivot already exists, go on the loop 
              j = pair -> second;
              continue;
            } else { // If the pivot is new, 
              // I record this wc into recorded_wc, and 
              recorded_wc.insert(make_pair(i, working_coboundary));
              // I output PP as WritePairs
              double death = pivot.getBirthday();
              outputPP(dim, birth, death);
              pivot_column_index.insert(make_pair(pivot.getIndex(), i));
              break;
            }
          } else { // If wc is empty, I output a PP as [birth,) 
            outputPP(-1, birth, dcg -> threshold);
            break;
          }
        } else { // (B) I have a new pivot and output PP as Writepairs 
          double death = pivot.getBirthday();
          outputPP(dim, birth, death);
          pivot_column_index.insert(make_pair(pivot.getIndex(), i));
          break;
        }			
        
      } while (true);
    }
  }
  
  void outputPP(int _dim, double _birth, double _death)
  {
    if(_birth != _death)
    {
      if(_death != dcg -> threshold)
        wp -> push_back(WritePairs3(_dim, _birth, _death));
      else
        wp -> push_back(WritePairs3(-1, _birth, dcg -> threshold));
    }
  }
  
  BirthdayIndex3 pop_pivot(priority_queue<BirthdayIndex3, vector<BirthdayIndex3>, BirthdayIndex3Comparator>& column)
  {
    if (column.empty())
      return BirthdayIndex3(0, -1, 0);
    else
    {
      auto pivot = column.top();
      column.pop();
      
      while (!column.empty() && column.top().index == pivot.getIndex())
      {
        column.pop();
        if (column.empty())
          return BirthdayIndex3(0, -1, 0);
        else
        {
          pivot = column.top();
          column.pop();
        }
      }
      return pivot;
    }
  }
  
  BirthdayIndex3 get_pivot(priority_queue<BirthdayIndex3, vector<BirthdayIndex3>, BirthdayIndex3Comparator>& column)
  {
    BirthdayIndex3 result = pop_pivot(column);
    
    if (result.getIndex() != -1)
      column.push(result);
    
    return result;
  }
  
  void assemble_columns_to_reduce()
  {
    ++dim;
    ctr -> dim = dim;
    
    if (dim == 1)
    { 
      ctr -> columns_to_reduce.clear();
      for(int z = 1; z <= az; ++z)
        for (int y = 1; y <= ay; ++y)
          for (int x = 1; x <= ax; ++x)
            for (int m = 0; m < 3; ++m) // the number of type
            {
              double index = x | (y << 9) | (z << 18) | (m << 27);
              if (pivot_column_index.find(index) == pivot_column_index.end())
              {
                double birthday = dcg -> getBirthday(index, 1);
                if (birthday != dcg -> threshold)
                  ctr -> columns_to_reduce.push_back(BirthdayIndex3(birthday, index, 1));
              }
            }
    }
    else if (dim == 2)
    { 
      ctr -> columns_to_reduce.clear();
      for(int z = 1; z <= az; ++z)
        for (int y = 1; y <= ay; ++y)
          for (int x = 1; x <= ax; ++x)
            for (int m = 0; m < 3; ++m) // the number of type
            {
              double index = x | (y << 9) | (z << 18) | (m << 27);
              if (pivot_column_index.find(index) == pivot_column_index.end())
              {
                double birthday = dcg -> getBirthday(index, 2);
                if (birthday != dcg -> threshold)
                  ctr -> columns_to_reduce.push_back(BirthdayIndex3(birthday, index, 2));
              }
            }
    }
    sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndex3Comparator());
  }
};

// method == 0 --> LINKFIND
// method == 1 --> COMPUTEPAIRS
// [[Rcpp::export]]
Rcpp::NumericMatrix cubical_3dim(Rcpp::NumericVector& image, double threshold, int method, int nx, int ny, int nz)
{
  vector<WritePairs3> writepairs; // dim birth death
  writepairs.clear();
  
  DenseCubicalGrids3* dcg = new DenseCubicalGrids3(image, threshold, nx, ny, nz);
  ColumnsToReduce3* ctr = new ColumnsToReduce3(dcg);
  
  switch (method)
  {
    case 0:
    {
      JointPairs3* jp = new JointPairs3(dcg, ctr, writepairs);
      jp -> joint_pairs_main(); // dim0
      
      ComputePairs3* cp = new ComputePairs3(dcg, ctr, writepairs);
      cp -> compute_pairs_main(); // dim1
      cp -> assemble_columns_to_reduce();
      
      cp -> compute_pairs_main(); // dim2
      
      // free pointers
      delete jp;
      delete cp;
      
      break;
    }
      
    case 1:
    {
      ComputePairs3* cp = new ComputePairs3(dcg, ctr, writepairs);
      cp -> compute_pairs_main(); // dim0
      cp -> assemble_columns_to_reduce();
      
      cp -> compute_pairs_main(); // dim1
      cp -> assemble_columns_to_reduce();
      
      cp -> compute_pairs_main(); // dim2
      
      // free pointers
      delete cp;
      
      break;
    }
  }
  
  // free pointers
  delete dcg;
  delete ctr;
  
  Rcpp::NumericMatrix ans(writepairs.size(), 3);
  for (int i = 0; i < ans.nrow(); i++)
  {
    ans(i, 0) = writepairs[i].getDimension();
    ans(i, 1) = writepairs[i].getBirth();
    ans(i, 2) = writepairs[i].getDeath();
  }
  
  return ans;
}
