/*
 Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
 
 This file is part of CubicalRipser_4dim.
 
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
 4 dimensional voxel data. The input data format conforms to that of DIPHA or of PERSEUS.
 See more descriptions in README.
 
 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
 You should have received a copy of the GNU Lesser General Public License along
 with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdint>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <queue>
#include <Rcpp.h>

using namespace std;

/*****birthday_index*****/
class BirthdayIndex4
{
  
public:
  double birthday;
  int index;
  int dim;
  
  // constructors
  BirthdayIndex4() : birthday(0), index(-1), dim(1) {};
  BirthdayIndex4(const BirthdayIndex4& b) : birthday(b.birthday), index(b.index), dim(b.dim) {};
  BirthdayIndex4(double _b, int _i, int _d) : birthday(_b), index(_i), dim(_d) {};
  
  // copy into this object
  void copyBirthdayIndex4(BirthdayIndex4 v)
  {
    birthday = v.birthday;
    index = v.index;
    dim = v.dim;
  }
  
  // getters
  double getBirthday() { return birthday; }
  long getIndex() { return index; }
  int getDimension() { return dim; }
};

// utility function to shorten code
bool compareBday(const BirthdayIndex4& o1, const BirthdayIndex4& o2)
{
  return (o1.birthday == o2.birthday ? o1.index < o2.index : o1.birthday > o2.birthday);
}

struct BirthdayIndex4Comparator
{
  bool operator()(const BirthdayIndex4& o1, const BirthdayIndex4& o2) const
  {
    return compareBday(o1, o2);
  }
};

struct BirthdayIndex4InverseComparator
{
  bool operator()(const BirthdayIndex4& o1, const BirthdayIndex4& o2) const
  {
    return !compareBday(o1, o2);
  }
};

/*****write_pairs*****/
class WritePairs4
{
public:
  int64_t dim;
  double birth;
  double death;
  
  WritePairs4(int64_t _dim, double _birth, double _death) : dim(_dim), birth(_birth), death(_death) {};
  
  int64_t getDimension() { return dim; }
  double getBirth() { return birth; }
  double getDeath() { return death; }
};

/*****dense_cubical_grids*****/
const int MAX_SIZE = 64;
const int EXPONENT = 6;

class DenseCubicalGrids4
{
public:
  double threshold;
  int dim;
  int ax, ay, az, aw;
  double dense4[MAX_SIZE][MAX_SIZE][MAX_SIZE][MAX_SIZE];
  
  DenseCubicalGrids4(const Rcpp::NumericVector& image, double _threshold, int nx, int ny, int nz, int nt) : threshold(_threshold), ax(nx), ay(ny), az(nz), aw(nt)
  {
    dim = 4;
    
    // set everything to threshold
    for (int i = 0; i < MAX_SIZE; i++)
      for (int j = 0; j < MAX_SIZE; j++)
        for (int k = 0; k < MAX_SIZE; k++)
          for (int l = 0; l < MAX_SIZE; l++)
            dense4[i][j][k][l] = threshold;
    
    // set values based on image
    for (int i = 0; i < ax * ay * az * aw; i++)
      dense4[i % ax + 1][i / ax % ay + 1][i / (ax * ay) % az + 1][i / (ax * ay * az) % aw + 1] = image(i);
  }
  double getBirthday(int index, int dim)
  {
    int cx = index & (MAX_SIZE - 1),
        cy = (index >> EXPONENT) & (MAX_SIZE - 1),
        cz = (index >> (2 * EXPONENT)) & (MAX_SIZE - 1),
        cw = (index >> (3 * EXPONENT)) & (MAX_SIZE - 1),
        cm = (index >> (4 * EXPONENT)) & 0x0f;
    
    switch (dim)
    {
      case 0:
        return dense4[cx][cy][cz][cw];
      case 1:
        switch (cm)
        {
          case 0:
            return max(dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw]);
          case 1:
            return max(dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw]);
          case 2:
            return max(dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw]);
          default:
            return max(dense4[cx][cy][cz][cw], dense4[cx][cy][cz][cw + 1]);
        }
      case 2:
        switch (cm)
        {
          case 0: // x - y (fix z, w)
            return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw]});
          case 1: // z - x (fix y, w)
            return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], dense4[cx + 1][cy][cz][cw]});
          case 2: // y - z (fix x, w)
            return max({dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw], dense4[cx][cy + 1][cz + 1][cw], dense4[cx][cy][cz + 1][cw]});
          case 3: // x - w
            return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], dense4[cx + 1][cy][cz][cw + 1], dense4[cx][cy][cz][cw + 1]});
          case 4: // y - w
            return max({dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw + 1], dense4[cx][cy][cz][cw + 1]});
          case 5: // z - w
            return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw], dense4[cx][cy][cz + 1][cw + 1], dense4[cx][cy][cz][cw + 1]});
        }
      case 3:
        switch (cm)
        {
          case 0: // x - y - z
            return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
                        dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], dense4[cx + 1][cy + 1][cz + 1][cw], dense4[cx][cy + 1][cz + 1][cw]});
          case 1: // x - y - w
            return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
                        dense4[cx][cy][cz][cw + 1], dense4[cx + 1][cy][cz][cw + 1], dense4[cx + 1][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1]});
          case 2: // x - z - w
            return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], dense4[cx + 1][cy][cz][cw + 1], dense4[cx][cy][cz][cw + 1],
                        dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw + 1], dense4[cx][cy][cz + 1][cw + 1]});
          case 3: // y - z - w
            return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw],
                        dense4[cx][cy][cz + 1][cw], dense4[cx][cy][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw]});
        }
      case 4:
        return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
                    dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], dense4[cx + 1][cy + 1][cz + 1][cw], dense4[cx][cy + 1][cz + 1][cw],
                    dense4[cx][cy][cz][cw + 1], dense4[cx + 1][cy][cz][cw + 1], dense4[cx + 1][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1],
                    dense4[cx][cy][cz + 1][cw + 1], dense4[cx + 1][cy][cz + 1][cw + 1], dense4[cx + 1][cy + 1][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw + 1]});
    }
    return threshold;
  }
};

/*****columns_to_reduce*****/
class ColumnsToReduce4
{
public:
  vector<BirthdayIndex4> columns_to_reduce;
  int dim;
  int max_of_index;
  
  ColumnsToReduce4(DenseCubicalGrids4* _dcg)
  {
    dim = 0;
    int ax = _dcg->ax;
    int ay = _dcg->ay;
    int az = _dcg->az;
    int aw = _dcg->aw;
    max_of_index = MAX_SIZE * MAX_SIZE * MAX_SIZE * (aw + 2);
    int index;
    double birthday;
    
    for(int w = aw; w > 0; --w)
      for(int z = az; z > 0; --z)
        for (int y = ay; y > 0; --y)
          for (int x = ax; x > 0; --x)
          {
            birthday = _dcg -> dense4[x][y][z][w];
            index = x | (y << EXPONENT) | (z << (2 * EXPONENT)) | (w << (3 * EXPONENT));
            if (birthday != _dcg -> threshold) columns_to_reduce.push_back(BirthdayIndex4(birthday, index, 0));
          }
    
    sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndex4Comparator());
  }
  int size() { return columns_to_reduce.size(); }
};

/*****simplex_coboundary_estimator*****/
class SimplexCoboundaryEnumerator4
{
public:
  BirthdayIndex4 simplex;
  DenseCubicalGrids4* dcg;
  int dim;
  double birthtime;
  int ax, ay, az, aw;
  int cx, cy, cz, cw, cm;
  int count;
  BirthdayIndex4 nextCoface;
  double threshold; 
  
  SimplexCoboundaryEnumerator4() { nextCoface = BirthdayIndex4(0, -1, 1); }
  void setSimplexCoboundaryEnumerator4(BirthdayIndex4 _s, DenseCubicalGrids4* _dcg)
  {
    simplex = _s;
    dcg = _dcg;
    dim = simplex.dim;
    birthtime = simplex.birthday;
    ax = _dcg->ax;
    ay = _dcg->ay;
    az = _dcg->az;
    aw = _dcg->aw;
    
    cx = simplex.index & (MAX_SIZE - 1);
    cy = (simplex.index >> EXPONENT) & (MAX_SIZE - 1);
    cz = (simplex.index >> (2 * EXPONENT)) & (MAX_SIZE - 1);
    cw = (simplex.index >> (3 * EXPONENT)) & (MAX_SIZE - 1);
    cm = (simplex.index >> (4 * EXPONENT)) & 0x0f;
    
    threshold = _dcg->threshold;
    count = 0;
  }
  bool hasNextCoface()
  {
    int index = 0;
    double birthday = 0;
    switch (dim)
    {
      case 0: // dim0
        for (int i = count; i < 8; ++i)
        {
          switch (i)
          {
            case 0: // w +
              index = (3 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx][cy][cz][cw + 1]);
              break;
          
            case 1: // w -
              index = (3 << (4 * EXPONENT)) | ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx][cy][cz][cw - 1]);
              break;
          
            case 2: // z +
              index = (2 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx][cy][cz + 1][cw]);
              break;
          
            case 3: // z -
              index = (2 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | ((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx][cy][cz - 1][cw]);
              break;
          
            case 4: // y +
              index = (1 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx][cy + 1][cz][cw]);
              break;
          
            case 5: // y -
              index = (1 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx][cy - 1][cz][cw]);
              break;
          
            case 6: // x +
              index = (0 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
              birthday = max(birthtime, dcg -> dense4[cx + 1][cy][cz][cw]);
              break;
          
            case 7: // x -
              index = (0 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
              birthday = max(birthtime, dcg -> dense4[cx - 1][cy][cz][cw]);
              break;
          }
          if (birthday != threshold)
          {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 1);
            return true;
          }
        }
        return false;
      case 1:// dim1
        switch (cm)
        {
          case 0: // dim1 type0 (x-axis -> )
            for(int i = count; i < 6; ++i)
            {
              switch (i)
              {
                case 0: // x - w +
                  index = (3 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1]});
                  break;

                case 1: // x - w -
                  index = (3 << (4 * EXPONENT)) | ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1]});
                  break;
            
                case 2: // x - z +
                  index = (1 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw]});
                  break;
            
                case 3: // x - z -
                  index = (1 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw]});
                  break;
            
                case 4: // x - y +
                  index = (0 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw]});
                  break;
            
                case 5: // x - y -
                  index = (0 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw]});
                  break;
              }
          
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex4(birthday, index, 2);
                return true;
              }
            }
            return false;
        
          case 1: // dim1 type1 (y-axis -> )
            for (int i = count; i < 6; ++i)
            {
              switch (i)
              {
                case 0: // y - w +
                  index = (4 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy + 1][cz][cw + 1]});
                  break;
            
                case 1: // y - w -
                  index = (4 << (4 * EXPONENT)) | ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy + 1][cz][cw - 1]});
                  break;
            
                case 2: // y - z +
                  index = (2 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw]});
                  break;
            
                case 3: // y - z -
                  index = (2 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy + 1][cz - 1][cw]});
                  break;
            
                case 4: // y - x +
                  index = (0 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw]});
                  break;
            
                case 5: // y - x -
                  index = (0 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
                  birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw]});
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex4(birthday, index, 2);
                return true;
              }
            }
            return false;
        
          case 2: // dim1 type2 (z-axis -> )
            for (int i = count; i < 6; ++i)
            {
              switch (i)
              {
                case 0: // z - w +
                  index = (5 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy][cz + 1][cw + 1]});
                  break;
            
                case 1: // z - w -
                  index = (5 << (4 * EXPONENT)) | ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy][cz + 1][cw - 1]});
                  break;
            
                case 2: // z - y +
                  index = (2 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw]});
                  break;
            
                case 3: // z - y -
                  index = (2 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz + 1][cw]});
                  break;
            
                case 4: // z - x +
                  index = (1 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw]});
                  break;
            
                case 5: // z - x -
                  index = (1 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
                  birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz + 1][cw]});
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex4(birthday, index, 2);
                return true;
              }
            }
            return false;
        
          case 3: // dim1 type3 (w-axis -> )
            for (int i = count; i < 6; ++i)
            {
              switch (i)
              {
                case 0: // w - z +
                  index = (5 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy][cz + 1][cw + 1]});
                  break;
            
                case 1: // w - z -
                  index = (5 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) | ((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy][cz - 1][cw + 1]});
                  break;
            
                case 2: // w - y +
                  index = (4 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz][cw + 1]});
                  break;
            
                case 3: // w - y -
                  index = (4 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz][cw + 1]});
                  break;
            
                case 4: // w - x +
                  index = (3 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
                  birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz][cw + 1]});
                  break;
            
                case 5: // w - x -
                  index = (3 << (4 * EXPONENT)) | (cw << (3 * EXPONENT)) |(cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
                  birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz][cw + 1]});
                  break;
              }
              if (birthday != threshold)
              {
                count = i + 1;
                nextCoface = BirthdayIndex4(birthday, index, 2);
                return true;
              }
            }
            return false;
      }
      return false;
    case 2:// dim2
      switch (cm) {
      case 0: // dim2 type0 (fix x - y)
        for(int i = count; i < 4; ++i){
          switch(i){
          case 0: // w +
            index = (1 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
                           dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
            break;
            
          case 1: // w -
            index = (1 << (4 * EXPONENT))| ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
                           dcg -> dense4[cx][cy + 1][cz][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz][cw - 1]});
            break;
            
          case 2: // z +
            index = (0 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                           dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
            break;
            
          case 3: // z -
            index = (0 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | ((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
                           dcg -> dense4[cx][cy + 1][cz - 1][cw],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw]});
            break;
            
          }
          
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 1: // dim2 type1 (fix x - z)
        for(int i = count; i < 4; ++i){
          switch(i){
          case 0: // w +
            index = (2 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
                           dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
            break;
            
          case 1: // w -
            index = (2 << (4 * EXPONENT))| ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
                           dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy][cz + 1][cw - 1]});
            break;
            
          case 2: // y +
            index = (0 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                           dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
            break;
            
          case 3: // y -
            index = (0 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
                           dcg -> dense4[cx][cy - 1][cz + 1][cw],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw]});
            break;
            
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 2: // dim2 type3 (fix y - z)
        for(int i = count; i < 4; ++i){
          switch(i){
          case 0: // w +
            index = (3 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy + 1][cz][cw + 1], 
                           dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // w -
            index = (3 << (4 * EXPONENT))| ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy + 1][cz][cw - 1], 
                           dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx][cy + 1][cz + 1][cw - 1]});
            break;
            
          case 2: // x +
            index = (0 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                           dcg -> dense4[cx + 1][cy][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
            break;
            
          case 3: // x -
            index = (0 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
            birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
                           dcg -> dense4[cx - 1][cy][cz + 1][cw],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw]});
            break;
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 3: // dim2 type2 (fix x - w)
        for(int i = count; i < 4; ++i){
          switch(i){
          case 0: // z +
            index = (2 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                           dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
            break;
            
          case 1: // z -
            index = (2 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | ((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
                           dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy][cz - 1][cw + 1]});
            break;
            
          case 2: // y +
            index = (1 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                           dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
            break;
            
          case 3: // y -
            index = (1 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
                           dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz][cw + 1]});
            break;
            
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 4: // dim2 type4 (fix y - w)
        for(int i = count; i < 4; ++i){
          switch(i){
          case 0: // z +
            index = (3 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw], 
                           dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // z -
            index = (3 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | ((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy + 1][cz - 1][cw], 
                           dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx][cy + 1][cz - 1][cw + 1]});
            break;
            
          case 2: // x +
            index = (1 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                           dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
            break;
            
          case 3: // x -
            index = (1 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
            birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
                           dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz][cw + 1]});
            break;
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 3);
            return true;
          }
        }
        return false;
        
      case 5: // dim2 type5 (fix z - w)
        for(int i = count; i < 4; ++i){
          switch(i){
          case 0: // y +
            index = (3 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw], 
                           dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // y -
            index = (3 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz + 1][cw], 
                           dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx][cy - 1][cz + 1][cw + 1]});
            break;
            
          case 2: // x +
            index = (2 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                           dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
            break;
            
          case 3: // x -
            index = (2 << (4 * EXPONENT))| (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
            birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz + 1][cw], 
                           dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy][cz + 1][cw + 1]});
            break;
          }
          if (birthday != threshold) {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 3);
            return true;
          }
        }
        return false;
      }
      return false;
    case 3:// dim3
      switch (cm) {
      case 0: // dim3 type0 (x - y - z)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // w +
            index = (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
                            dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
                            dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],
                            dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // w -
            index = ((cw - 1) << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
                            dcg -> dense4[cx][cy + 1][cz][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz][cw - 1],
                            dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy][cz + 1][cw - 1],
                            dcg -> dense4[cx][cy + 1][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw - 1]});
            break;
          }
          if (birthday != threshold)
          {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 4);
            return true;
          }
        }
        return false;
        
      case 1: // dim3 type1 (x - y - w)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // z +
            index = (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                            dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
                            dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],
                            dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // z -
            index = (cw << (3 * EXPONENT)) | ((cz - 1) << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
                            dcg -> dense4[cx][cy + 1][cz - 1][cw],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw],
                            dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy][cz - 1][cw + 1],
                            dcg -> dense4[cx][cy + 1][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw + 1]});
            break;
          }
          if (birthday != threshold)
          {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 4);
            return true;
          }
        }
        return false;
        
      case 2: // dim3 type2 (x - z - w)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // y +
            index = (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                            dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
                            dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
                            dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // y -
            index = (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | ((cy - 1) << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
                            dcg -> dense4[cx][cy - 1][cz + 1][cw],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw],
                            dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz][cw + 1],
                            dcg -> dense4[cx][cy - 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw + 1]});
            break;
          }
          if (birthday != threshold)
          {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 4);
            return true;
          }
        }
        return false;
        
      case 3: // dim3 type3 (y - z - w)
        for(int i = count; i < 2; ++i){
          switch(i){
          case 0: // x +
            index = (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | cx;
            birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                            dcg -> dense4[cx + 1][cy][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
                            dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
                            dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
            break;
            
          case 1: // x -
            index = (cw << (3 * EXPONENT)) | (cz << (2 * EXPONENT)) | (cy << EXPONENT) | (cx - 1);
            birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
                            dcg -> dense4[cx - 1][cy][cz + 1][cw],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw],
                            dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz][cw + 1],
                            dcg -> dense4[cx - 1][cy][cz + 1][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw + 1]});
            break;
          }
          if (birthday != threshold)
          {
            count = i + 1;
            nextCoface = BirthdayIndex4(birthday, index, 4);
            return true;
          }
        }
        return false;
      }
      return false;
    }
    return false;
  }
  BirthdayIndex4 getNextCoface() { return nextCoface; }
};

/*****union_find*****/
class UnionFind4
{
public:
  int max_of_index;
  vector<int> parent;
  vector<double> birthtime;
  vector<double> time_max;
  DenseCubicalGrids4* dcg;
  
  UnionFind4(int moi, DenseCubicalGrids4* _dcg) : parent(moi), birthtime(moi), time_max(moi) // "n" is the number of cubes.
  {
    dcg = _dcg;
    max_of_index = moi;
    
    for(int i = 0; i < moi; ++i){
      parent[i] = i;
      birthtime[i] = dcg->getBirthday(i, 0);
      time_max[i] = dcg->getBirthday(i, 0);
    }
  }
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
    if (birthtime[x] > birthtime[y]) {
      parent[x] = y; 
      birthtime[y] = min(birthtime[x], birthtime[y]);
      time_max[y] = max(time_max[x], time_max[y]);
    } else if(birthtime[x] < birthtime[y]) {
      parent[y] = x;
      birthtime[x] = min(birthtime[x], birthtime[y]);
      time_max[x] = max(time_max[x], time_max[y]);
    } else { //birthtime[x] == birthtime[y]
      parent[x] = y;
      time_max[y] = max(time_max[x], time_max[y]);
    }
  }
};

class JointPairs4
{
  int n; // the number of cubes
  int ctr_moi;
  int ax, ay, az, aw;
  DenseCubicalGrids4* dcg;
  ColumnsToReduce4* ctr;
  vector<WritePairs4> *wp;
  double u, v;
  vector<int64_t> cubes_edges;
  vector<BirthdayIndex4> dim1_simplex_list;
  
public:
  JointPairs4(DenseCubicalGrids4* _dcg, ColumnsToReduce4* _ctr, vector<WritePairs4> &_wp)
  {
    dcg = _dcg;
    ax = dcg -> ax;
    ay = dcg -> ay;
    az = dcg -> az;
    aw = dcg -> aw;
    ctr = _ctr; // ctr is "dim0"simplex list.
    ctr_moi = ctr -> max_of_index;
    n = ctr -> columns_to_reduce.size();
    
    wp = &_wp;
    
    for(int x = 1; x <= ax; ++x)
      for(int y = 1; y <= ay; ++y)
        for(int z = 1; z <= az; ++z)
          for(int w = 1; w <= aw; ++w)
            for(int type = 0; type < 4; ++type) // change
            {
              int index = x | (y << EXPONENT) | (z << (2 * EXPONENT)) | (w << (3 * EXPONENT)) | (type << (4 * EXPONENT));
              double birthday = dcg -> getBirthday(index, 1);
              if(birthday < dcg -> threshold) dim1_simplex_list.push_back(BirthdayIndex4(birthday, index, 1));
            }
            
    sort(dim1_simplex_list.rbegin(), dim1_simplex_list.rend(), BirthdayIndex4Comparator());
  }
  void joint_pairs_main()
  {
    UnionFind4 dset(ctr_moi, dcg);
    ctr -> columns_to_reduce.clear();
    ctr -> dim = 1;
    double min_birth = dcg -> threshold;
    
    for (BirthdayIndex4 e : dim1_simplex_list)
    {
      int index = e.getIndex();
      int cx = index & (MAX_SIZE - 1);
      int cy = (index >> EXPONENT) & (MAX_SIZE - 1);
      int cz = (index >> (2 * EXPONENT)) & (MAX_SIZE - 1);
      int cw = (index >> (3 * EXPONENT)) & (MAX_SIZE - 1);
      int cm = (index >> (4 * EXPONENT)) & 0x0f;
      
      int ce0=0, ce1 =0;
      switch (cm)
      {
      case 0:
        ce0 =  cx | (cy << EXPONENT) | (cz << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        ce1 =  (cx + 1) | (cy << EXPONENT) | (cz << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        break;
      case 1:
        ce0 =  cx | (cy << EXPONENT) | (cz << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        ce1 =  cx | ((cy + 1) << EXPONENT) | (cz << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        break;
      case 2:
        ce0 =  cx | (cy << EXPONENT) | (cz << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        ce1 =  cx | (cy << EXPONENT) | ((cz + 1) << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        break;
      case 3:
        ce0 =  cx | (cy << EXPONENT) | (cz << (2 * EXPONENT)) | (cw << (3 * EXPONENT));
        ce1 =  cx | (cy << EXPONENT) | (cz << (2 * EXPONENT)) | ((cw + 1) << (3 * EXPONENT));
        break;
      }
      u = dset.find(ce0);
      v = dset.find(ce1);
      if(min_birth >= min(dset.birthtime[u], dset.birthtime[v]))
      {
        min_birth = min(dset.birthtime[u], dset.birthtime[v]);
      }
      
      if(u != v)
      {
        double birth = max(dset.birthtime[u], dset.birthtime[v]);
        double death = max(dset.time_max[u], dset.time_max[v]);
        
        if(birth == death){
          dset.link(u, v);
        } else {
          wp -> push_back(WritePairs4(0, birth, death));
          dset.link(u, v);
        }
      } else { // If two values have same "parent", these are potential edges which make a 2-simplex.
        ctr -> columns_to_reduce.push_back(e);
      }
    }
    
    wp -> push_back(WritePairs4(-1, min_birth, dcg->threshold));
    sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndex4Comparator());
  }
};

/*****compute_pairs*****/
template <class Key, class T> class hash_map4 : public std::unordered_map<Key, T> {};

class ComputePairs4
{
public:
  DenseCubicalGrids4* dcg;
  ColumnsToReduce4* ctr;
  hash_map4<int, int> pivot_column_index;
  int ax, ay, az, aw;
  int dim;
  vector<WritePairs4> *wp;
  
  ComputePairs4(DenseCubicalGrids4* _dcg, ColumnsToReduce4* _ctr, vector<WritePairs4> &_wp)
  {
    dcg = _dcg;
    ctr = _ctr;
    dim = _ctr -> dim;
    wp = &_wp;
    
    ax = _dcg -> ax;
    ay = _dcg -> ay;
    az = _dcg -> az;
    aw = _dcg -> aw;
  }
  void compute_pairs_main()
  {
    vector<BirthdayIndex4> coface_entries;
    SimplexCoboundaryEnumerator4 cofaces;
    unordered_map<int, priority_queue<BirthdayIndex4, vector<BirthdayIndex4>, BirthdayIndex4Comparator>> recorded_wc;
    
    pivot_column_index = hash_map4<int, int>();
    auto ctl_size = ctr -> columns_to_reduce.size();
    pivot_column_index.reserve(ctl_size);
    recorded_wc.reserve(ctl_size);
    
    
    for(int i = 0; i < ctl_size; ++i) {
      if (i % 10000 == 0) {
        Rcpp::checkUserInterrupt();
      }
      
      auto column_to_reduce = ctr -> columns_to_reduce[i]; 
      priority_queue<BirthdayIndex4, vector<BirthdayIndex4>, BirthdayIndex4Comparator> working_coboundary;
      double birth = column_to_reduce.getBirthday();
      
      int j = i;
      BirthdayIndex4 pivot(0, -1, 0);
      bool might_be_apparent_pair = true;
      bool goto_found_persistence_pair = false;
      
      do {
        auto simplex = ctr -> columns_to_reduce[j];
        coface_entries.clear();
        cofaces.setSimplexCoboundaryEnumerator4(simplex, dcg);
        
        while (cofaces.hasNextCoface() && !goto_found_persistence_pair) {
          BirthdayIndex4 coface = cofaces.getNextCoface();
          coface_entries.push_back(coface);
          if (might_be_apparent_pair && (simplex.getBirthday() == coface.getBirthday())) {
            if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) {
              pivot.copyBirthdayIndex4(coface);
              goto_found_persistence_pair = true;// goto (B)
            } else {
              might_be_apparent_pair = false;// goto(A)
            }
          }
        }
        
        if (!goto_found_persistence_pair) {// (A) I haven't had a pivot
          auto findWc = recorded_wc.find(j);
          if(findWc != recorded_wc.end()){// if the pivot is old,
            auto wc = findWc->second;
            while (!wc.empty()){// we push the data of the old pivot's wc
              auto e = wc.top();
              working_coboundary.push(e);
              wc.pop();
            }
          } else {
            for (auto e : coface_entries) {// making wc here
              working_coboundary.push(e);
            }
          }
          pivot = get_pivot(working_coboundary);// getting a pivot from wc
          
          if (pivot.getIndex() != -1) {//When I have a pivot, ...
            auto pair = pivot_column_index.find(pivot.getIndex());
            if (pair != pivot_column_index.end()) {	// if the pivot already exists, go on the loop 
              j = pair->second;
              continue;
            } else {// if the pivot is new, 
              // I record this wc into recorded_wc, and 
              recorded_wc.insert(make_pair(i, working_coboundary));
              // I output PP as WritePairs
              double death = pivot.getBirthday();
              outputPP(dim,birth,death);
              pivot_column_index.insert(make_pair(pivot.getIndex(), i));
              break;
            }
          } else {// if wc is empty, I output a PP as [birth,) 
            outputPP(-1, birth, dcg->threshold);
            break;
          }
        } else {// (B) I have a pivot and output PP as WritePairs 
          double death = pivot.getBirthday();
          outputPP(dim,birth,death);
          pivot_column_index.insert(make_pair(pivot.getIndex(), i));
          break;
        }			
        
      } while (true);
    }
  }
  void outputPP(int _dim, double _birth, double _death)
  {
    if(_birth != _death){
      if(_death != dcg -> threshold) {
        wp -> push_back(WritePairs4(_dim, _birth, _death));
      } else {
        wp -> push_back(WritePairs4(-1, _birth, dcg -> threshold));
      }
    }
  }
  BirthdayIndex4 pop_pivot(priority_queue<BirthdayIndex4, vector<BirthdayIndex4>, BirthdayIndex4Comparator>& column)
  {
    if (column.empty()) {
      return BirthdayIndex4(0, -1, 0);
    } else {
      auto pivot = column.top();
      column.pop();
      
      while (!column.empty() && column.top().index == pivot.getIndex()) {
        column.pop();
        if (column.empty())
          return BirthdayIndex4(0, -1, 0);
        else {
          pivot = column.top();
          column.pop();
        }
      }
      return pivot;
    }
  }
  BirthdayIndex4 get_pivot(priority_queue<BirthdayIndex4, vector<BirthdayIndex4>, BirthdayIndex4Comparator>& column)
  {
    BirthdayIndex4 result = pop_pivot(column);
    if (result.getIndex() != -1)
      column.push(result);
    return result;
  }
  void assemble_columns_to_reduce()
  {
    ++dim;
    ctr -> dim = dim;
    
    if (dim == 1) { 
      ctr -> columns_to_reduce.clear();
      for(int w = 1; w <= aw; ++w){
        for(int z = 1; z <= az; ++z){
          for (int y = 1; y <= ay; ++y) {
            for (int x = 1; x <= ax; ++x) {
              for (int m = 0; m < 4; ++m) { // the number of type
                double index = x | (y << EXPONENT) | (z << (2 * EXPONENT)) | (w << (3 * EXPONENT)) | (m << (4 * EXPONENT));
                if (pivot_column_index.find(index) == pivot_column_index.end()) {
                  double birthday = dcg -> getBirthday(index, 1);
                  if (birthday != dcg -> threshold) {
                    ctr->columns_to_reduce.push_back(BirthdayIndex4(birthday, index, 1));
                  }
                }
              }
            }
          }
        }
      }
    } else if(dim == 2){ 
      ctr -> columns_to_reduce.clear();
      for(int w = 1; w <= aw; ++w){
        for(int z = 1; z <= az; ++z){
          for (int y = 1; y <= ay; ++y) {
            for (int x = 1; x <= ax; ++x) {
              for (int m = 0; m < 6; ++m) { // the number of type
                double index = x | (y << EXPONENT) | (z << (2 * EXPONENT)) | (w << (3 * EXPONENT)) | (m << (4 * EXPONENT));
                if (pivot_column_index.find(index) == pivot_column_index.end()) {
                  double birthday = dcg -> getBirthday(index, 2);
                  if (birthday != dcg -> threshold) {
                    ctr->columns_to_reduce.push_back(BirthdayIndex4(birthday, index, 2));
                  }
                }
              }
            }
          }
        }
      }
    } else if(dim == 3){
      ctr -> columns_to_reduce.clear();
      for(int w = 1; w <= aw; ++w){
        for(int z = 1; z <= az; ++z){
          for (int y = 1; y <= ay; ++y) {
            for (int x = 1; x <= ax; ++x) {
              for (int m = 0; m < 4; ++m) { // the number of type
                double index = x | (y << EXPONENT) | (z << (2 * EXPONENT)) | (w << (3 * EXPONENT)) | (m << (4 * EXPONENT));
                if (pivot_column_index.find(index) == pivot_column_index.end()) {
                  double birthday = dcg -> getBirthday(index, 3);
                  if (birthday != dcg -> threshold) {
                    ctr -> columns_to_reduce.push_back(BirthdayIndex4(birthday, index, 3));
                  }
                }
              }
            }
          }
        }
      }
    }
    sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndex4Comparator());
  }
};

// method == 0 --> LINKFIND algorithm
// method == 1 --> COMPUTEPAIRS algorithm
// [[Rcpp::export]]
Rcpp::NumericMatrix cubical_4dim(Rcpp::NumericVector& image, double threshold, int method, int nx, int ny, int nz, int nt)
{
  vector<WritePairs4> writepairs; // dim birth death
  writepairs.clear();
  
  DenseCubicalGrids4* dcg = new DenseCubicalGrids4(image, threshold, nx, ny, nz, nt);
  ColumnsToReduce4* ctr = new ColumnsToReduce4(dcg); 
  
  switch(method){
  case 0:
  {
    JointPairs4* jp = new JointPairs4(dcg, ctr, writepairs);
    jp -> joint_pairs_main();
    
    ComputePairs4* cp = new ComputePairs4(dcg, ctr, writepairs);
    cp -> compute_pairs_main(); // dim1
    
    cp -> assemble_columns_to_reduce();
    cp -> compute_pairs_main(); // dim2
    
    cp -> assemble_columns_to_reduce();
    cp -> compute_pairs_main(); // dim3
    
    // free pointers
    delete jp;
    delete cp;
    
    break;		
  }
  case 1:
  {	
    ComputePairs4* cp = new ComputePairs4(dcg, ctr, writepairs);
    cp -> compute_pairs_main(); // dim0
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim1
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim2
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim3
    
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
