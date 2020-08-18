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
  
public:
  double birthday;
  int index;
  int dim;
  
  BirthdayIndex();
  BirthdayIndex(const BirthdayIndex& b);
  BirthdayIndex(double _b, int _i, int _d);
  void copyBirthdayIndex(BirthdayIndex v);
  double getBirthday();
  long getIndex();
  int getDimension();
  void print();
};

struct BirthdayIndexComparator
{
  bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const;
};

struct BirthdayIndexInverseComparator
{
  bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const;
};

BirthdayIndex::BirthdayIndex(){
  birthday = 0;
  index = -1;
  dim = 1;
}

BirthdayIndex::BirthdayIndex(const BirthdayIndex& b){
  birthday = b.birthday;
  index = b.index;
  dim = b.dim;
}

BirthdayIndex::BirthdayIndex(double _b, int _i, int _d){
  birthday = _b;
  index = _i;
  dim = _d;
}

void BirthdayIndex::copyBirthdayIndex(BirthdayIndex v){
  birthday = v.birthday;
  index = v.index;
  dim = v.dim;
}

double BirthdayIndex::getBirthday(){
  return birthday;
}

long BirthdayIndex::getIndex(){
  return index;
}

int BirthdayIndex::getDimension(){
  return dim;
}

void BirthdayIndex::print(){
  Rcpp::Rcout << "(dob:" << birthday << "," << index << ")" << std::endl;
}

bool BirthdayIndexComparator::operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const{
  if(o1.birthday == o2.birthday){
    if(o1.index < o2.index){
      return true;
    } else{
      return false;
    }
  } else {
    if(o1.birthday > o2.birthday){
      return true;
    } else {
      return false;
    }
  }
}

bool BirthdayIndexInverseComparator::operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const{
  if(o1.birthday == o2.birthday){
    if(o1.index < o2.index){
      return false;
    } else {
      return true;
    }
  } else {
    if(o1.birthday > o2.birthday){
      return false;
    } else {
      return true;
    }
  }
}

/*****write_pairs*****/
class WritePairs
{
  
public:
  int64_t dim;
  double birth;
  double death;
  
  WritePairs(int64_t _dim, double _birth, double _death);
  int64_t getDimension();
  double getBirth();
  double getDeath();
};

WritePairs::WritePairs(int64_t _dim, double _birth, double _death){
  dim = _dim;
  birth = _birth;
  death = _death;
}

int64_t WritePairs::getDimension(){
  return dim;
}

double WritePairs::getBirth(){
  return birth;
}

double WritePairs::getDeath(){
  return death;
}

/*****dense_cubical_grids*****/
enum file_format { DIPHA, PERSEUS };

class DenseCubicalGrids { // file_read
  
public:
  double threshold;
  int dim;
  int ax, ay, az, aw;
  double dense4[128][128][128][128];
  file_format format;
  
  DenseCubicalGrids(const Rcpp::NumericVector& image, double threshold, int nx, int ny, int nz, int nt) ; 
  double getBirthday(int index, int dim);
};

DenseCubicalGrids::DenseCubicalGrids(const Rcpp::NumericVector& image, double threshold, int nx, int ny, int nz, int nt)
{
  ax = nx; ay = ny; az = nz; aw = nt;
  Rcpp::Rcout << "HERE4\n";
  // convert NumericVector to 4d array
  double ****dat = new double***[nx];
  Rcpp::Rcout << "HERE5\n";
  for (int i = 0; i < nx; i++)
  {
    dat[i] = new double**[ny];
    for (int j = 0; j < ny; j++)
    {
      dat[i][j] = new double*[nz];
      for (int k = 0; k < nz; k++)
        dat[i][j][k] = new double[nt];
    }
  }
  Rcpp::Rcout << "HERE6\n";
  int tempX, tempY, tempZ, tempT;
  for (int i = 0; i < nx * ny * nz * nt; i++)
  {
    tempX = i % nx;
    tempY = i / nx % ny;
    tempZ = i / (nx * ny) % nz;
    tempT = i / (nx * ny * nz) % nt;
    dat[tempX][tempY][tempZ][tempT] = image(i);
  }
  
  Rcpp::Rcout << "HERE7\n";
  for (int x = 0; x <= nx + 1; x++)
    for (int y = 0; y <= ny + 1; y++)
      for (int z = 0; z <= nz + 1; z++)
        for (int t = 0; t <= nt + 1; t++)
        {
          if (0 < x && x <= nx &&
              0 < y && y <= ny &&
              0 < z && z <= nz &&
              0 < t && t <= nt)
            dense4[x][y][z][t] = dat[x - 1][y - 1][z - 1][t - 1];
          else
            dense4[x][y][z][t] = threshold;
          
          // debug
          Rcpp::Rcout << "(" << x << "," << y << "," << z << "," << t << ") = " << dense4[x][y][z][t] << "\n";
        }
        
        for (int i = 0; i < nx; i++)
        {
          for (int j = 0; j < ny; j++)
          {
            for (int k = 0; k < nz; k++)
            {
              delete[] dat[i][j][k];
            }
            delete[] dat[i][j];
          }
          delete[] dat[i];
        }
        delete[] dat;
}


double DenseCubicalGrids::getBirthday(int index, int dim){
  int cx = index & 0x7f;
  int cy = (index >> 7) & 0x7f;
  int cz = (index >> 14) & 0x7f;
  int cw = (index >> 21) & 0x7f;
  int cm = (index >> 28) & 0x0f;
  
  switch(dim){
  case 0:
    return dense4[cx][cy][cz][cw];
  case 1:
    switch(cm){
    case 0:
      return max(dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw]);
    case 1:
      return max(dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw]);
    case 2:
      return max(dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw]);
    default:
      return max(dense4[cx][cy][cz][cw], dense4[cx][cy][cz][cw + 1]);
    break;
    }
  case 2:
    switch(cm){
    case 0: // x - y (fix z, w)
      return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
                 dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw]});
    case 1: // z - x (fix y, w)
      return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw], 
                 dense4[cx + 1][cy][cz + 1][cw], dense4[cx + 1][cy][cz][cw]});
    case 2: // y - z (fix x, w)
      return max({dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw], 
                 dense4[cx][cy + 1][cz + 1][cw], dense4[cx][cy][cz + 1][cw]});
    case 3: // x - w
      return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
                 dense4[cx + 1][cy][cz][cw + 1], dense4[cx][cy][cz][cw + 1]});
    case 4: // y - w
      return max({dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw], 
                 dense4[cx][cy + 1][cz][cw + 1], dense4[cx][cy][cz][cw + 1]});
    case 5: // z - w
      return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw], 
                 dense4[cx][cy][cz + 1][cw + 1], dense4[cx][cy][cz][cw + 1]});
    }
  case 3:
    switch(cm){
    case 0: // x - y - z
      return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
                 dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
                                                                       dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], 
                                                                                                                             dense4[cx + 1][cy + 1][cz + 1][cw], dense4[cx][cy + 1][cz + 1][cw]
      });
    case 1: // x - y - w
      return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
                 dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
                                                                       dense4[cx][cy][cz][cw + 1], dense4[cx + 1][cy][cz][cw + 1], 
                                                                                                                         dense4[cx + 1][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1]
      });
    case 2: // x - z - w
      return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
                 dense4[cx + 1][cy][cz][cw + 1], dense4[cx][cy][cz][cw + 1],
                                                                   dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], 
                                                                                                                         dense4[cx + 1][cy][cz + 1][cw + 1], dense4[cx][cy][cz + 1][cw + 1]
      });
    case 3: // y - z - w
      return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz][cw + 1], 
                 dense4[cx][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw],
                                                                       dense4[cx][cy][cz + 1][cw], dense4[cx][cy][cz + 1][cw + 1], 
                                                                                                                         dense4[cx][cy + 1][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw]
      });
    }
  case 4:
    return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
               dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
                                                                     dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], 
                                                                                                                           dense4[cx + 1][cy + 1][cz + 1][cw], dense4[cx][cy + 1][cz + 1][cw],
                                                                                                                                                                                         dense4[cx][cy][cz][cw + 1], dense4[cx + 1][cy][cz][cw + 1], 
                                                                                                                                                                                                                                           dense4[cx + 1][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1],
                                                                                                                                                                                                                                                                                                     dense4[cx][cy][cz + 1][cw + 1], dense4[cx + 1][cy][cz + 1][cw + 1], 
                                                                                                                                                                                                                                                                                                                                                               dense4[cx + 1][cy + 1][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw + 1]
    });
  }
  return threshold;
}

/*****columns_to_reduce*****/
class ColumnsToReduce
{
  
public:
  vector<BirthdayIndex> columns_to_reduce;
  int dim;
  int max_of_index;
  
  ColumnsToReduce(DenseCubicalGrids* _dcg);
  int size() ;
  
};

ColumnsToReduce::ColumnsToReduce(DenseCubicalGrids* _dcg) {
  dim = 0;
  int ax = _dcg -> ax;
  int ay = _dcg -> ay;
  int az = _dcg -> az;
  int aw = _dcg -> aw;
  max_of_index = 128*128*128*(aw + 2);
  int index;
  double birthday;
  
  for(int w = aw; w > 0; --w){
    for(int z = az; z > 0; --z){
      for (int y = ay; y > 0; --y) {
        for (int x = ax; x > 0; --x) {
          birthday = _dcg -> dense4[x][y][z][w];
          index = x | (y << 7) | (z << 14) | (w << 21);
          if (birthday != _dcg -> threshold) {
            columns_to_reduce.push_back(BirthdayIndex(birthday, index, 0));
          }
        }
      }
    }
  }
  sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndexComparator());
}

int ColumnsToReduce::size() {
  return columns_to_reduce.size();
}

/*****simplex_coboundary_estimator*****/
class SimplexCoboundaryEnumerator
{
public:
  BirthdayIndex simplex;
  DenseCubicalGrids* dcg;
  int dim;
  double birthtime;
  int ax, ay, az, aw;
  int cx, cy, cz, cw, cm;
  int count;
  BirthdayIndex nextCoface;
  double threshold; 
  
  SimplexCoboundaryEnumerator();
  void setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg);
  bool hasNextCoface() ;
  BirthdayIndex getNextCoface();
};

SimplexCoboundaryEnumerator::SimplexCoboundaryEnumerator(){
  nextCoface = BirthdayIndex(0, -1, 1);
}


void SimplexCoboundaryEnumerator::setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg) {
  simplex = _s;
  dcg = _dcg;
  dim = simplex.dim;
  birthtime = simplex.birthday;
  ax = _dcg -> ax;
  ay = _dcg -> ay;
  az = _dcg -> az;
  aw = _dcg -> aw;
  
  cx = simplex.index & 0x7f;
  cy = (simplex.index >> 7) & 0x7f;
  cz = (simplex.index >> 14) & 0x7f;
  cw = (simplex.index >> 21) & 0x7f;
  cm = (simplex.index >> 28) & 0x0f;
  
  threshold = _dcg->threshold;
  count = 0;
}


bool SimplexCoboundaryEnumerator::hasNextCoface() {
  int index = 0;
  double birthday = 0;
  switch (dim) {
  case 0: // dim0
    for (int i = count; i < 8; ++i) {
      switch (i){
      case 0: // w +
        index = (3 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx][cy][cz][cw + 1]);
        break;
        
      case 1: // w -
        index = (3 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx][cy][cz][cw - 1]);
        break;
        
      case 2: // z +
        index = (2 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx][cy][cz + 1][cw]);
        break;
        
      case 3: // z -
        index = (2 << 28) | (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx][cy][cz - 1][cw]);
        break;
        
      case 4: // y +
        index = (1 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx][cy + 1][cz][cw]);
        break;
        
      case 5: // y -
        index = (1 << 28) | (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx][cy - 1][cz][cw]);
        break;
        
      case 6: // x +
        index = (0 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
        birthday = max(birthtime, dcg -> dense4[cx + 1][cy][cz][cw]);
        break;
        
      case 7: // x -
        index = (0 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
        birthday = max(birthtime, dcg -> dense4[cx - 1][cy][cz][cw]);
        break;
      }
      if (birthday != threshold) {
        count = i + 1;
        nextCoface = BirthdayIndex(birthday, index, 1);
        return true;
      }
    }
    return false;
  case 1:// dim1
    switch (cm) {
    case 0: // dim1 type0 (x-axis -> )
      for(int i = count; i < 6; ++i){
        switch(i){
        case 0: // x - w +
          index = (3 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1]});
          break;
          
        case 1: // x - w -
          index = (3 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1]});
          break;
          
        case 2: // x - z +
          index = (1 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw]});
          break;
          
        case 3: // x - z -
          index = (1 << 28) | (cw << 21) |((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw]});
          break;
          
        case 4: // x - y +
          index = (0 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw]});
          break;
          
        case 5: // x - y -
          index = (0 << 28) | (cw << 21) |(cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw]});
          break;
        }
        
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 2);
          return true;
        }
      }
      return false;
      
    case 1: // dim1 type1 (y-axis -> )
      for(int i = count; i < 6; ++i){
        switch(i){
        case 0: // y - w +
          index = (4 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy + 1][cz][cw + 1]});
          break;
          
        case 1: // y - w -
          index = (4 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy + 1][cz][cw - 1]});
          break;
          
        case 2: // y - z +
          index = (2 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw]});
          break;
          
        case 3: // y - z -
          index = (2 << 28) | (cw << 21) |((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy + 1][cz - 1][cw]});
          break;
          
        case 4: // y - x +
          index = (0 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw]});
          break;
          
        case 5: // y - x -
          index = (0 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 2);
          return true;
        }
      }
      return false;
      
    case 2: // dim1 type2 (z-axis -> )
      for(int i = count; i < 6; ++i){
        switch(i){
        case 0: // z - w +
          index = (5 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy][cz + 1][cw + 1]});
          break;
          
        case 1: // z - w -
          index = (5 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy][cz + 1][cw - 1]});
          break;
          
        case 2: // z - y +
          index = (2 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw]});
          break;
          
        case 3: // z - y -
          index = (2 << 28) | (cw << 21) |(cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz + 1][cw]});
          break;
          
        case 4: // z - x +
          index = (1 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw]});
          break;
          
        case 5: // z - x -
          index = (1 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz + 1][cw]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 2);
          return true;
        }
      }
      return false;
      
    case 3: // dim1 type3 (w-axis -> )
      for(int i = count; i < 6; ++i){
        switch(i){
        case 0: // w - z +
          index = (5 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy][cz + 1][cw + 1]});
          break;
          
        case 1: // w - z -
          index = (5 << 28) | (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy][cz - 1][cw + 1]});
          break;
          
        case 2: // w - y +
          index = (4 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz][cw + 1]});
          break;
          
        case 3: // w - y -
          index = (4 << 28) | (cw << 21) |(cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz][cw + 1]});
          break;
          
        case 4: // w - x +
          index = (3 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz][cw + 1]});
          break;
          
        case 5: // w - x -
          index = (3 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz][cw + 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 2);
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
          index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
                         dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
          break;
          
        case 1: // w -
          index = (1 << 28)| ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
                         dcg -> dense4[cx][cy + 1][cz][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz][cw - 1]});
          break;
          
        case 2: // z +
          index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                         dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
          break;
          
        case 3: // z -
          index = (0 << 28)| (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
                         dcg -> dense4[cx][cy + 1][cz - 1][cw],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw]});
          break;
          
        }
        
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 3);
          return true;
        }
      }
      return false;
      
    case 1: // dim2 type1 (fix x - z)
      for(int i = count; i < 4; ++i){
        switch(i){
        case 0: // w +
          index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
                         dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
          break;
          
        case 1: // w -
          index = (2 << 28)| ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
                         dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy][cz + 1][cw - 1]});
          break;
          
        case 2: // y +
          index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
          break;
          
        case 3: // y -
          index = (0 << 28)| (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
                         dcg -> dense4[cx][cy - 1][cz + 1][cw],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw]});
          break;
          
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 3);
          return true;
        }
      }
      return false;
      
    case 2: // dim2 type3 (fix y - z)
      for(int i = count; i < 4; ++i){
        switch(i){
        case 0: // w +
          index = (3 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy + 1][cz][cw + 1], 
                         dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // w -
          index = (3 << 28)| ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy + 1][cz][cw - 1], 
                         dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx][cy + 1][cz + 1][cw - 1]});
          break;
          
        case 2: // x +
          index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx + 1][cy][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
          break;
          
        case 3: // x -
          index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx - 1][cy][cz + 1][cw],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 3);
          return true;
        }
      }
      return false;
      
    case 3: // dim2 type2 (fix x - w)
      for(int i = count; i < 4; ++i){
        switch(i){
        case 0: // z +
          index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                         dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
          break;
          
        case 1: // z -
          index = (2 << 28)| (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
                         dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy][cz - 1][cw + 1]});
          break;
          
        case 2: // y +
          index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
          break;
          
        case 3: // y -
          index = (1 << 28)| (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
                         dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz][cw + 1]});
          break;
          
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 3);
          return true;
        }
      }
      return false;
      
    case 4: // dim2 type4 (fix y - w)
      for(int i = count; i < 4; ++i){
        switch(i){
        case 0: // z +
          index = (3 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw], 
                         dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // z -
          index = (3 << 28)| (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy + 1][cz - 1][cw], 
                         dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx][cy + 1][cz - 1][cw + 1]});
          break;
          
        case 2: // x +
          index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
          break;
          
        case 3: // x -
          index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz][cw + 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 3);
          return true;
        }
      }
      return false;
      
    case 5: // dim2 type5 (fix z - w)
      for(int i = count; i < 4; ++i){
        switch(i){
        case 0: // y +
          index = (3 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw], 
                         dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // y -
          index = (3 << 28)| (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz + 1][cw], 
                         dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx][cy - 1][cz + 1][cw + 1]});
          break;
          
        case 2: // x +
          index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                         dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
          break;
          
        case 3: // x -
          index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz + 1][cw], 
                         dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy][cz + 1][cw + 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 3);
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
          index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
                         dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
                                                                                                dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],
                                                                                                                                                                       dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // w -
          index = ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
                         dcg -> dense4[cx][cy + 1][cz][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz][cw - 1],
                                                                                                dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy][cz + 1][cw - 1],
                                                                                                                                                                       dcg -> dense4[cx][cy + 1][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw - 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 4);
          return true;
        }
      }
      return false;
      
    case 1: // dim3 type1 (x - y - w)
      for(int i = count; i < 2; ++i){
        switch(i){
        case 0: // z +
          index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
                         dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
                                                                                                    dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],
                                                                                                                                                                           dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // z -
          index = (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
                         dcg -> dense4[cx][cy + 1][cz - 1][cw],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw],
                                                                                                    dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy][cz - 1][cw + 1],
                                                                                                                                                                           dcg -> dense4[cx][cy + 1][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw + 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 4);
          return true;
        }
      }
      return false;
      
    case 2: // dim3 type2 (x - z - w)
      for(int i = count; i < 2; ++i){
        switch(i){
        case 0: // y +
          index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
                                                                                                    dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
                                                                                                                                                                           dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // y -
          index = (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
                         dcg -> dense4[cx][cy - 1][cz + 1][cw],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw],
                                                                                                    dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz][cw + 1],
                                                                                                                                                                           dcg -> dense4[cx][cy - 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw + 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 4);
          return true;
        }
      }
      return false;
      
    case 3: // dim3 type3 (y - z - w)
      for(int i = count; i < 2; ++i){
        switch(i){
        case 0: // x +
          index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
          birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx + 1][cy][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
                                                                                                    dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
                                                                                                                                                                           dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
          break;
          
        case 1: // x -
          index = (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
          birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
                         dcg -> dense4[cx - 1][cy][cz + 1][cw],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw],
                                                                                                    dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz][cw + 1],
                                                                                                                                                                           dcg -> dense4[cx - 1][cy][cz + 1][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw + 1]});
          break;
        }
        if (birthday != threshold) {
          count = i + 1;
          nextCoface = BirthdayIndex(birthday, index, 4);
          return true;
        }
      }
      return false;
    }
    return false;
  }
  return false;
}

BirthdayIndex SimplexCoboundaryEnumerator::getNextCoface() {
  return nextCoface;
}

/*****union_find*****/
class UnionFind{
public:
  int max_of_index;
  vector<int> parent;
  vector<double> birthtime;
  vector<double> time_max;
  DenseCubicalGrids* dcg;
  
  UnionFind(int moi, DenseCubicalGrids* _dcg);
  int find(int x);
  void link(int x, int y);
};


class JointPairs{
  
  int n; // the number of cubes
  int ctr_moi;
  int ax, ay, az, aw;
  DenseCubicalGrids* dcg;
  ColumnsToReduce* ctr;
  vector<WritePairs> *wp;
  bool print;
  double u, v;
  vector<int64_t> cubes_edges;
  vector<BirthdayIndex> dim1_simplex_list;
  
public:
  JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print);
  void joint_pairs_main();
};

UnionFind::UnionFind(int moi, DenseCubicalGrids* _dcg) : parent(moi), birthtime(moi), time_max(moi) { // Thie "n" is the number of cubes.
  dcg = _dcg;
  max_of_index = moi;
  
  for(int i = 0; i < moi; ++i){
    parent[i] = i;
    birthtime[i] = dcg->getBirthday(i, 0);
    time_max[i] = dcg->getBirthday(i, 0);
  }
}

int UnionFind::find(int x){ // Thie "x" is Index.
  int y = x, z = parent[y];
  while (z != y) {
    y = z;
    z = parent[y];
  }
  y = parent[x];
  while (z != y) {
    parent[x] = z;
    x = y;
    y = parent[x];
  }
  return z;
}

void UnionFind::link(int x, int y){
  x = find(x);
  y = find(y);
  if (x == y) return;
  if (birthtime[x] > birthtime[y]){
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


JointPairs::JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print){
  dcg = _dcg;
  ax = dcg -> ax;
  ay = dcg -> ay;
  az = dcg -> az;
  aw = dcg -> aw;
  ctr = _ctr; // ctr is "dim0"simplex list.
  ctr_moi = ctr -> max_of_index;
  n = ctr -> columns_to_reduce.size();
  print = _print;
  
  wp = &_wp;
  
  for(int x = 1; x <= ax; ++x){
    for(int y = 1; y <= ay; ++y){
      for(int z = 1; z <= az; ++z){
        for(int w = 1; w <= aw; ++w){
          for(int type = 0; type < 4; ++type){ // change
            int index = x | (y << 7) | (z << 14) | (w << 21) | (type << 28);
            double birthday = dcg -> getBirthday(index, 1);
            if(birthday < dcg -> threshold){
              dim1_simplex_list.push_back(BirthdayIndex(birthday, index, 1));
            }
          }
        }
      }
    }
  }
  sort(dim1_simplex_list.rbegin(), dim1_simplex_list.rend(), BirthdayIndexComparator());
}

void JointPairs::joint_pairs_main(){
  UnionFind dset(ctr_moi, dcg);
  ctr -> columns_to_reduce.clear();
  ctr -> dim = 1;
  double min_birth = dcg -> threshold;
  
  if(print == true){
    cout << "persistence intervals in dim " << 0 << ":" << endl;
  }
  
  for(BirthdayIndex e : dim1_simplex_list){
    int index = e.getIndex();
    int cx = index & 0x7f;
    int cy = (index >> 7) & 0x7f;
    int cz = (index >> 14) & 0x7f;
    int cw = (index >> 21) & 0x7f;
    int cm = (index >> 28) & 0x0f;
    
    int ce0=0, ce1 =0;
    switch(cm){
    case 0:
      ce0 =  cx | (cy << 7) | (cz << 14) | (cw << 21);
      ce1 =  (cx + 1) | (cy << 7) | (cz << 14) | (cw << 21);
      break;
    case 1:
      ce0 =  cx | (cy << 7) | (cz << 14) | (cw << 21);
      ce1 =  cx | ((cy + 1) << 7) | (cz << 14) | (cw << 21);
      break;
    case 2:
      ce0 =  cx | (cy << 7) | (cz << 14) | (cw << 21);
      ce1 =  cx | (cy << 7) | ((cz + 1) << 14) | (cw << 21);
      break;
    case 3:
      ce0 =  cx | (cy << 7) | (cz << 14) | (cw << 21);
      ce1 =  cx | (cy << 7) | (cz << 14) | ((cw + 1) << 21);
      break;
    }
    u = dset.find(ce0);
    v = dset.find(ce1);
    if(min_birth >= min(dset.birthtime[u], dset.birthtime[v])){
      min_birth = min(dset.birthtime[u], dset.birthtime[v]);
    }
    
    if(u != v){
      double birth = max(dset.birthtime[u], dset.birthtime[v]);
      double death = max(dset.time_max[u], dset.time_max[v]);
      
      if(birth == death){
        dset.link(u, v);
      } else {
        if(print == true){
          cout << "[" << birth << "," << death << ")" << endl;
        }
        
        wp -> push_back(WritePairs(0, birth, death));
        dset.link(u, v);
      }
    } else { // If two values have same "parent", these are potential edges which make a 2-simplex.
      ctr -> columns_to_reduce.push_back(e);
    }
  }
  
  if(print == true){
    cout << "[" << min_birth << ", )" << endl;
  }
  
  wp -> push_back(WritePairs(-1, min_birth, dcg->threshold));
  sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator());
}

/*****compute_pairs*****/
template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};


class ComputePairs
{
  
public:
  DenseCubicalGrids* dcg;
  ColumnsToReduce* ctr;
  hash_map<int, int> pivot_column_index;
  int ax, ay, az, aw;
  int dim;
  vector<WritePairs> *wp;
  bool print;
  
  ComputePairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print);
  void compute_pairs_main();
  void outputPP(int _dim, double _birth, double _death);
  BirthdayIndex pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>& column);
  BirthdayIndex get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>& column);
  void assemble_columns_to_reduce();
};

ComputePairs::ComputePairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print){
  dcg = _dcg;
  ctr = _ctr;
  dim = _ctr -> dim;
  wp = &_wp;
  print = _print;
  
  ax = _dcg -> ax;
  ay = _dcg -> ay;
  az = _dcg -> az;
  aw = _dcg -> aw;
}

void ComputePairs::compute_pairs_main(){
  if(print == true){
    cout << "persistence intervals in dim " << dim << ":" << endl;
  }
  
  vector<BirthdayIndex> coface_entries;
  SimplexCoboundaryEnumerator cofaces;
  unordered_map<int, priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>> recorded_wc;
  
  pivot_column_index = hash_map<int, int>();
  auto ctl_size = ctr -> columns_to_reduce.size();
  pivot_column_index.reserve(ctl_size);
  recorded_wc.reserve(ctl_size);
  
  
  for(int i = 0; i < ctl_size; ++i){ 
    auto column_to_reduce = ctr -> columns_to_reduce[i]; 
    priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator> working_coboundary;
    double birth = column_to_reduce.getBirthday();
    
    int j = i;
    BirthdayIndex pivot(0, -1, 0);
    bool might_be_apparent_pair = true;
    bool goto_found_persistence_pair = false;
    
    do {
      auto simplex = ctr -> columns_to_reduce[j];
      coface_entries.clear();
      cofaces.setSimplexCoboundaryEnumerator(simplex, dcg);
      
      while (cofaces.hasNextCoface() && !goto_found_persistence_pair) {
        BirthdayIndex coface = cofaces.getNextCoface();
        coface_entries.push_back(coface);
        if (might_be_apparent_pair && (simplex.getBirthday() == coface.getBirthday())) {
          if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) {
            pivot.copyBirthdayIndex(coface);
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

void ComputePairs::outputPP(int _dim, double _birth, double _death){
  if(_birth != _death){
    if(_death != dcg -> threshold){
      if(print == true){
        cout << "[" <<_birth << "," << _death << ")" << endl;
      }
      
      wp -> push_back(WritePairs(_dim, _birth, _death));
    } else {
      if(print == true){
        cout << "[" << _birth << ", )" << endl;
      }
      
      wp -> push_back(WritePairs(-1, _birth, dcg -> threshold));
    }
  }
}

BirthdayIndex ComputePairs::pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>& column){
  if (column.empty()) {
    return BirthdayIndex(0, -1, 0);
  } else {
    auto pivot = column.top();
    column.pop();
    
    while (!column.empty() && column.top().index == pivot.getIndex()) {
      column.pop();
      if (column.empty())
        return BirthdayIndex(0, -1, 0);
      else {
        pivot = column.top();
        column.pop();
      }
    }
    return pivot;
  }
}

BirthdayIndex ComputePairs::get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
  column) {
  BirthdayIndex result = pop_pivot(column);
  if (result.getIndex() != -1) {
    column.push(result);
  }
  return result;
}

void ComputePairs::assemble_columns_to_reduce() {
  ++dim;
  ctr -> dim = dim;
  
  if (dim == 1) { 
    ctr -> columns_to_reduce.clear();
    for(int w = 1; w <= aw; ++w){
      for(int z = 1; z <= az; ++z){
        for (int y = 1; y <= ay; ++y) {
          for (int x = 1; x <= ax; ++x) {
            for (int m = 0; m < 4; ++m) { // the number of type
              double index = x | (y << 7) | (z << 14) | (w << 21) | (m << 28);
              if (pivot_column_index.find(index) == pivot_column_index.end()) {
                double birthday = dcg -> getBirthday(index, 1);
                if (birthday != dcg -> threshold) {
                  ctr->columns_to_reduce.push_back(BirthdayIndex(birthday, index, 1));
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
              double index = x | (y << 7) | (z << 14) | (w << 21) | (m << 28);
              if (pivot_column_index.find(index) == pivot_column_index.end()) {
                double birthday = dcg -> getBirthday(index, 2);
                if (birthday != dcg -> threshold) {
                  ctr->columns_to_reduce.push_back(BirthdayIndex(birthday, index, 2));
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
              double index = x | (y << 7) | (z << 14) | (w << 21) | (m << 28);
              if (pivot_column_index.find(index) == pivot_column_index.end()) {
                double birthday = dcg -> getBirthday(index, 3);
                if (birthday != dcg -> threshold) {
                  ctr -> columns_to_reduce.push_back(BirthdayIndex(birthday, index, 3));
                }
              }
            }
          }
        }
      }
    }
  }
  sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator());
}

/*****cubical_ripser4dim*****/
// method == 0 --> link find algo
// method == 1 --> compute pairs algo
// 
Rcpp::NumericMatrix cubical_4dim(Rcpp::NumericVector& image, double threshold, int method, int nx, int ny, int nz, int nt)
{
  Rcpp::Rcout << "HERE3\n";
  bool print = false;
  
  vector<WritePairs> writepairs; // dim birth death
  writepairs.clear();
  
  DenseCubicalGrids* dcg = new DenseCubicalGrids(image, threshold, nx, ny, nz, nt);
  ColumnsToReduce* ctr = new ColumnsToReduce(dcg); 
  
  
  
  switch(method)
  {
    // link find algo
  case 0:
  {
    JointPairs* jp = new JointPairs(dcg, ctr, writepairs, print);
    jp -> joint_pairs_main();
    
    ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, print);
    cp -> compute_pairs_main(); // dim1
    
    cp -> assemble_columns_to_reduce();
    cp -> compute_pairs_main(); // dim2
    
    cp -> assemble_columns_to_reduce();
    cp -> compute_pairs_main(); // dim3
    break;		
  }
    // compute pairs algo
  case 1:
  {	
    ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, print);
    cp -> compute_pairs_main(); // dim0
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim1
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim2
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim3
    break;
  }
  }
  
  Rcpp::NumericMatrix ans(writepairs.size(), 3);
  for (int i = 0; i < ans.nrow(); i++)
  {
    ans(i, 0) = writepairs[i].getDimension();
    ans(i, 1) = writepairs[i].getBirth();
    ans(i, 2) = writepairs[i].getDeath();
    
    // debug
    Rcpp::Rcout << ans(i, 0) << " " << ans(i, 1) << " " << ans(i, 2) << "\n";
  }
  
  return ans;
}