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

#define FILE_OUTPUT

#include <iostream>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>
#include <queue>

using namespace std;

/*****birthday_index*****/
class BirthdayIndex
{
  
public:
  double birthday;
  int index;
  int dim;
  
  BirthdayIndex();
  
  BirthdayIndex(double _b, int _i, int _d);
  
  BirthdayIndex(const BirthdayIndex& b);
  
  void copyBirthdayIndex(BirthdayIndex v);
  
  double getBirthday();
  
  long getIndex();
  
  int getDimension();
  
  void print();
  
  void VertexPrint();
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

BirthdayIndex::BirthdayIndex(double _b, int _i, int _d){
  birthday = _b;
  index = _i;
  dim = _d;
}

BirthdayIndex::BirthdayIndex(const BirthdayIndex& b){
  birthday = b.birthday;
  index = b.index;
  dim = b.dim;
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
  std::cout << "(dob:" << birthday << "," << index << ")" << std::endl;
}

void BirthdayIndex::VertexPrint(){
  int px = index & 0x01ff;
  int py = (index >> 9) & 0x01ff;
  int pz = (index >> 18) & 0x01ff;
  int pm = (index >> 27) & 0xff;
  
  cout << "birthday : (m, z, y, x) = " << birthday << " : (" << pm << ", " << pz << ", " << py << ", " << px << ")" << endl; 
}

bool BirthdayIndexComparator::operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const{
  if(o1.birthday == o2.birthday){
    if(o1.index < o2.index){
      return true;
    } else {
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

/*****coeff*****/
class Coeff
{
public:
  int cx, cy, cz, cm;
  
  Coeff();
  
  void setXYZ(int _cx, int _cy, int _cz);
  
  void setXYZM(int _cx, int _cy, int _cz, int _cm);
  
  void setIndex(int index);
  
  int getIndex();
  
};

Coeff::Coeff(){
  cx = 0;
  cy = 0;
  cz = 0;
  cm = 0;
}

void Coeff::setXYZ(int _cx, int _cy, int _cz){
  cx = _cx;
  cy = _cy;
  cz = _cz;
  cm = 0;
}

void Coeff::setXYZM(int _cx, int _cy, int _cz, int _cm){
  cx = _cx;
  cy = _cy;
  cz = _cz;
  cm = _cm;
}

void Coeff::setIndex(int index){
  cx = index & 0x01ff;
  cy = (index >> 9) & 0x01ff;
  cz = (index >> 18) & 0x01ff;
  cm = (index >> 27) & 0xff;
}

int Coeff::getIndex(){
  return cx | cy << 9 | cz << 18 | cm << 27;
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

/*****vertices*****/
class Vertices
{
public:	
  Coeff** vertex;
  int dim; 
  int ox, oy, oz;
  int type;
  
  Vertices();
  
  void setVertices(int _dim, int _ox, int _oy, int _oz, int _om);
  
};

Vertices::Vertices(){
  dim = 0;
  vertex = new Coeff*[8];
  for(int d = 0; d < 8; ++d){
    vertex[d] = new Coeff();
  }
}

void Vertices::setVertices(int _dim, int _ox, int _oy, int _oz, int _om){ // 0 cell
  dim = _dim;
  ox = _ox;
  oy = _oy;
  oz = _oz;
  type = _om;
  
  if(dim == 0){
    vertex[0] -> setXYZ(_ox, _oy, _oz);
  } else if(dim == 1){
    switch(_om){
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
    
  } else if(dim == 2){
    switch(_om){
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
    
  } else if(dim == 3){ // cube
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

/*****dense_cubical_grids*****/
enum file_format { DIPHA, PERSEUS };

class DenseCubicalGrids { // file_read
public:
  double threshold;
  int dim;
  int ax, ay, az;
  double dense3[512][512][512];
  file_format format;
  
  DenseCubicalGrids(const std::string& filename, double _threshold, file_format _format);
  
  double getBirthday(int index, int dim);
  
  void GetSimplexVertices(int index, int dim, Vertices* v);
  
};

DenseCubicalGrids::DenseCubicalGrids(const string& filename, double _threshold, file_format _format)  {
  
  threshold = _threshold;
  format = _format;
  
  switch(format){
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
    for(int z = 0; z < az + 2; ++z){
      for (int y = 0; y < ay + 2; ++y) {
        for (int x = 0; x < ax + 2; ++x) {
          if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az){
            if (!fin.eof()) {
              fin.read( ( char * ) &dou, sizeof( double ) );
              dense3[x][y][z] = dou;
            } else {
              cerr << "file endof error " << endl;
            }
          }
          else {
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
    
    for(int z = 0; z < az + 2; ++z){
      for (int y = 0; y <ay + 2; ++y) { 
        for (int x = 0; x < ax + 2; ++x) { 
          if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az){ 
            if (!reading_file.eof()) { 
              getline(reading_file, reading_line_buffer); 
              dense3[x][y][z] = atoi(reading_line_buffer.c_str()); 
              if (dense3[x][y][z] == -1) { 
                dense3[x][y][z] = threshold; 
              } 
            } 
          }
          else { 
            dense3[x][y][z] = threshold; 
          } 
        } 
      }
    }
    break;
  }
  }
}


double DenseCubicalGrids::getBirthday(int index, int dim){
  int cx = index & 0x01ff;
  int cy = (index >> 9) & 0x01ff;
  int cz = (index >> 18) & 0x01ff;
  int cm = (index >> 27) & 0xff;
  
  switch(dim){
  case 0:
    return dense3[cx][cy][cz];
  case 1:
    switch(cm){
    case 0:
      return max(dense3[cx][cy][cz], dense3[cx + 1][cy][cz]);
    case 1:
      return max(dense3[cx][cy][cz], dense3[cx][cy + 1][cz]);
    case 2:
      return max(dense3[cx][cy][cz], dense3[cx][cy][cz + 1]);
    }
  case 2:
    switch(cm){
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


void DenseCubicalGrids::GetSimplexVertices(int index, int dim, Vertices* v){
  int cx = index & 0x01ff;
  int cy = (index >> 9) & 0x01ff;
  int cz = (index >> 18) & 0x01ff;
  int cm = (index >> 27) & 0xff;
  
  v -> setVertices(dim ,cx, cy, cz , cm);
}

/*****union_find*****/
class UnionFind{
public:
  int max_of_index;
  vector<int> parent;
  vector<double> birthtime;
  vector<double> time_max;
  DenseCubicalGrids* dcg;
  
  UnionFind(int moi, DenseCubicalGrids* _dcg); // Thie "n" is the number of cubes.
  
  int find(int x); // Thie "x" is Index.
  
  void link(int x, int y);
};

UnionFind::UnionFind(int moi, DenseCubicalGrids* _dcg) : parent(moi), birthtime(moi), time_max(moi) { // Thie "n" is the number of cubes.
  dcg = _dcg;
  max_of_index = moi;
  
  for(int i = 0; i < moi; ++i){
    parent[i] = i;
    birthtime[i] = dcg -> getBirthday(i, 0);
    time_max[i] = dcg -> getBirthday(i, 0);
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

/*****columns_to_reduce*****/
class ColumnsToReduce{
public:
  
  vector<BirthdayIndex> columns_to_reduce;
  int dim;
  int max_of_index;
  
  ColumnsToReduce(DenseCubicalGrids* _dcg); 
  
  int size(); 
};

ColumnsToReduce::ColumnsToReduce(DenseCubicalGrids* _dcg) { 
  dim = 0;
  int ax = _dcg -> ax;
  int ay = _dcg -> ay;
  int az = _dcg -> az;
  max_of_index = 512 * 512 * (az + 2);
  int index;
  double birthday;
  
  for(int z = az; z > 0; --z){
    for (int y = ay; y > 0; --y) {
      for (int x = ax; x > 0; --x) {
        birthday = _dcg -> dense3[x][y][z];
        index = x | (y << 9) | (z << 18);
        if (birthday != _dcg -> threshold) {
          columns_to_reduce.push_back(BirthdayIndex(birthday, index, 0));
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
  Vertices* vtx;
  double birthtime;
  int ax, ay, az;
  int cx, cy, cz;
  int count;
  BirthdayIndex nextCoface;
  double threshold;
  
  SimplexCoboundaryEnumerator();
  
  void setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg); 
  
  bool hasNextCoface(); 
  
  BirthdayIndex getNextCoface();
};

SimplexCoboundaryEnumerator::SimplexCoboundaryEnumerator(){
  vtx = new Vertices();
  nextCoface = BirthdayIndex(0, -1, 1);
}

void SimplexCoboundaryEnumerator::setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg) {
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

bool SimplexCoboundaryEnumerator::hasNextCoface() {
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
        nextCoface = BirthdayIndex(birthday, index, 1);
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
          nextCoface = BirthdayIndex(birthday, index, 2);
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
          nextCoface = BirthdayIndex(birthday, index, 2);
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
          nextCoface = BirthdayIndex(birthday, index, 2);
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
          nextCoface = BirthdayIndex(birthday, index, 3);
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
          nextCoface = BirthdayIndex(birthday, index, 3);
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
          nextCoface = BirthdayIndex(birthday, index, 3);
          return true;
        }
      }
      return false;
    }
  return false;
  }
}

BirthdayIndex SimplexCoboundaryEnumerator::getNextCoface() {
  return nextCoface;
}

/*****joint_pairs*****/
class JointPairs{
  
  int n; // the number of cubes
  int ctr_moi;
  int ax, ay, az;
  DenseCubicalGrids* dcg;
  ColumnsToReduce* ctr;
  vector<WritePairs> *wp;
  bool print;
  Vertices* vtx;
  double u, v;
  vector<int64_t> cubes_edges;
  vector<BirthdayIndex> dim1_simplex_list;
  
public:
  JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print);
  
  void joint_pairs_main();
};

JointPairs::JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print){
  dcg = _dcg;
  ax = dcg -> ax;
  ay = dcg -> ay;
  az = dcg -> az;
  ctr = _ctr; // ctr is "0-dim"simplex list.
  ctr_moi = ctr -> max_of_index;
  n = ctr -> columns_to_reduce.size();
  print = _print;
  
  wp = &_wp;
  vtx = new Vertices();
  
  for(int x = 1; x <= ax; ++x){
    for(int y = 1; y <= ay; ++y){
      for(int z = 1; z <= az; ++z){
        for(int type = 0; type < 3; ++type){
          int index = x | (y << 9) | (z << 18) | (type << 27);
          double birthday = dcg -> getBirthday(index, 1);
          if(birthday < dcg -> threshold){
            dim1_simplex_list.push_back(BirthdayIndex(birthday, index, 1));
          }
        }
      }
    }
  }
  sort(dim1_simplex_list.rbegin(), dim1_simplex_list.rend(), BirthdayIndexComparator());
}

void JointPairs::joint_pairs_main(){
  cubes_edges.reserve(2);
  UnionFind dset(ctr_moi, dcg);
  ctr -> columns_to_reduce.clear();
  ctr -> dim = 1;
  double min_birth = dcg -> threshold;
  
  if(print == true){
    cout << "persistence intervals in dim " << 0 << ":" << endl;
  }
  
  for(auto e : dim1_simplex_list){
    cubes_edges.clear();
    dcg -> GetSimplexVertices(e.getIndex(), 1, vtx);
    
    cubes_edges[0] = vtx -> vertex[0] -> getIndex();
    cubes_edges[1] = vtx -> vertex[1] -> getIndex();
    
    u = dset.find(cubes_edges[0]);
    v = dset.find(cubes_edges[1]);
    
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
  
  wp -> push_back(WritePairs(-1, min_birth, dcg -> threshold));
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
  int ax, ay, az;
  int dim;
  vector<WritePairs> *wp;
  bool print;
  
  ComputePairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print);
  
  void compute_pairs_main();
  
  void outputPP(int _dim, double _birth, double _death);
  
  BirthdayIndex pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
    column);
  
  BirthdayIndex get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
    column);
  
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
}

void ComputePairs::compute_pairs_main(){
  if(print == true){
    cout << "persistence intervals in dim " << dim << ":" << endl;
  }
  
  pivot_column_index = hash_map<int, int>();
  vector<BirthdayIndex> coface_entries;
  auto ctl_size = ctr -> columns_to_reduce.size();
  SimplexCoboundaryEnumerator cofaces;
  unordered_map<int, priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>> recorded_wc;
  
  pivot_column_index.reserve(ctl_size);
  recorded_wc.reserve(ctl_size);
  
  for(int i = 0; i < ctl_size; ++i){ 
    auto column_to_reduce = ctr -> columns_to_reduce[i]; 
    priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator> 
      working_coboundary;
    double birth = column_to_reduce.getBirthday();
    
    int j = i;
    BirthdayIndex pivot(0, -1, 0);
    bool might_be_apparent_pair = true;
    bool goto_found_persistence_pair = false;
    
    do {
      auto simplex = ctr -> columns_to_reduce[j]; // get CTR[i] 
      coface_entries.clear();
      cofaces.setSimplexCoboundaryEnumerator(simplex, dcg);// make coface data
      
      while (cofaces.hasNextCoface() && !goto_found_persistence_pair) { // repeat there remains a coface
        BirthdayIndex coface = cofaces.getNextCoface();
        coface_entries.push_back(coface);
        if (might_be_apparent_pair && (simplex.getBirthday() == coface.getBirthday())) { // If bt is the same, go thru
          if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) { // If coface is not in pivot list
            pivot.copyBirthdayIndex(coface); // I have a new pivot
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
            // I output PP as Writepairs
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

BirthdayIndex ComputePairs::pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
  column){
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
    for(int z = 1; z <= az; ++z){
      for (int y = 1; y <= ay; ++y) {
        for (int x = 1; x <= ax; ++x) {
          for (int m = 0; m < 3; ++m) { // the number of type
            double index = x | (y << 9) | (z << 18) | (m << 27);
            if (pivot_column_index.find(index) == pivot_column_index.end()) {
              double birthday = dcg -> getBirthday(index, 1);
              if (birthday != dcg -> threshold) {
                ctr -> columns_to_reduce.push_back(BirthdayIndex(birthday, index, 1));
              }
            }
          }
        }
      }
    }
  } else if(dim == 2){ 
    ctr -> columns_to_reduce.clear();
    for(int z = 1; z <= az; ++z){
      for (int y = 1; y <= ay; ++y) {
        for (int x = 1; x <= ax; ++x) {
          for (int m = 0; m < 3; ++m) { // the number of type
            double index = x | (y << 9) | (z << 18) | (m << 27);
            if (pivot_column_index.find(index) == pivot_column_index.end()) {
              double birthday = dcg -> getBirthday(index, 2);
              if (birthday != dcg -> threshold) {
                ctr -> columns_to_reduce.push_back(BirthdayIndex(birthday, index, 2));
              }
            }
          }
        }
      }
    }
  }
  sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator());
}

enum calculation_method { LINKFIND, COMPUTEPAIRS};

void print_usage_and_exit(int exit_code) {
  cerr << "Usage: "
       << "CR3 "
       << "[options] [input_filename]" << endl
       << endl
       << "Options:" << endl
       << endl
       << "  --help           print this screen" << endl
       << "  --format         use the specified file format for the input. Options are:" << endl
       << "                     dipha          (voxel matrix in DIPHA file format; default)" << endl
       << "                     perseus        (voxel matrix in Perseus file format)" << endl
       << "  --threshold <t>  compute cubical complexes up to birth time <t>" << endl
       << "  --method         method to compute the persistent homology of the cubical complexes. Options are" << endl
       << "                     link_find      (calculating the 0-dim PP, use 'link_find' algorithm; default)" << endl
       << "                     compute_pairs  (calculating the 0-dim PP, use 'compute_pairs' algorithm)" << endl
       << "  --output         name of file that will contain the persistence diagram " << endl
       << "  --print          print persistence pairs on your console" << endl
       << endl;
  
  exit(exit_code);
}

int main(int argc, char** argv){
  
  const char* filename = nullptr;
  string output_filename = "answer_3dim.diagram"; //default name
  file_format format = PERSEUS;
  calculation_method method = LINKFIND;
  double threshold = 99999;
  bool print = false;
  
  for (int i = 1; i < argc; ++i) {
    const string arg(argv[i]);
    if (arg == "--help") {
      print_usage_and_exit(0);
    } else if (arg == "--threshold") {
      string parameter = string(argv[++i]);
      size_t next_pos;
      threshold = stod(parameter, &next_pos);
      if (next_pos != parameter.size()) print_usage_and_exit(-1);
    } else if (arg == "--format") {
      string parameter = string(argv[++i]);
      if (parameter == "dipha") {
        format = DIPHA;
      } else if (parameter == "perseus") {
        format = PERSEUS;
      } else {
        print_usage_and_exit(-1);
      }
    } else if(arg == "--method") {
      string parameter = string(argv[++i]);
      if (parameter == "link_find") {
        method = LINKFIND;
      } else if (parameter == "compute_pairs") {
        method = COMPUTEPAIRS;
      } else {
        print_usage_and_exit(-1);
      }
    } else if (arg == "--output") {
      output_filename = string(argv[++i]);
    } else if (arg == "--print"){
      print = true;
    } else {
      if (filename) { print_usage_and_exit(-1); }
      filename = argv[i];
    }
  }
  
  ifstream file_stream(filename);
  if (filename && file_stream.fail()) {
    cerr << "couldn't open file " << filename << endl;
    exit(-1);
  }
  
  vector<WritePairs> writepairs; // dim birth death
  writepairs.clear();
  
  DenseCubicalGrids* dcg = new DenseCubicalGrids(filename, threshold, format);
  ColumnsToReduce* ctr = new ColumnsToReduce(dcg);
  
  switch(method){
  case LINKFIND:
  {
    JointPairs* jp = new JointPairs(dcg, ctr, writepairs, print);
    jp -> joint_pairs_main(); // dim0
    
    ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, print);
    cp -> compute_pairs_main(); // dim1
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim2
    break;
  }
    
  case COMPUTEPAIRS:
  {
    ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, print);
    cp -> compute_pairs_main(); // dim0
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim1
    cp -> assemble_columns_to_reduce();
    
    cp -> compute_pairs_main(); // dim2
    break;
  }
  }
  
#ifdef FILE_OUTPUT
  ofstream writing_file;
  
  string extension = ".csv";
  if(equal(extension.rbegin(), extension.rend(), output_filename.rbegin()) == true){
    
    string outname = output_filename;// .csv file
    writing_file.open(outname, ios::out);
    if(!writing_file.is_open()){
      cerr << " error: open file for output failed! " << endl;
    }
    
    int64_t p = writepairs.size();
    for(int64_t i = 0; i < p; ++i){
      writing_file << writepairs[i].getDimension() << ",";
      
      writing_file << writepairs[i].getBirth() << ",";
      writing_file << writepairs[i].getDeath() << endl;
    }
    writing_file.close();
  } else {
    
    writing_file.open(output_filename, ios::out | ios::binary);
    
    if(!writing_file.is_open()){
      cerr << " error: open file for output failed! " << endl;
    }
    
    int64_t mn = 8067171840;
    writing_file.write((char *) &mn, sizeof( int64_t )); // magic number
    int64_t type = 2;
    writing_file.write((char *) &type, sizeof( int64_t )); // type number of PERSISTENCE_DIAGRAM
    int64_t p = writepairs.size();
    cout << "the number of pairs : " << p << endl;
    writing_file.write((char *) &p, sizeof( int64_t )); // number of points in the diagram p
    for(int64_t i = 0; i < p; ++i){
      int64_t writedim = writepairs[i].getDimension();
      writing_file.write((char *) &writedim, sizeof( int64_t )); // dim
      
      double writebirth = writepairs[i].getBirth();
      writing_file.write((char *) &writebirth, sizeof( double )); // birth
      
      double writedeath = writepairs[i].getDeath();
      writing_file.write((char *) &writedeath, sizeof( double )); // death
    }
    writing_file.close();
  }
#endif
  
  return 0;
}
