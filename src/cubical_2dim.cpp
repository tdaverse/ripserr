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
class BirthdayIndex2
{
  //member vars
public:
  double birthday;
  int index;
  int dim;
  
  // constructors
  BirthdayIndex2(double _b, int _i, int _d) : birthday(_b), index(_i), dim(_d) {}
  BirthdayIndex2() : BirthdayIndex2(0, -1, 1) {}
  BirthdayIndex2(const BirthdayIndex2& b) : BirthdayIndex2(b.birthday, b.index, b.dim) {}

  // copy method
  void copyBirthdayIndex(BirthdayIndex2 v) { birthday = v.birthday; index = v.index; dim = v.dim; }
  
  // getters
  double getBirthday() { return birthday; }
  long getIndex() { return index; }
  int getDimension() { return dim; }
};

bool cmp(const BirthdayIndex2& o1, const BirthdayIndex2& o2) { return (o1.birthday == o2.birthday ? o1.index < o2.index : o1.birthday > o2.birthday); }

struct BirthdayIndex2Comparator
{
  bool operator()(const BirthdayIndex2& o1, const BirthdayIndex2& o2) const
  { return cmp(o1, o2); }
};

struct BirthdayIndex2InverseComparator
{
  bool operator()(const BirthdayIndex2& o1, const BirthdayIndex2& o2) const
  { return !cmp(o1, o2); }
};

/*****dense_cubical_grids*****/
class DenseCubicalGrids2 // file_read
{
public:
  double threshold;
  int dim;
  int ax, ay;
  double dense2[2048][1024];

  // constructor (w/ file read)
  DenseCubicalGrids2(const Rcpp::NumericMatrix& image, double _threshold) : threshold(_threshold), ax(image.nrow()), ay(image.ncol())
  {
    // assert that dimensions are not too big
    assert(0 < ax && ax < 2000 && 0 < ay && ay < 1000);

    // copy over data from NumericMatrix into DenseCubicalGrids member var
    for (int y = 0; y < ay + 2; y++)
      for (int x = 0; x < ax + 2; x++)
        if (0 < x && x <= ax && 0 < y && y <= ay) dense2[x][y] = image(x - 1, y - 1);
        else dense2[x][y] = threshold;
  }

  // getter
  double getBirthday(int index, int dim)
  {
    int cx = index & 0x07ff,
        cy = (index >> 11) & 0x03ff,
        cm = (index >> 21) & 0xff;

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
class WritePairs2
{
  // member vals
public:
  int64_t dim;
  double birth;
  double death;

  // constructor
  WritePairs2(int64_t _dim, double _birth, double _death) : dim(_dim), birth(_birth), death(_death) {}

  // getters
  int64_t getDimension() { return dim; }
  double getBirth() { return birth; }
  double getDeath() { return death; }
};

/*****columns_to_reduce*****/
class ColumnsToReduce2
{
  // member vars
public:
  vector<BirthdayIndex2> columns_to_reduce;
  int dim;
  int max_of_index;

  // constructor
  ColumnsToReduce2(DenseCubicalGrids2* _dcg) : dim(0)
  {
    int ax = _dcg->ax,
        ay = _dcg->ay,
        index;
    max_of_index = 2048 * (ay + 2);
    double birthday;
    for (int y = ay; y > 0; --y)
      for (int x = ax; x > 0; --x)
      {
        birthday = _dcg->dense2[x][y];
        index = x | (y << 11);
        if (birthday != _dcg -> threshold) columns_to_reduce.push_back(BirthdayIndex2(birthday, index, 0));
      }
    sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndex2Comparator());
  }

  // getter (length of member vector)
  int size() { return columns_to_reduce.size(); }
};

/*****simplex_coboundary_enumerator*****/
class SimplexCoboundaryEnumerator2
{
  // member vars
public:
  BirthdayIndex2 simplex;
  DenseCubicalGrids2* dcg;
  int dim;
  double birthtime;
  int ax, ay;
  int cx, cy, cm;
  int count;
  BirthdayIndex2 nextCoface;
  double threshold;

  // constructor
  SimplexCoboundaryEnumerator2() : nextCoface(BirthdayIndex2(0, -1, 1)) {}

  // member methods
  void setSimplexCoboundaryEnumerator2(BirthdayIndex2 _s, DenseCubicalGrids2* _dcg)
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
          nextCoface = BirthdayIndex2(birthday, index, 1);
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
            nextCoface = BirthdayIndex2(birthday, index, 2);
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
            nextCoface = BirthdayIndex2(birthday, index, 2);
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
            nextCoface = BirthdayIndex2(birthday, index, 2);
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
            nextCoface = BirthdayIndex2(birthday, index, 2);
            return true;
          }
        }
        return false;
      }
    }
    return false;
  }

  // getter
  BirthdayIndex2 getNextCoface() { return nextCoface; }
};

/*****union_find*****/
class UnionFind2
{
  // member vars
public:
  int max_of_index;
  vector<int> parent;
  vector<double> birthtime;
  vector<double> time_max;
  DenseCubicalGrids2* dcg;

  // constructor
  UnionFind2(int moi, DenseCubicalGrids2* _dcg) : max_of_index(moi) // Thie "n" is the number of cubes.
  {
    parent = vector<int>(moi);
    birthtime = vector<double>(moi);
    time_max = vector<double>(moi);
    
    dcg = _dcg;

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

/*****joint_pairs*****/
class JointPairs2
{
  int n; // the number of cubes
  int ctr_moi;
  int ax, ay;
  DenseCubicalGrids2* dcg;
  ColumnsToReduce2* ctr;
  vector<WritePairs2> *wp;
  bool print;
  double u, v;
  vector<int64_t> cubes_edges;
  vector<BirthdayIndex2> dim1_simplex_list;

public:
  // constructor
  JointPairs2(DenseCubicalGrids2* _dcg, ColumnsToReduce2* _ctr, vector<WritePairs2> &_wp, const bool _print)
  {
    dcg = _dcg;
    ax = dcg -> ax;
    ay = dcg -> ay;
    ctr = _ctr; // ctr is "dim0" simplex list.
    ctr_moi = ctr -> max_of_index;
    n = ctr -> columns_to_reduce.size();
    print = _print;

    wp = &_wp;

    for (int x = 1; x <= ax; ++x)
    {
      for (int y = 1; y <= ay; ++y)
      {
        for (int type = 0; type < 2; ++type)
        {
          int index = x | (y << 11) | (type << 21);
          double birthday = dcg -> getBirthday(index, 1);
          if (birthday < dcg -> threshold)
          {
            dim1_simplex_list.push_back(BirthdayIndex2(birthday, index, 1));
          }
        }
      }
    }

    sort(dim1_simplex_list.rbegin(), dim1_simplex_list.rend(), BirthdayIndex2Comparator());
  }

  // member method - workhorse
  void joint_pairs_main()
  {
    UnionFind2 dset(ctr_moi, dcg);
    ctr->columns_to_reduce.clear();
    ctr->dim = 1;
    double min_birth = dcg->threshold;

    for (BirthdayIndex2 e : dim1_simplex_list)
    {
      int index = e.getIndex();
      int cx = index & 0x07ff;
      int cy = (index >> 11) & 0x03ff;
      int cm = (index >> 21) & 0xff;
      int ce0=0, ce1 =0;

      switch (cm)
      {
      case 0:
        ce0 = ((cy) << 11) | cx;
        ce1 = ((cy) << 11) | (cx + 1);
        break;
      default:
        ce0 = ((cy) << 11) | cx;
      ce1 = ((cy + 1) << 11) | cx;
      break;
      }

      u = dset.find(ce0);
      v = dset.find(ce1);
      if (min_birth >= min(dset.birthtime[u], dset.birthtime[v]))
      {
        min_birth = min(dset.birthtime[u], dset.birthtime[v]);
      }

      if (u != v)
      {
        double birth = max(dset.birthtime[u], dset.birthtime[v]);
        double death = max(dset.time_max[u], dset.time_max[v]);
        if (birth == death)
        {
          dset.link(u, v);
        }
        else
        {
          wp->push_back(WritePairs2(0, birth, death));
          dset.link(u, v);
        }
      }
      else // If two values have same "parent", these are potential edges which make a 2-simplex.
      {
        ctr->columns_to_reduce.push_back(e);
      }
    }

    wp->push_back(WritePairs2(-1, min_birth, dcg->threshold));
    sort(ctr->columns_to_reduce.begin(), ctr->columns_to_reduce.end(), BirthdayIndex2Comparator());
  }
};

/*****compute_pairs*****/
template <class Key, class T> class hash_map2 : public unordered_map<Key, T> {};

class ComputePairs2
{
  //member vars
public:
  DenseCubicalGrids2* dcg;
  ColumnsToReduce2* ctr;
  hash_map2<int, int> pivot_column_index;
  int ax, ay;
  int dim;
  vector<WritePairs2> *wp;
  bool print;

  // constructor
  ComputePairs2(DenseCubicalGrids2* _dcg, ColumnsToReduce2* _ctr, vector<WritePairs2> &_wp, const bool _print)
  {
    dcg = _dcg;
    ctr = _ctr;
    dim = _ctr -> dim;
    wp = &_wp;
    print = _print;

    ax = _dcg -> ax;
    ay = _dcg -> ay;
  }

  // member methods
  //   workhorse
  void compute_pairs_main()
  {
    vector<BirthdayIndex2> coface_entries;
    SimplexCoboundaryEnumerator2 cofaces;
    unordered_map<int, priority_queue<BirthdayIndex2, vector<BirthdayIndex2>, BirthdayIndex2Comparator>> recorded_wc;

    pivot_column_index = hash_map2<int, int>();
    auto ctl_size = ctr->columns_to_reduce.size();
    pivot_column_index.reserve(ctl_size);
    recorded_wc.reserve(ctl_size);

    for (int i = 0; i < ctl_size; ++i)
    {
      if (i % 2500 == 0) {
        Rcpp::checkUserInterrupt();
      }
      
      auto column_to_reduce = ctr->columns_to_reduce[i];
      priority_queue<BirthdayIndex2, vector<BirthdayIndex2>, BirthdayIndex2Comparator> working_coboundary;
      double birth = column_to_reduce.getBirthday();

      int j = i;
      BirthdayIndex2 pivot(0, -1, 0);
      bool might_be_apparent_pair = true;
      bool goto_found_persistence_pair = false;

      do {
        auto simplex = ctr->columns_to_reduce[j];// get CTR[i]
        coface_entries.clear();
        cofaces.setSimplexCoboundaryEnumerator2(simplex, dcg);// make cofaces data

        while (cofaces.hasNextCoface() && !goto_found_persistence_pair) // repeat there remains a coface
        {
          BirthdayIndex2 coface = cofaces.getNextCoface();
          coface_entries.push_back(coface);
          if (might_be_apparent_pair && (simplex.getBirthday() == coface.getBirthday())) // if bt is the same, go thru
          {
            if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) // if coface is not in pivot list
            {
              pivot.copyBirthdayIndex(coface);// I have a new pivot
              goto_found_persistence_pair = true;// goto (B)
            }
            else // if pivot list contains this coface,
            {
              might_be_apparent_pair = false;// goto(A)
            }
          }
        }

        if (!goto_found_persistence_pair) // (A) if pivot list contains this coface
        {
          auto findWc = recorded_wc.find(j); // we seek wc list by 'j'
          if (findWc != recorded_wc.end()) // if the pivot is old,
          {
            auto wc = findWc->second;
            while (!wc.empty()) // we push the data of the old pivot's wc
            {
              auto e = wc.top();
              working_coboundary.push(e);
              wc.pop();
            }
          }
          else // if the pivot is new,
          {
            for (auto e : coface_entries) // making wc here
            {
              working_coboundary.push(e);
            }
          }
          pivot = get_pivot(working_coboundary); // getting a pivot from wc

          if (pivot.getIndex() != -1) //When I have a pivot, ...
          {
            auto pair = pivot_column_index.find(pivot.getIndex());
            if (pair != pivot_column_index.end()) // if the pivot already exists, go on the loop
            {
              j = pair->second;
              continue;
            }
            else // if the pivot is new,
            {
              // I record this wc into recorded_wc, and
              recorded_wc.insert(make_pair(i, working_coboundary));
              // I output PP as Writepairs
              double death = pivot.getBirthday();
              outputPP(dim, birth, death);
              pivot_column_index.insert(make_pair(pivot.getIndex(), i));
              break;
            }
          }
          else // if wc is empty, I output a PP as [birth,)
          {
            outputPP(-1, birth, dcg->threshold);
            break;
          }
        }
        else // (B) I have a new pivot and output PP as Writepairs
        {
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
    if (_birth != _death)
    {
      if (_death != dcg-> threshold)
      {
        wp->push_back(WritePairs2(_dim, _birth, _death));
      }
      else
      {
        wp->push_back(WritePairs2(-1, _birth, dcg -> threshold));
      }
    }
  }

  BirthdayIndex2 pop_pivot(priority_queue<BirthdayIndex2, vector<BirthdayIndex2>, BirthdayIndex2Comparator>& column)
  {
    if (column.empty())
    {
      return BirthdayIndex2(0, -1, 0);
    }
    else
    {
      auto pivot = column.top();
      column.pop();

      while (!column.empty() && column.top().index == pivot.getIndex())
      {
        column.pop();
        if (column.empty())
          return BirthdayIndex2(0, -1, 0);
        else
        {
          pivot = column.top();
          column.pop();
        }
      }
      return pivot;
    }
  }

  BirthdayIndex2 get_pivot(priority_queue<BirthdayIndex2, vector<BirthdayIndex2>, BirthdayIndex2Comparator>& column)
  {
    BirthdayIndex2 result = pop_pivot(column);
    if (result.getIndex() != -1)
    {
      column.push(result);
    }
    return result;
  }

  void assemble_columns_to_reduce()
  {
    ++dim;
    ctr->dim = dim;
    const int typenum = 2;
    if (dim == 1)
    {
      ctr->columns_to_reduce.clear();
      for (int y = 1; y <= ay; ++y)
      {
        for (int x = 1; x <= ax; ++x)
        {
          for (int m = 0; m < typenum; ++m)
          {
            double index = x | (y << 11) | (m << 21);
            if (pivot_column_index.find(index) == pivot_column_index.end())
            {
              double birthday = dcg -> getBirthday(index, 1);
              if (birthday != dcg -> threshold)
              {
                ctr -> columns_to_reduce.push_back(BirthdayIndex2(birthday, index, 1));
              }
            }
          }
        }
      }
    }
    sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndex2Comparator());
  }
};

// method = 0 --> link find algo (default)
// method = 1 --> compute pairs algo
// [[Rcpp::export]]
Rcpp::NumericMatrix cubical_2dim(const Rcpp::NumericMatrix& image, double threshold, int method)
{
  bool print = false;

  vector<WritePairs2> writepairs; // dim birth death
  writepairs.clear();

  DenseCubicalGrids2* dcg = new DenseCubicalGrids2(image, threshold);
  ColumnsToReduce2* ctr = new ColumnsToReduce2(dcg);

  switch(method)
  {
    case 0:
    {
      JointPairs2* jp = new JointPairs2(dcg, ctr, writepairs, print);
      jp->joint_pairs_main(); // dim0

      ComputePairs2* cp = new ComputePairs2(dcg, ctr, writepairs, print);
      cp->compute_pairs_main(); // dim1
      
      // free pointers
      delete jp;
      delete cp;

      break;
    }

    case 1:
    {
      ComputePairs2* cp = new ComputePairs2(dcg, ctr, writepairs, print);
      cp->compute_pairs_main(); // dim0
      cp->assemble_columns_to_reduce();

      cp->compute_pairs_main(); // dim1
      
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
