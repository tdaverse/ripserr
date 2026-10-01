#include <queue>
#include <Rcpp.h>

#include "cube.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
#include "config.h"
#include "dense_cubical_grids.h"
#include "ph_2d.h"

using namespace std;

// [[Rcpp::export]]
Rcpp::NumericMatrix cubical_compute(Rcpp::NumericVector data, 
    Rcpp::IntegerVector dims, 
    int maxdim = 3,
    bool top_dim = false,
    bool embedded = false, 
    double threshold = 9999 ){

    // modified logic from cubicalripser_pybind.h

    // cubicalripserq: no need for representatives part

    Config config;
    // if (threshold == -1){
    //     threshold = DBL_MAX; //compileAttributes doesn't know what DBL_MAX is?
    // }
    config.threshold = threshold;

    // cubicalripserq: format only needed for CLI + deleting any representatives references to avoid confusing myself
    //config.format = NUMPY;

    vector<WritePairs> writepairs; // (dim birth death x y z)
    writepairs.reserve(1000);

    std::unique_ptr<DenseCubicalGrids> dcg;
    vector<Cube> ctr;

    const size_t nd = dims.size();
    if (nd < 1 || nd > 4) {
        throw std::invalid_argument("computePH: input array must have 1 to 4 dimensions");
    }

    // cubicalripserq: R arrays are already Fortran order; keep variable as it's used in line 60
    bool fortran_order = true;


    const uint8_t ndim = static_cast<uint8_t>(nd);
    config.maxdim = maxdim;
    const uint32_t sx = static_cast<uint32_t>(dims[0]);
    const uint32_t sy = (nd > 1) ? static_cast<uint32_t>(dims[1]) : 1u;
    const uint32_t sz = (nd > 2) ? static_cast<uint32_t>(dims[2]) : 1u;
    const uint32_t sw = (nd > 3) ? static_cast<uint32_t>(dims[3]) : 1u;
    dcg = std::make_unique<DenseCubicalGrids>(config, ndim, sx, sy, sz, sw);
    config.maxdim = std::min<uint8_t>(config.maxdim, dcg->dim - 1);
    if (top_dim && dcg->dim > 1) {
        config.method = ALEXANDER;
        config.embedded = !embedded;
    } else {
        config.embedded = embedded;
    }

    dcg->gridFromArray(data.begin(), embedded, fortran_order);
    dcg->finalisePadding();

    // compute PH
    if (config.method == ALEXANDER) {
        auto jp = std::make_unique<JointPairs>(dcg.get(), writepairs, config);
        if (dcg->dim == 1) {
            jp->enum_edges({0}, ctr);
            jp->joint_pairs_main(ctr, 0); // dim0
        } else if (dcg->dim == 2) {
            jp->enum_edges({0, 1, 3, 4}, ctr);
            jp->joint_pairs_main(ctr, 1); // dim1
        } else if (dcg->dim == 3) {
            jp->enum_edges({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, ctr);
            jp->joint_pairs_main(ctr, 2); // dim2
        }
    } else {

        bool fastpath_handled = false;
        if (dcg->dim <= 2 && dcg->az == 1 && dcg->aw == 1) {
            fastpath_handled = compute_PH_2d(dcg.get(), writepairs, config);
        }
        if (!fastpath_handled) {
            auto jp = std::make_unique<JointPairs>(dcg.get(), writepairs, config);
            std::vector<uint32_t> betti;
            if (dcg->dim == 1) {
                jp->enum_edges({0}, ctr);
            } else if (dcg->dim == 2) {
                jp->enum_edges({0, 1}, ctr);
            } else if (dcg->dim == 3) {
                jp->enum_edges({0, 1, 2}, ctr);
            } else if (dcg->dim == 4) {
                jp->enum_edges({0, 1, 2, 3}, ctr);
            }
            jp->joint_pairs_main(ctr, 0); // dim0
            betti.push_back(writepairs.size());
            if (config.maxdim > 0) {
                ComputePairs cp(dcg.get(), writepairs, config);
                cp.compute_pairs_main(ctr); // dim1
                betti.push_back(writepairs.size() - betti[0]);
                if (config.maxdim > 1) {
                    cp.assemble_columns_to_reduce(ctr, 2);
                    cp.compute_pairs_main(ctr); // dim2
                    betti.push_back(writepairs.size() - betti[0] - betti[1]);
                    if (config.maxdim > 2) {
                        cp.assemble_columns_to_reduce(ctr, 3);
                        cp.compute_pairs_main(ctr); // dim3
                        betti.push_back(writepairs.size() - betti[0] - betti[1] - betti[2]);
                    }
                }
            }
        }
    }

    // result
    // determine shift between dcg and the voxel coordinates
    auto pad_x = (dcg->ax - dcg->img_x) / 2;
    auto pad_y = (dcg->ay - dcg->img_y) / 2;
    auto pad_z = (dcg->az - dcg->img_z) / 2;
    auto pad_w = (dcg->aw - dcg->img_w) / 2;
    const int64_t p = static_cast<int64_t>(writepairs.size());
    const int num_column = (dcg->dim > 3) ? 11 : 9;


    // Redo result construction for Rcpp instead of nanobind
    Rcpp::NumericMatrix result(p, num_column);
    for (int64_t i = 0; i < p; i++){
        result(i, 0) = writepairs[i].dim;
        result(i, 1) = writepairs[i].birth;
        result(i, 2) = writepairs[i].death;
        result(i, 3) = writepairs[i].birth_x - pad_x;
        result(i, 4) = writepairs[i].birth_y - pad_y;
        result(i, 5) = writepairs[i].birth_z - pad_z;
        if (dcg->dim > 3){
            result(i, 6) = writepairs[i].birth_w - pad_w;
            result(i, 7) = writepairs[i].death_x - pad_x;
            result(i, 8) = writepairs[i].death_y - pad_y;
            result(i, 9) = writepairs[i].death_z - pad_z;
            result(i, 10) = writepairs[i].death_w - pad_w;
        } else {
            result(i, 6) = writepairs[i].death_x - pad_x;
            result(i, 7) = writepairs[i].death_y - pad_y;
            result(i, 8) = writepairs[i].death_z - pad_z;
        }
    }
    
    return result;

}