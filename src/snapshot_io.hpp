#pragma once

#include "basis_params.hpp"

#include <string>
#include <vector>

namespace ecg1d {

// One buffered snapshots.csv row. The current writer is for N=1 runs:
// A/B are 1x1 and R is scalar.
struct SnapshotRow {
    int snapshot;
    std::string event;
    int step;
    double t;
    double phi;
    int basis_index;
    double u_re, u_im;
    double A_re, A_im;
    double B_re, B_im;
    double R_re, R_im;
    int name;
};

// Buffer snapshots in memory and write snapshots.csv once at the end.
struct SnapshotSaver {
    std::vector<SnapshotRow> rows;
    int event_idx = 0;
    int row_count = 0;
    std::string path;

    explicit SnapshotSaver(const std::string& task_dir);

    void save(const std::string& event,
              int step,
              double t,
              double phi,
              const std::vector<BasisParams>& basis);

    void write() const;
};

}  // namespace ecg1d
