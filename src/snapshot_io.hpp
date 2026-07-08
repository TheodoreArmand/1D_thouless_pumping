#pragma once

#include "basis_params.hpp"

#include <string>
#include <vector>

namespace ecg1d {

// One buffered snapshots.csv row, N-general.
// vals holds, in order: u_re,u_im; A(i,j) for i<=j (upper triangle, re/im);
// B(a,a) diagonal (re/im); R(a) (re/im). For N=1 this is exactly the legacy
// sequence u,A,B,R, so N=1 output stays byte-identical to the old writer.
struct SnapshotRow {
    int snapshot;
    std::string event;
    int step;
    double t;
    double phi;
    int basis_index;
    std::vector<double> vals;
    int name;
};

// Buffer snapshots in memory and write snapshots.csv once at the end.
// The particle number N is captured from the first save() and fixes the
// header layout (legacy columns for N=1; A00,A01,A11,B00,B11,R0,R1,... for
// N>=2). B must be diagonal (engine convention); nonzero off-diagonals throw.
struct SnapshotSaver {
    std::vector<SnapshotRow> rows;
    int event_idx = 0;
    int row_count = 0;
    int N = -1;
    std::string path;

    explicit SnapshotSaver(const std::string& task_dir);

    void save(const std::string& event,
              int step,
              double t,
              double phi,
              const std::vector<BasisParams>& basis);

    // Write the buffered rows to snapshots.csv (single flush).
    void write() const;
};

}  // namespace ecg1d
