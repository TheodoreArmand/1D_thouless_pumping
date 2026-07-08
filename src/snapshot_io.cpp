#include "snapshot_io.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ecg1d {

SnapshotSaver::SnapshotSaver(const std::string& task_dir)
    : path(task_dir + "/snapshots.csv") {}

void SnapshotSaver::save(const std::string& event,
                         int step,
                         double t,
                         double phi,
                         const std::vector<BasisParams>& basis) {
    const int this_event = event_idx++;
    for (int k = 0; k < static_cast<int>(basis.size()); ++k) {
        const auto& b = basis[k];
        const int n = b.N();
        if (N < 0) {
            N = n;
        } else if (N != n) {
            throw std::runtime_error("snapshots.csv writer: mixed N across basis functions");
        }
        std::vector<double> vals;
        vals.reserve(2 * (1 + n * (n + 1) / 2 + 2 * n));
        vals.push_back(b.u.real());
        vals.push_back(b.u.imag());
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                vals.push_back(b.A(i, j).real());
                vals.push_back(b.A(i, j).imag());
            }
        }
        for (int a = 0; a < n; ++a) {
            for (int c = 0; c < n; ++c) {
                if (c != a && std::abs(b.B(a, c)) != 0.0) {
                    throw std::runtime_error(
                        "snapshots.csv writer: nonzero off-diagonal B (engine convention is diagonal B)");
                }
            }
            vals.push_back(b.B(a, a).real());
            vals.push_back(b.B(a, a).imag());
        }
        for (int a = 0; a < n; ++a) {
            vals.push_back(b.R(a).real());
            vals.push_back(b.R(a).imag());
        }
        rows.push_back(SnapshotRow{this_event, event, step, t, phi, k, std::move(vals), b.name});
        row_count++;
    }
}

void SnapshotSaver::write() const {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    if (N <= 1) {
        // Legacy N=1 header, byte-identical to the original writer
        // (also used when no rows were ever saved).
        f << "snapshot,event,step,t,phi,basis_index,u_re,u_im,"
             "A_re,A_im,B_re,B_im,R_re,R_im,name\n";
    } else {
        std::ostringstream h;
        h << "snapshot,event,step,t,phi,basis_index,u_re,u_im";
        for (int i = 0; i < N; ++i)
            for (int j = i; j < N; ++j)
                h << ",A" << i << j << "_re,A" << i << j << "_im";
        for (int a = 0; a < N; ++a) h << ",B" << a << a << "_re,B" << a << a << "_im";
        for (int a = 0; a < N; ++a) h << ",R" << a << "_re,R" << a << "_im";
        h << ",name\n";
        f << h.str();
    }
    f << std::setprecision(17);
    for (const auto& r : rows) {
        f << r.snapshot << ","
          << r.event << ","
          << r.step << ","
          << r.t << ","
          << r.phi << ","
          << r.basis_index;
        for (const double v : r.vals) f << "," << v;
        f << "," << r.name << "\n";
    }
}

}  // namespace ecg1d
