#include "snapshot_io.hpp"

#include <fstream>
#include <iomanip>
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
        if (b.N() != 1) {
            throw std::runtime_error("snapshots.csv writer currently expects N=1");
        }
        rows.push_back(SnapshotRow{
            this_event, event, step, t, phi, k,
            b.u.real(), b.u.imag(),
            b.A(0, 0).real(), b.A(0, 0).imag(),
            b.B(0, 0).real(), b.B(0, 0).imag(),
            b.R(0).real(), b.R(0).imag(),
            b.name});
        row_count++;
    }
}

void SnapshotSaver::write() const {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    f << "snapshot,event,step,t,phi,basis_index,u_re,u_im,"
         "A_re,A_im,B_re,B_im,R_re,R_im,name\n";
    f << std::setprecision(17);
    for (const auto& r : rows) {
        f << r.snapshot << ","
          << r.event << ","
          << r.step << ","
          << r.t << ","
          << r.phi << ","
          << r.basis_index << ","
          << r.u_re << "," << r.u_im << ","
          << r.A_re << "," << r.A_im << ","
          << r.B_re << "," << r.B_im << ","
          << r.R_re << "," << r.R_im << ","
          << r.name << "\n";
    }
}

}  // namespace ecg1d
