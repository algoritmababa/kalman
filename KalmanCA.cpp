// ============================================================================
//  kalmanca.cpp   (C++11)   --  2-durumlu CA Kalman filtresi (2 eksen: AZ, EL)
// ============================================================================
#include "KalmanCA.h"
#include <cmath>

namespace ptz {

KalmanCA::KalmanCA(double dt, double qJerk, double rMeas)
    : dt_(dt), q_(qJerk), r_(rMeas), init_(false),
      x0_{0.0, 0.0}, x1_{0.0, 0.0},
      p00_{0.0, 0.0}, p01_{0.0, 0.0}, p11_{0.0, 0.0} {}

void KalmanCA::predict() {
    if (!init_) return;

    const double dt = dt_;
    // Q : beyaz-jerk CA surec gurultusu
    const double dt2 = dt * dt, dt3 = dt2 * dt;

    for (int i = 0; i < 2; ++i) {          // AZ ve EL eksenleri
        // Durum: x = F x ,  F = [1 dt ; 0 1]
        x0_[i] += dt * x1_[i];             // w += a*dt
        // x1_[i] degismez

        // Kovaryans: P = F P F^T + Q
        const double n00 = p00_[i] + 2.0 * dt * p01_[i] + dt * dt * p11_[i];
        const double n01 = p01_[i] + dt * p11_[i];
        const double n11 = p11_[i];

        p00_[i] = n00 + q_ * dt3 / 3.0;
        p01_[i] = n01 + q_ * dt2 / 2.0;
        p11_[i] = n11 + q_ * dt;
    }
}

void KalmanCA::update(double zAz, double zEl) {
    const double z[2] = { zAz, zEl };

    if (!init_) {
        // Ilk olcumle baslat: w=z, a=0, ivmede buyuk belirsizlik.
        for (int i = 0; i < 2; ++i) {
            x0_[i]  = z[i];  x1_[i] = 0.0;
            p00_[i] = r_;    p01_[i] = 0.0;  p11_[i] = 1.0e3;
        }
        init_ = true;
        return;
    }

    for (int i = 0; i < 2; ++i) {          // AZ ve EL eksenleri
        // S = H P H^T + R = p00 + r   (H = [1 0])
        const double S  = p00_[i] + r_;
        const double K0 = p00_[i] / S;     // Kalman kazanci (w)
        const double K1 = p01_[i] / S;     // Kalman kazanci (a)

        const double y = z[i] - x0_[i];    // yenilik (innovation)
        x0_[i] += K0 * y;
        x1_[i] += K1 * y;

        // P = (I - K H) P   (H=[1 0]) -- bu yapida sonuc tam simetriktir
        const double p00 = p00_[i], p01 = p01_[i], p11 = p11_[i];
        p00_[i] = (1.0 - K0) * p00;
        p01_[i] = (1.0 - K0) * p01;
        p11_[i] = p11 - K1 * p01;
    }
}

double KalmanCA::predictAheadStd(int axis, double H) const {
    // var = [1 H] P [1 H]^T = p00 + 2 H p01 + H^2 p11
    const double var = p00_[axis] + 2.0 * H * p01_[axis] + H * H * p11_[axis];
    return var > 0.0 ? std::sqrt(var) : 0.0;
}

} // namespace ptz
