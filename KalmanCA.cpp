// ============================================================================
//  kalmanca.cpp   (C++11)   --  2-durumlu CA Kalman filtresi
// ============================================================================
#include "KalmanCA.h"
#include <cmath>

namespace ptz {

KalmanCA::KalmanCA(double dt, double qJerk, double rMeas)
    : dt_(dt), q_(qJerk), r_(rMeas), init_(false),
      x0_(0.0), x1_(0.0), p00_(0.0), p01_(0.0), p11_(0.0) {}

void KalmanCA::predict() {
    if (!init_) return;

    const double dt = dt_;
    // Durum: x = F x ,  F = [1 dt ; 0 1]
    x0_ += dt * x1_;            // w += a*dt
    // x1_ degismez

    // Kovaryans: P = F P F^T + Q
    const double n00 = p00_ + 2.0 * dt * p01_ + dt * dt * p11_;
    const double n01 = p01_ + dt * p11_;
    const double n11 = p11_;

    // Q : beyaz-jerk CA surec gurultusu
    const double dt2 = dt * dt, dt3 = dt2 * dt;
    p00_ = n00 + q_ * dt3 / 3.0;
    p01_ = n01 + q_ * dt2 / 2.0;
    p11_ = n11 + q_ * dt;
}

void KalmanCA::update(double z) {
    if (!init_) {
        // Ilk olcumle baslat: w=z, a=0, ivmede buyuk belirsizlik.
        x0_  = z;   x1_ = 0.0;
        p00_ = r_;  p01_ = 0.0;  p11_ = 1.0e3;
        init_ = true;
        return;
    }

    // S = H P H^T + R = p00 + r   (H = [1 0])
    const double S  = p00_ + r_;
    const double K0 = p00_ / S;        // Kalman kazanci (w)
    const double K1 = p01_ / S;        // Kalman kazanci (a)

    const double y = z - x0_;          // yenilik (innovation)
    x0_ += K0 * y;
    x1_ += K1 * y;

    // P = (I - K H) P   (H=[1 0]) -- bu yapida sonuc tam simetriktir
    const double p00 = p00_, p01 = p01_, p11 = p11_;
    p00_ = (1.0 - K0) * p00;
    p01_ = (1.0 - K0) * p01;
    p11_ = p11 - K1 * p01;
}

double KalmanCA::predictAhead(double H) const {
    return x0_ + x1_ * H;              // w(t+H) = w + a*H
}

double KalmanCA::predictAheadStd(double H) const {
    // var = [1 H] P [1 H]^T = p00 + 2 H p01 + H^2 p11
    const double var = p00_ + 2.0 * H * p01_ + H * H * p11_;
    return var > 0.0 ? std::sqrt(var) : 0.0;
}

} // namespace ptz
