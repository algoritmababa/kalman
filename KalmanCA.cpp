// ============================================================================
//  kalmanca.cpp   (C++11)   --  2-durumlu CA Kalman filtresi (2 eksen: AZ, EL)
//                               DEGISKEN dt: her adimda gecen sure olculur.
// ============================================================================
#include "KalmanCA.h"
#include <cmath>

namespace ptz {

KalmanCA::KalmanCA(double qJerk, double rMeas)
    : q_(qJerk), r_(rMeas), init_(false),
      x0_{0.0, 0.0}, x1_{0.0, 0.0},
      p00_{0.0, 0.0}, p01_{0.0, 0.0}, p11_{0.0, 0.0},
      lastUpdateSec_(0.0), lastDt_(0.0),
      lastHorizon_{0.0, 0.0} {}

// Durumu dt kadar ileri tasi:  x = F x ,  P = F P F^T + Q
void KalmanCA::predictStep(double dt) {
    const double dt2 = dt * dt, dt3 = dt2 * dt;

    for (int i = 0; i < 2; ++i) {          // AZ ve EL eksenleri
        // Durum: F = [1 dt ; 0 1]
        x0_[i] += dt * x1_[i];             // w += a*dt
        // x1_[i] degismez

        // Kovaryans: P = F P F^T + Q  (Q : beyaz-jerk CA surec gurultusu)
        const double n00 = p00_[i] + 2.0 * dt * p01_[i] + dt2 * p11_[i];
        const double n01 = p01_[i] + dt * p11_[i];
        const double n11 = p11_[i];

        p00_[i] = n00 + q_ * dt3 / 3.0;
        p01_[i] = n01 + q_ * dt2 / 2.0;
        p11_[i] = n11 + q_ * dt;
    }
    lastDt_ = dt;
}

void KalmanCA::update(double zAz, double zEl) {
    updateAt(zAz, zEl, nowSec());
}

void KalmanCA::update(double zAz, double zEl, double tSec) {
    updateAt(zAz, zEl, tSec);
}

void KalmanCA::updateAt(double zAz, double zEl, double tSec) {
    const double z[2] = { zAz, zEl };

    if (!init_) {
        // Ilk olcumle baslat: w=z, a=0, ivmede buyuk belirsizlik.
        for (int i = 0; i < 2; ++i) {
            x0_[i]  = z[i];  x1_[i] = 0.0;
            p00_[i] = r_;    p01_[i] = 0.0;  p11_[i] = 1.0e3;
        }
        lastUpdateSec_ = tSec;
        lastDt_ = 0.0;
        init_ = true;
        return;
    }

    // dt = son update'ten bu yana gecen sure (filtre icinde hesaplanir).
    double dt = tSec - lastUpdateSec_;
    lastUpdateSec_ = tSec;
    if (dt < 0.0) dt = 0.0;                // monoton olmayan zaman damgasina karsi koruma
    predictStep(dt);

    // Olcum duzeltmesi
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

// Ufuk = SON OLCUM anindan (lastUpdateSec_) bu yana gecen sure.
// Sensor kesildiginde bu sure buyur -> cikti zamana bagli olu-hesap yapar.
double KalmanCA::predictAheadAz(double tSec) {
    if (!init_) { lastHorizon_[AZ] = 0.0; return 0.0; }
    double H = tSec - lastUpdateSec_;
    if (H < 0.0) H = 0.0;
    lastHorizon_[AZ] = H;
    return x0_[AZ] + x1_[AZ] * H;          // w(now) = w + a*H
}

double KalmanCA::predictAheadEl(double tSec) {
    if (!init_) { lastHorizon_[EL] = 0.0; return 0.0; }
    double H = tSec - lastUpdateSec_;
    if (H < 0.0) H = 0.0;
    lastHorizon_[EL] = H;
    return x0_[EL] + x1_[EL] * H;          // w(now) = w + a*H
}

double KalmanCA::predictAheadStd(int axis, double H) const {
    // P'nin H kadar ileri tasinmis hali:  P(H) = F(H) P F(H)^T + Q(H)
    // var(w) = p00 + 2 H p01 + H^2 p11 + q*H^3/3
    // Son terim (Q birikimi) kesinti uzadikca belirsizligin dogru sekilde
    // buyumesini saglar.
    const double var = p00_[axis] + 2.0 * H * p01_[axis] + H * H * p11_[axis]
                     + q_ * H * H * H / 3.0;
    return var > 0.0 ? std::sqrt(var) : 0.0;
}

} // namespace ptz
