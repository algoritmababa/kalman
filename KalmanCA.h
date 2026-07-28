// ============================================================================
//  kalmanca.h   (C++11)
//
//  Acisal hiz tahmini icin Sabit-Ivme (Constant Acceleration, CA) Kalman
//  filtresi. Harici bagimlilik yok (Eigen vb. gerekmez); 2x2 islemler elle
//  acilmistir -> embedded uygun.
//
//  Durum :  x = [ w , a ]^T     w = acisal hiz (deg/s), a = acisal ivme (deg/s^2)
//  Model :  w' = a ,  a' = beyaz gurultu (jerk)         (theta uzerinde CA)
//           F = [ 1  dt ; 0  1 ]
//  Olcum :  z = w                                       H = [ 1  0 ]
//  Cikti :  w(t+H) = w + a*H                            (predictAhead)
// ============================================================================
#ifndef KALMANCA_H
#define KALMANCA_H

namespace ptz {

class KalmanCA {
public:
    // dt     : adim periyodu (s)
    // qJerk  : surec gurultusu (jerk) spektral yogunlugu  [ (deg/s^2)^2 / s ]
    // rMeas  : omega olcum gurultu varyansi               [ (deg/s)^2 ]
    KalmanCA(double dt, double qJerk, double rMeas);

    void predict();             // bir adim ileri (dt)
    void update(double zOmega); // olcum guncellemesi (z = w); ilk cagri filtreyi baslatir

    double omega() const { return x0_; }   // tahmini acisal hiz   (deg/s)
    double alpha() const { return x1_; }    // tahmini acisal ivme  (deg/s^2)

    // H saniye ilerisi icin acisal hiz tahmini ve belirsizligi.
    double predictAhead(double horizonSec) const;     // w(t+H) = w + a*H
    double predictAheadStd(double horizonSec) const;   // 1-sigma (deg/s)

    bool   isInitialized() const { return init_; }

    void reset() { init_ = false; x0_ = x1_ = 0.0; p00_ = p01_ = p11_ = 0.0; }
    double getx0() { return x0_;    }

private:
    double dt_, q_, r_;
    bool   init_;
    double x0_, x1_;            // w, a
    double p00_, p01_, p11_;    // simetrik kovaryans
};

} // namespace ptz

#endif // KALMANCA_H
