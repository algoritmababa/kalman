// ============================================================================
//  kalmanca.h   (C++11)
//
//  Acisal hiz tahmini icin Sabit-Ivme (Constant Acceleration, CA) Kalman
//  filtresi -- 2 BOYUTLU: azimuth ve elevation eksenleri.
//  Harici bagimlilik yok (Eigen vb. gerekmez); 2x2 islemler elle acilmistir
//  -> embedded uygun.
//
//  Iki eksen dinamik olarak bagimsizdir (capraz kuplaj yok); her eksende
//  ayni CA modeli, ayni F/Q/R semasi ve ayni parametreler kullanilir.
//
//  Durum :  x = [ w , a ]^T     w = acisal hiz (deg/s), a = acisal ivme (deg/s^2)
//           her eksen (az, el) icin ayri
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
    // Eksen indeksleri (iceride dizi erisimi ve isteğe bagli genel API icin)
    enum Axis { AZ = 0, EL = 1 };

    // dt     : adim periyodu (s)
    // qJerk  : surec gurultusu (jerk) spektral yogunlugu  [ (deg/s^2)^2 / s ]
    // rMeas  : omega olcum gurultu varyansi               [ (deg/s)^2 ]
    KalmanCA(double dt, double qJerk, double rMeas);

    void predict();                         // her iki eksen bir adim ileri (dt)
    void update(double zAz, double zEl);    // olcum guncellemesi (z = w);
                                            // ilk cagri filtreyi baslatir

    // Tahmini acisal hiz (deg/s)
    double omegaAz() const { return x0_[AZ]; }
    double omegaEl() const { return x0_[EL]; }
    // Tahmini acisal ivme (deg/s^2)
    double alphaAz() const { return x1_[AZ]; }
    double alphaEl() const { return x1_[EL]; }

    // H saniye ilerisi icin acisal hiz tahmini ve belirsizligi.
    double predictAheadAz(double horizonSec) const { return x0_[AZ] + x1_[AZ] * horizonSec; }
    double predictAheadEl(double horizonSec) const { return x0_[EL] + x1_[EL] * horizonSec; }
    double predictAheadStdAz(double horizonSec) const { return predictAheadStd(AZ, horizonSec); }
    double predictAheadStdEl(double horizonSec) const { return predictAheadStd(EL, horizonSec); }

    bool   isInitialized() const { return init_; }

    void reset() {
        init_ = false;
        for (int i = 0; i < 2; ++i) {
            x0_[i] = x1_[i] = 0.0;
            p00_[i] = p01_[i] = p11_[i] = 0.0;
        }
    }

    double getx0Az() const { return x0_[AZ]; }
    double getx0El() const { return x0_[EL]; }

private:
    double predictAheadStd(int axis, double H) const;   // 1-sigma (deg/s)

    double dt_, q_, r_;
    bool   init_;
    double x0_[2], x1_[2];            // [AZ, EL] :  w, a
    double p00_[2], p01_[2], p11_[2]; // [AZ, EL] :  simetrik kovaryans
};

} // namespace ptz

#endif // KALMANCA_H
