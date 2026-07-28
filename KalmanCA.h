// ============================================================================
//  kalmanca.h   (C++11)
//
//  Acisal hiz tahmini icin Sabit-Ivme (Constant Acceleration, CA) Kalman
//  filtresi -- 2 BOYUTLU: azimuth ve elevation eksenleri.
//  2x2 islemler elle acilmistir.  Tek harici baglilik: <chrono> (dahili
//  zaman olcumu icin standart kutuphane).
//
//  DEGISKEN dt: adim periyodu sabit degildir; her cagride gecen sure
//  filtre icinde saatten (steady_clock) olculur. Isterseniz olcumun kendi
//  zaman damgasini (tSec) da verebilirsiniz (degisken-dt icin onerilir).
//
//  Iki eksen dinamik olarak bagimsizdir (capraz kuplaj yok); her eksende
//  ayni CA modeli, ayni F/Q/R semasi ve ayni parametreler kullanilir.
//
//  Durum :  x = [ w , a ]^T     w = acisal hiz (deg/s), a = acisal ivme (deg/s^2)
//           her eksen (az, el) icin ayri
//  Model :  w' = a ,  a' = beyaz gurultu (jerk)         (theta uzerinde CA)
//           F = [ 1  dt ; 0  1 ]   (dt her adimda olculur)
//  Olcum :  z = w                                       H = [ 1  0 ]
//  Cikti :  w(t+H) = w + a*H                            (predictAhead)
// ============================================================================
#ifndef KALMANCA_H
#define KALMANCA_H

#include <chrono>

namespace ptz {

class KalmanCA {
public:
    // Eksen indeksleri (iceride dizi erisimi ve isteğe bagli genel API icin)
    enum Axis { AZ = 0, EL = 1 };

    // qJerk  : surec gurultusu (jerk) spektral yogunlugu  [ (deg/s^2)^2 / s ]
    // rMeas  : omega olcum gurultu varyansi               [ (deg/s)^2 ]
    // dt artik parametre DEGIL: her adimda gecen sure filtre icinde olculur.
    KalmanCA(double qJerk, double rMeas);

    // Olcum guncellemesi (z = w). dt = son update'ten bu yana gecen sure,
    // saatten (steady_clock) olculur. Ilk cagri filtreyi baslatir (predict yok).
    void update(double zAz, double zEl);

    // Ayni guncelleme, ancak dt'yi cagiran verir: tSec = olcumun monoton
    // artan zaman damgasi (saniye). Degisken-dt sistemlerde onerilen yol
    // (sensor zaman damgasi, saat okumasindan daha isabetlidir).
    void update(double zAz, double zEl, double tSec);

    // Tahmini acisal hiz (deg/s)
    double omegaAz() const { return x0_[AZ]; }
    double omegaEl() const { return x0_[EL]; }
    // Tahmini acisal ivme (deg/s^2)
    double alphaAz() const { return x1_[AZ]; }
    double alphaEl() const { return x1_[EL]; }

    // Ileri tahmin (olu-hesap / dead reckoning). Ufuk (H) SON OLCUM anindan
    // (son update) bu yana gecen sure olarak filtre icinde olculur:
    //     H = simdi - lastUpdateSec_ ,   w(now) = w + a*H
    // Sensor kesildiginde H zamanla buyur -> cikti zamana bagli ilerler.
    // Olcum her adimda geliyorsa H~0 -> mevcut hiz. Kendi zaman damgasini
    // (lastHorizon_) yazar -> const degil.
    double predictAheadAz() { return predictAheadAz(nowSec()); }
    double predictAheadEl() { return predictAheadEl(nowSec()); }
    // Kendi zaman damganizla (test / harici saat / sensor zaman damgasi):
    double predictAheadAz(double tSec);
    double predictAheadEl(double tSec);

    // Son predictAhead ufkundaki 1-sigma belirsizlik (deg/s). Kesinti boyunca
    // surec gurultusu (Q) birikerek buyur. Once ilgili predictAhead* cagrilir.
    double predictAheadStdAz() const { return predictAheadStd(AZ, lastHorizon_[AZ]); }
    double predictAheadStdEl() const { return predictAheadStd(EL, lastHorizon_[EL]); }

    // En son olcumden bu yana gecen sure (s). Kesinti (bayatlik) olcusu.
    double timeSinceUpdate() const { return init_ ? (nowSec() - lastUpdateSec_) : 0.0; }

    double lastDt() const { return lastDt_; }   // son update adiminin dt'si (s)

    bool   isInitialized() const { return init_; }

    void reset() {
        init_ = false;
        for (int i = 0; i < 2; ++i) {
            x0_[i] = x1_[i] = 0.0;
            p00_[i] = p01_[i] = p11_[i] = 0.0;
            lastHorizon_[i]  = 0.0;
        }
        lastUpdateSec_ = 0.0;
        lastDt_ = 0.0;
    }

    double getx0Az() const { return x0_[AZ]; }
    double getx0El() const { return x0_[EL]; }

private:
    typedef std::chrono::steady_clock Clock;
    static double nowSec() {
        return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();
    }

    void   updateAt(double zAz, double zEl, double tSec);
    void   predictStep(double dt);                      // durumu dt kadar ilerlet
    double predictAheadStd(int axis, double H) const;   // 1-sigma (deg/s)

    double q_, r_;
    bool   init_;
    double x0_[2], x1_[2];            // [AZ, EL] :  w, a
    double p00_[2], p01_[2], p11_[2]; // [AZ, EL] :  simetrik kovaryans

    double lastUpdateSec_;            // son update zaman damgasi (s)
    double lastDt_;                   // son update adiminin dt'si (s)
    double lastHorizon_[2];           // eksen basina son kullanilan ufuk (std icin)
};

} // namespace ptz

#endif // KALMANCA_H
