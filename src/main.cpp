#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <memory>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

typedef std::vector<double> Vector;
typedef std::vector<std::vector<double>> VecVector;


/*
 * Print Array in a form to copy-past in PyBlastAfterglow
 */
template<class T>
static void print_xy_as_numpy(T * x_arr, double * y_arr, int nxy, int ever_i=10) {
    int i = 0;
    std::cout << "x_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0) {
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << x_arr[i];
            else
                std::cout << x_arr[i];
        }
    }
    std::cout << "]) \n";

    std::cout << "y_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0){
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << y_arr[i];
            else
                std::cout << y_arr[i];
        }
    }
    std::cout <<  "]) \n";

    std::cout << "plt.loglog(x_arr, y_arr, ls='--', label='afg')" << "\n";
    std::cout << "\n";
}
template<class T>
static void print_xy_as_numpy(T & x_arr, T & y_arr, int nxy, int ever_i=10) {
    int i = 0;
    std::cout << "x_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0) {
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << x_arr[i];
            else if (i > 0)
                std::cout << ", " << x_arr[i];
            else
                std::cout << x_arr[i];
        }
    }
    std::cout << "]) \n";

    std::cout << "y_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0){
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << y_arr[i];
            else if (i > 0)
                std::cout << ", " << y_arr[i];
            else
                std::cout << y_arr[i];
        }
    }
    std::cout <<  "]) \n";

    std::cout << "plt.loglog(x_arr, y_arr, ls='--', label='afg')" << "\n";
    std::cout << "\n";
}
template<class T>
static void print_x_as_numpy(T & x_arr, int ever_i, std::string name="x_arr", std::string label="") {
    int i = 0;
    if (!label.empty()){
        std::cout << label << "\n";
    }
    std::cout << name << " = np.array([ ";
    for (i = 0; i < x_arr.size(); i++) {
        if (i % ever_i == 0) {
            if ( (i != x_arr.size() - 1) && i != 0 )
                std::cout << ", " << x_arr[i];
            else
                std::cout << x_arr[i];
        }
    }
    std::cout << "]) \n";
}


inline namespace CGS{
    const double &c  = 2.99792458e10;     // speed of light in cm / s
    const double &pi = 3.141592653589793; // pi
    const double &mp = 1.6726e-24;        // proton mass
    const double &me = 9.1094e-28;        //electron mass
    const double &mppme = mp + me;
    const double &qe = 4.803204e-10;
    const double &sigmaT = 6.6524e-25;
    const double &gamma_c_w_fac = 6 * pi * me * c / sigmaT; // Used in gamma_c calculation
    const double &kB = 1.380658e-16; // Boltzmann constant
    const double &hcgs = 6.6260755e-27; // Planck constant in cgs
    const double &h    = 6.6260755e-27; // erg s
    const double &lambda_c = (hcgs / (me * c)); // Compton wavelength
    const double &mec2 = 8.187105649650028e-07;  // erg # electron mass in erg, mass_energy equivalence
    const double &pc = 3.0857e18; // cm
    const double &year= 3.154e+7; // sec
    const double &day = 86400;
    const double &solar_m = 1.989e+33; // solar mass in g
    const double &deg2rad = 0.017453292;
    const double &cgs2mJy = 1.e26;
    const double &gravconst = 6.67259e-8;// # cm^3 g^-1 s^-2
    /// --------------------------------
    const double &SI_pc   =3.0856E+16; // m
    const double &RtD = 57.295779513082325;
    const double &A_RAD = 7.565767e-15; // radiation constant in unit of [erg/cm^3/K^4]
    const double &EV_TO_ERG = 1.602e-12; // [erg/eV]
    const double &M_ELEC = 9.11e-28; // electron mass in unit of [g]
    const double &ELEC = 4.803e-10; // electron charge in unit of [cgs]
    const double &SIGMA_T = 6.65e-25; // Thomson cross section in unit of [cm^2]
    const double &K_B = 1.38e-16; // Boltzman constant in unit of [erg/K]
    const double &H = 6.626e-27; // Planck constant in unit of [erg s]
    const double &M_PRO = 1.0/6.02e23; /* atomic mass in unit of [g] */
    const double &MeC2 = 8.1871398e-7; // in unit of [erg]
    const double &GAMMA13 = 2.67893; /* Gamma(1/3) */
};

namespace TOOLS{
    // --- make a vector of logarithmically spaced points
    template<typename T>
    class Logspace {
    private:
        T curValue, base;

    public:
        Logspace(T first, T base) : curValue(first), base(base) {}

        T operator()() {
            T retval = curValue;
            curValue *= base;
            return retval;
        }
    };

    static void MakeLogspaceVec(std::vector<double> & vec, const double &start, const double &stop, const double &base = 10) {
        int num = (int)vec.size();
        double realStart = std::pow(base, start);
        double realBase = std::pow(base, (stop-start)/num);
        std::generate_n(std::back_inserter(vec), num, Logspace<double>(realStart,realBase));
    }

    static std::vector<double> MakeLogspaceVec(const double &start,
                                               const double &stop,
                                               const int &num = 50,
                                               const double &base = 10) {
        double realStart = pow(base, start);
        double realBase = pow(base, (stop-start)/num);

        std::vector<double> retval;
        retval.reserve(num);

        std::generate_n(std::back_inserter(retval), num, Logspace<double>(realStart,realBase));
        return std::move(retval);

    }

    //
//    static std::valarray<double> MakeLogspace(const double &start,
//                                              const double &stop,
//                                              const int &num = 50,
//                                              const double &base = 10) {
//        auto retval = MakeLogspaceVec(start, stop, num, base);
//        return std::move(std::valarray<double> (retval.data(), retval.size()));
//
//    }

    std::vector<double> linspace(double first, double last, int len) {
        std::vector<double> result(len);
        double step = (last-first) / (len - 1);
        for (int i=0; i<len; i++) { result[i] = first + i*step; }
        return result;
    }

}

/// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
void sort_indexes(const std::vector<T> &v, std::vector<size_t> & idx) {

    // initialize original index locations
//    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

//    return std::move( idx );
}
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(begin(idx), end(idx), 0);

    // sort indexes based on comparing values in vv
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when vv contains elements of equal values
    std::stable_sort(begin(idx), end(idx),
                     [&v](size_t i1, size_t i2)
                     {return v[i1] < v[i2];});

    return std::move( idx );
}
Vector sort_by_indexes(const Vector & array, const std::vector<size_t> & indexes){
    if (array.size() != indexes.size()){
        std::cerr << AT << " size mismatch\n";
        exit(1);
    }
    Vector sorted (array.size());
    for(size_t i = 0; i < array.size(); i++){
        sorted[i] = array[indexes[i]];
    }
    return std::move( sorted );
}

static size_t findIndex( const double & x, const Vector & arr, size_t N ) {
    if(x <= arr[0])
        return 0;
    else if(x >= arr[N-1])
        return N-2;

    unsigned int i = ((unsigned int) N) >> 1;
    unsigned int a = 0;
    unsigned int b = N-1;

    // https://stackoverflow.com/questions/4192440/is-there-any-difference-between-1u-and-1-in-c/4192469
    // untill the m_size of b-a > 1 continue shrinking the array, approaching the 'x'
    while (b-a > 1u) // 1U is an unsigned value with the single bit 0 set
    {
        i = (b+a) >> 1; // ???
        if (arr[i] > x)
            b = i;
        else
            a = i;
    }

    return (int)a;
}

static inline double interpSegLin( size_t & a, size_t & b, const double & x, Vector & X, Vector & Y) {
    // take two indexis, 'a' 'b' and two arrays 'X' and 'Y' and interpolate
    // between them 'Y' for at 'x" of "X'
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];
    return ya + (yb-ya) * (x-xa)/(xb-xa);
}

static inline double interpSegLog( size_t & a, size_t & b, double x, Vector & X, Vector & Y) {
//        std::cout << a << ' ' << b << ' '<< x << ' '<< X << ' '<< Y << ' '<< N << "\n";
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];

    return ya * std::pow(yb/ya, log(x/xa)/log(xb/xa));
}

double Gamma(const double & beta){
//        return sqrt(1.0 / (1.0 - (beta * beta)));
    return sqrt(1. + (beta * beta / (1. - beta * beta)));
}

/// Blandford-McKee self-simialr solution
class BlandfordMcKee2{
    struct Pars{
        double Etot = -1;
        double k = 0;
        double Rscale = -1;
        double gamAdi = 4/3.;
        double eta   = 1.e-5;         //          eta = p/(rho*c^2)
        double n_ext = 1.e0;            // external medium number density

        double n0      = 1.;            // cm-3:    CBM number density
        double rho0    = n0*CGS::mp;    // g.cm-3:  comoving CBM mass density
        double rhoNorm = rho0;          // n0*mp_;       // g.cm-3:  comoving CBM mass density
        double lNorm = CGS::c;          // distance normalised to c
        double vNorm = CGS::c;          // velocity normalised to c
        double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

        double rho_ = 0;
        double Gamma_ = 0;
    };
    Pars * p_pars = nullptr;
public:
    BlandfordMcKee2(){
        p_pars = new Pars;
    }

    ~BlandfordMcKee2(){ delete p_pars; }
    void calcBM(double r, double t, double & rho, double & u, double & p){

        auto & _p = * p_pars;

        // returns normalised primitive variables for BM solution at position r (denormalised).

        // repeating from above (these are local)
        double lfacShock2 = _p.Etot * (17.- 4.*_p.k)/(8.*CGS::pi*_p.n_ext*CGS::mp*std::pow(t,3)*std::pow(CGS::c,5));
        // Blandford&McKee(1976) eq. 69
        double RShock = CGS::c * t * (1. - 1./(1. + 2. * (4. - _p.k) * lfacShock2));
        // Blandford&McKee(1976) eq. 27
        double rhoa = _p.n_ext * CGS::mp * std::pow(RShock/_p.Rscale, -_p.k);

        double pa = rhoa * _p.eta * CGS::c * CGS::c;
        double ha = CGS::c * CGS::c + pa * _p.gamAdi / (_p.gamAdi - 1.) / rhoa;
        double pf = 2./3. * lfacShock2 * ha * rhoa;
        double lfacf = sqrt(fmax(1.,.5 * lfacShock2));
        double Df = 2. * lfacShock2 * rhoa;

        // Blandford&McKee(1976) eq. 8-10
        double chi = (1. + 2. * (4. - _p.k) * lfacShock2) * ( 1. - r / ( CGS::c * t) );
        p = pf*std::pow(chi, -(17.-4.*_p.k)/(12.-3.*_p.k)) / _p.pNorm;
        double lfac = sqrt(lfacf * lfacf/chi + 1);  // (+1) to ensure lfac>1
        double D = Df * std::pow(chi, -(7.-2.*_p.k)/(4.-_p.k));
        // Blandford&McKee(1976) eq. 28-30 / 65-67
        rho = D / lfac / _p.rhoNorm;
        u = CGS::c * sqrt(1.-1./(lfac*lfac))*lfac / _p.vNorm;
    }

    void eval(double R, double Rshock, double Gamma, double n_ism, double rho2){
        auto & _p = * p_pars;
//        _p.n0 = n_ism ;
        _p.n_ext = n_ism ;

        double lfacShock2 = Gamma * Gamma;
//        double RShock = CGS::c * t * (1. - 1./(1. + 2. * (4. - _p.k) * lfacShock2));
        double t = Rshock / ( CGS::c * (1. - 1./(1. + 2. * (4. - _p.k) * lfacShock2)) );
        double rhoa = _p.n_ext * CGS::mp;// * pow(RShock/_p.Rscale, -_p.k);
        double Df = 2. * lfacShock2 * rhoa;
        double lfacf = sqrt(fmax(1.,.5 * lfacShock2));
        double chi = (1. + 2. * (4. - _p.k) * lfacShock2) * ( 1. - R / ( CGS::c * t) );
        double lfac = sqrt(lfacf * lfacf / chi + 1);  // (+1) to ensure lfac>1
        double D = Df * std::pow(chi, -(7. - 2. * _p.k)/(4. - _p.k));
        _p.rho_ = D / lfac / _p.rhoNorm;
        _p.Gamma_ = lfac;//CGS::c * sqrt(1.-1./(lfac*lfac)) / p_pars->vNorm;
        if(_p.rho_ > rho2 / CGS::mp){
            std::cerr << " rho behind shock(" << _p.rho_ << ") > rho2(" << rho2 << ") [BM]\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if(_p.Gamma_ > Gamma){
            std::cerr << " Gamma behind shock(" << _p.Gamma_ << ") > Gamma(" << Gamma << ") [BM]\n";
            std::cerr << AT << "\n";
            exit(1);
        }

    }

    double rho_downstream(double R, double Rshock, double Gamma, double n_ism, double rho2){
        eval(R,Rshock,Gamma,n_ism,rho2);
        return p_pars->rho_;
    }

    double gamma_downstream(double R, double Rshock, double Gamma, double n_ism, double rho2){
        eval(R,Rshock,Gamma,n_ism,rho2);
        return p_pars->Gamma_;
    }

};

/// Sedov-Taylor self-similar solution
class SedovTaylor{
    struct Pars{
        double gamma = -1;
        size_t nu = 0;
        double w = -1;
        // ------------
        double w1{}, w2{}, w3{};
        double b0{}, b1{}, b2{}, b3{}, b4{}, b5{}, b6{}, b7{}, b8{};
        double c0{}, c1{}, c2{}, c3{}, c4{}, c5{}, c6{}, c7{}, c8{};
        Vector f{};//, 1e5);
        Vector eta{};
        Vector d{};
        Vector p{};
        Vector vv{};
    };
    Pars * p_pars{};
public:
    SedovTaylor() {
        p_pars = new Pars();
        p_pars->f.resize(1e5);
        p_pars->eta.resize( p_pars->f.size() );
    }
    ~SedovTaylor() { delete p_pars; }
    void setPars(double gamma, size_t nu, double w){
        p_pars->gamma = gamma;
        p_pars->nu = nu;
        p_pars->w = w;
    }
    void evaluate(){

//        std::cout << AT << " \n Computing Sedov-Taylor profile for gamma="
//                  << p_pars->gamma<<" and grid of " << p_pars->eta.size() << " \n";

        auto nu = (double)p_pars->nu;
        double gamma = p_pars->gamma;
        double w = p_pars->w;
        // Constants for the parametric equations:
        p_pars->w1 = (3 * nu - 2 + gamma * (2 - nu)) / (gamma + 1.);
        p_pars->w2 = (2. * (gamma - 1) + nu) / gamma;
        p_pars->w3 = nu * (2. - gamma);

        p_pars->b0 = 1. / (nu * gamma - nu + 2);
        p_pars->b2 = (gamma - 1.) / (gamma * (p_pars->w2 - w));
        p_pars->b3 = (nu - w) / (float(gamma) * (p_pars->w2 - w));
        p_pars->b5 = (2. * nu - w * (gamma + 1)) / (p_pars->w3 - w);
        p_pars->b6 = 2. / (nu + 2 - w);
        p_pars->b1 = p_pars->b2 + (gamma + 1.) * p_pars->b0 - p_pars->b6;
        p_pars->b4 = p_pars->b1 * (nu - w) * (nu + 2. - w) / (p_pars->w3 - w);
        p_pars->b7 = w * p_pars->b6;
        p_pars->b8 = nu * p_pars->b6;

        // simple interpolation of correct function (only for nu=1,2,3)
        p_pars->c0 = 2 * (nu - 1) * CGS::pi + (nu - 2) * (nu - 3);
        p_pars->c5 = 2. / (gamma - 1);
        p_pars->c6 = (gamma + 1) / 2.;
        p_pars->c1 = p_pars->c5 * gamma;
        p_pars->c2 = p_pars->c6 / gamma;
        p_pars->c3 = (nu * gamma - nu + 2.) / ((p_pars->w1 - w) * p_pars->c6);
        p_pars->c4 = (nu + 2. - w) * p_pars->b0 * p_pars->c6;

        // Characterize the solution
        double f_min =  p_pars->w1 > w ? p_pars->c2 : p_pars->c6;

        p_pars->f = TOOLS::MakeLogspaceVec(log10(f_min), 0., (int)p_pars->f.size()); // log10(f_min)
        print_x_as_numpy(p_pars->f, 100, "f");


//                p_pars->c2 ? p_pars->w1 > w : p_pars->c6
        bool use_uniform_log_grid = false;
        if (use_uniform_log_grid){
            p_pars->f = TOOLS::MakeLogspaceVec(log10(f_min), 0., (int)p_pars->f.size());
        }
        else {
            double tmp_fmin = log10(f_min);
            double tmp_fmax = tmp_fmin;
            std::vector<size_t> segment_lengths = {10000, 1000, 100, 100};
            VecVector segments{};
            size_t tot_size = 0;
            for (size_t i_segment = 0; i_segment < segment_lengths.size() - 1; i_segment++) {
                double _p = (double) i_segment - (double) segment_lengths.size() - 1;
                double step = std::pow(10, (double) _p);
                if (tmp_fmax == 1) tmp_fmax = 0;
                //            std::cout << " tmp_min="<<tmp_fmin<<" tmp_max="<<tmp_fmin + step<<"\n";
                segments.emplace_back(TOOLS::MakeLogspaceVec(tmp_fmin, tmp_fmin + step, (int) segment_lengths[i_segment]));
                //            std::cout << segments[i_segment] << "\n";
                tmp_fmin = tmp_fmin + step;
                tot_size = tot_size + segment_lengths[i_segment];
            }
            segments.emplace_back(TOOLS::MakeLogspaceVec(tmp_fmin, 0, (int) segment_lengths[segment_lengths.size() - 1]));
//            std::cout << segments[segments.size() - 1] << "\n";
            tot_size = tot_size + segment_lengths[segment_lengths.size() - 1];
            //        std::cout << " tmp_min="<<tmp_fmin<<" tmp_max="<<0<<"\n";
            Vector tmp(tot_size);
            size_t ii = 0;
            for (size_t i_segment = 0; i_segment < segments.size(); i_segment++) {
                for (size_t ip = 0; ip < segment_lengths[i_segment]; ip++) {
                    tmp[ii] = segments[i_segment][ip];
                    ii++;
                }
            }
//            std::cout << tmp << '\n';
            p_pars->f = tmp; // log10(f_min)
        }


        // Sort the etas for our interpolation function
        p_pars->eta = parametrized_eta(p_pars->f);
        std::vector<size_t> idx = sort_indexes(p_pars->eta);
        p_pars->f = sort_by_indexes(p_pars->f, idx);
        p_pars->eta = sort_by_indexes(p_pars->eta, idx);

        p_pars->d = parametrized_d(p_pars->f);
        p_pars->p = parametrized_p(p_pars->f);
        p_pars->vv = parametrized_v(p_pars->f);

        print_x_as_numpy(p_pars->eta,100,"eta");
        print_x_as_numpy(p_pars->d,100,"d");

        if (p_pars->eta[0] > 0.){
            std::cerr << " in Sedov Taylow, the extend to eta is too short. Addd more 0.00000..." << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        // Finally Calculate the normalization of R_s:
//        auto tmp = std::pow(p_pars->vv, 2);
        Vector integral_(p_pars->eta.size());// = std::pow(p_pars->eta, nu - 1) * (p_pars->d * std::pow(p_pars->vv, 2.) + p_pars->p);
        for (size_t i = 0; i < p_pars->eta.size(); i++)
            integral_[i] = std::pow(p_pars->eta[i], nu - 1) * (p_pars->d[i] * std::pow(p_pars->vv[i], 2.) + p_pars->p[i]);

        Vector integral( integral_.size() - 1 );
        Vector deta (integral_.size() - 1);
        double integ = 0;
        for (size_t i = 0; i < integral.size(); i++){
            integral[i] = 0.5 * (integral_[i+1] + integral_[i]);
            deta[i] = (p_pars->eta[i+1] - p_pars->eta[i]);
            integ += integral[i] * deta[i];
        }
        double alpha = integ * (8 * p_pars->c0) / ((std::pow(gamma, 2) - 1.) * std::pow(nu + 2. - w, 2));
//        p_pars->_c = pow(1. / alpha, 1. / (nu + 2 - w))


    }

    Vector parametrized_eta(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < var.size(); i++)
            res[i] = std::pow(var[i], -p_pars->b6)
                     * std::pow((p_pars->c1 * (var[i] - p_pars->c2)), p_pars->b2)
                     * std::pow( (p_pars->c3 * (p_pars->c4 - var[i])), -p_pars->b1);
        return res;
//        return std::pow(var, -p_pars->b6) * std::pow((p_pars->c1 * (var - p_pars->c2)), p_pars->b2)
//               * std::pow( (p_pars->c3 * (p_pars->c4 - var)), -p_pars->b1);
    }

    Vector parametrized_d(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < var.size(); i++)
            res[i] = std::pow(var[i], -p_pars->b7)
                     * ( std::pow(p_pars->c1 * (var[i] - p_pars->c2), p_pars->b3 - p_pars->w * p_pars->b2) )
                     * ( std::pow(p_pars->c3 * (p_pars->c4 - var[i]), p_pars->b4 + p_pars->w * p_pars->b1) )
                     * std::pow(p_pars->c5 * (p_pars->c6 - var[i]),  -p_pars->b5 );
        return res;
//        return std::pow(var, -p_pars->b7)
//               * ( std::pow(p_pars->c1 * (var - p_pars->c2), p_pars->b3 - p_pars->w * p_pars->b2) )
//               * ( std::pow(p_pars->c3 * (p_pars->c4 - var), p_pars->b4 + p_pars->w * p_pars->b1) )
//               * std::pow(p_pars->c5 * (p_pars->c6 - var),  -p_pars->b5 );
    }

    Vector parametrized_p(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < var.size(); i++)
            res[i] = std::pow(var[i], p_pars->b8)
                     * ( std::pow(p_pars->c3 * (p_pars->c4 - var[i]), p_pars->b4 + (p_pars->w - 2) * p_pars->b1))
                     * ( std::pow(p_pars->c5 * (p_pars->c6 - var[i]), 1 - p_pars->b5));
        return res;
//        return std::pow(var, p_pars->b8)
//               * ( std::pow(p_pars->c3 * (p_pars->c4 - var), p_pars->b4 + (p_pars->w - 2) * p_pars->b1))
//               * ( std::pow(p_pars->c5 * (p_pars->c6 - var), 1 - p_pars->b5));
    }

    Vector parametrized_v(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < p_pars->eta.size(); i++)
            res[i] = parametrized_eta(var)[i] * var[i];
        return res;
//        return parametrized_eta(var) * var;
    }

    double rho_profile_int(double r, double r_shock, double rho2) {
        // density at radius r
        double eta = r / r_shock;
        size_t ia = findIndex(eta, p_pars->eta, p_pars->eta.size());
        size_t ib = ia + 1;
        if ((ia == 0) || (ib == p_pars->eta.size()-1)){
//            std::cerr << "Eta[0,1,2,3,4,5]="<<p_pars->eta[0]<<", "<<p_pars->eta[1]<<", "<<p_pars->eta[2]<<", "<<p_pars->eta[3]<<", "<<p_pars->eta[4]<<"\n";
//            std::cerr << "d[0,1,2,3,4,5]="<<p_pars->d[0]<<", "<<p_pars->d[1]<<", "<<p_pars->d[2]<<", "<<p_pars->d[3]<<", "<<p_pars->d[4]<<"\n";
//            std::cerr << AT << " ia=0 or ib=n-1 for eta="<<eta<<"\n";
//            exit(1);
        }
        double dint = interpSegLin(ia, ib, eta, p_pars->eta, p_pars->d);
        if (dint > 1){
            std::cerr << " Sedov profile: rho > rho2 \n";
            std::cerr << AT << "\n";
        }
        if (!std::isfinite(dint)){
            std::cerr << " NAN in ST interpolation" << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        return rho2 * dint;
    }

    double Gamma_profile_int(double r, double r_shock, double vshock){
        // velocity at radius r,
        double eta = r / r_shock;
        size_t ia = findIndex(eta, p_pars->eta, p_pars->eta.size());
        size_t ib = ia + 1;
        double vint = interpSegLin(ia, ib, eta, p_pars->eta, p_pars->vv);
        if (vshock * (2. / (p_pars->gamma + 1.)) > 1){
            std::cerr << " Sedov profile: vv > vshock \n";
            std::cerr << AT "\n";
        }
        return Gamma( vint * vshock * (2. / (p_pars->gamma + 1.)) );
    }

    double pressure_profile_int(double r, double r_shock, double p_shock){ // p' = -(gamAdi -1)Eint2/V2
        double eta = r / r_shock;
        size_t ia = findIndex(eta, p_pars->eta, p_pars->eta.size());
        size_t ib = ia + 1;
//        std::cout<<p_pars->p<<"\n";
        double pint = interpSegLin(ia, ib, eta, p_pars->eta, p_pars->p);
        return p_shock * pint;
    }
};

class SelfSimilar{
    std::unique_ptr<SedovTaylor> p_st{};
    std::unique_ptr<BlandfordMcKee2> p_bm{};
public:
    SelfSimilar(){
        p_st = std::make_unique<SedovTaylor>();
        p_bm = std::make_unique<BlandfordMcKee2>();

        /// initilize
        p_st->setPars(1.5, 3, 0.); // TODO this should not be here and evaluated for EVERY bw...
        p_st->evaluate();
    }
    double getRho(double r, double r_shock, double rho2, double GammaSh){
        if (GammaSh < 2){
            return p_st->rho_profile_int(r, r_shock, rho2);
        }
        else{
            return p_bm->rho_downstream(r, r_shock, GammaSh, 1, rho2);
        }
    }
};


int main() {
    FILE    *fout1  = nullptr; // Final model output file
    fout1 = fopen("../lh1.out", "w");
    if (fout1 == nullptr){
        throw std::runtime_error("file not found");
    }

    SelfSimilar ss = SelfSimilar();
    Vector Rsh     { 1e10, 1e12, 1e13, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20 };
    Vector GammaSh {  2.,   3.,   4.,   5.,   5.,   5.,   4.,   3.,   2.   };
    Vector rho2sh  { 1e-1, 4e-1, 8e-1, 1e0,  2e0,  3e0,  4e0,  5e0,  1e1 };
    Vector radii  = TOOLS::MakeLogspaceVec(10,20,100);
    size_t ii = 0;
    for (size_t i = 0; i < Rsh.size(); i++){
        Vector rho(radii.size());
        /// for each set of shock parameters compute the self-similar solution and save it into a file
        for (size_t j = 0; j < radii.size(); j++){
            rho[j] = ss.getRho(radii[j],Rsh[i],rho2sh[i],GammaSh[i]);
            fprintf(fout1, "%d %lf %d %lf %lf %lf %lf\n",
                    (int)i,
                    radii[j],
                    (int)j,
                    Rsh[i],
                    GammaSh[i],
                    rho2sh[i],
                    rho[j]
            );
            ii++;
        }
    }
    fclose(fout1);
    fout1 = nullptr;
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
