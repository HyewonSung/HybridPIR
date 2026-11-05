#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>
#include <sys/time.h>
#include "generic_utils.h"
#include "spqlios/lagrangehalfc_impl.h"
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstring>
#include "64header.h"
#include "poc_64types.h"

Random* global_random = new Random();

static inline std::string pretty_bytes(uint64_t bytes) {
    const double KiB = 1024.0;
    const double MiB = KiB * 1024.0;
    const double GiB = MiB * 1024.0;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    if (bytes >= (uint64_t)GiB) {
        oss << (bytes / GiB) << " GiB";
    } else if (bytes >= (uint64_t)MiB) {
        oss << (bytes / MiB) << " MiB";
    } else if (bytes >= (uint64_t)KiB) {
        oss << (bytes / KiB) << " KiB";
    } else {
        oss << bytes << " B";
    }
    return oss.str();
}

static inline void torusPolyClear(Torus64Polynomial* p, int N) {
    std::memset(p->coefs, 0, sizeof(Torus64) * N);
}
static inline void torusPolySubTo(Torus64Polynomial* dst, const Torus64Polynomial* src, int N) {
    for (int i = 0; i < N; i++) dst->coefs[i] -= src->coefs[i];
}

static inline void lagrangeCopy(LagrangeHalfCPolynomiala* dst, const LagrangeHalfCPolynomiala* src, int N) {
    std::memcpy(dst->values, src->values, sizeof(double) * (size_t)N);
}

static inline void tlweFFTClear(TLweSampleFFTa* ct, int N, const Globals* env) {
    (void)N;

    LagrangeHalfCPolynomialClear_lvl2(&ct->a[0], env);
    LagrangeHalfCPolynomialClear_lvl2(ct->b, env);
}

static inline void tlweFFTCopy(TLweSampleFFTa* dst, const TLweSampleFFTa* src, int N) {
    lagrangeCopy(&dst->a[0], &src->a[0], N);
    lagrangeCopy(dst->b, src->b, N);
}

static inline uint64_t torus_to_uint_round(Torus64 t, int MM1) {
    const int shift = 64 - MM1;
    if (shift == 0) return (uint64_t)t;
    const Torus64 half = (Torus64(1) << (shift - 1));
    const uint64_t mask = (MM1 == 64) ? ~UINT64_C(0) : ((UINT64_C(1) << MM1) - 1);
    return (uint64_t)((t + half) >> shift) & mask;
}
static inline Torus64 uint_to_torus_center(uint64_t u, int MM1) {
    const uint64_t p = (MM1 == 64) ? ~UINT64_C(0) : (UINT64_C(1) << MM1);
    const uint64_t half = (MM1 == 0) ? 0 : (p >> 1);
    int64_t s = (int64_t)u;
    if (MM1 != 64 && u >= half) s -= (int64_t)p;
    const int shift = 64 - MM1;
    return (Torus64)(s << shift);
}

static inline void print_record_bytes_from_torus(const char* label,
                                                 const Torus64* coeffs,
                                                 int coeff_count,
                                                 int stride,
                                                 int MM1,
                                                 uint64_t max_bytes) {
    std::cout << label;
    const uint64_t mask = (MM1 >= 64) ? ~UINT64_C(0) : ((UINT64_C(1) << MM1) - 1);
    uint64_t acc = 0, printed = 0;
    int accbits = 0;
    for (int i = 0; i < coeff_count && printed < max_bytes; i += stride) {
        acc |= (torus_to_uint_round(coeffs[i], MM1) & mask) << accbits;
        accbits += MM1;
        while (accbits >= 8 && printed < max_bytes) {
            std::cout << (int)(acc & 0xFF) << " ";
            acc >>= 8; accbits -= 8; ++printed;
        }
    }
    std::cout << "\n";
}

static int rs_env_int(const char* key, int dflt) {
    const char* v = getenv(key);
    return v ? atoi(v) : dflt;
}

static const int    RS_ENABLE    = rs_env_int("RS_ENABLE",      1);
static const int    SEED_ENABLE  = rs_env_int("SEED_ENABLE",    1);
static const int    RS_N0        = RS_ENABLE ? rs_env_int("RS_N0", 512) : 2048;
static const int    RS_LOG_Q2    = rs_env_int("RS_LOG_Q2",      24);
static const int    RS_LOG_Q3    = rs_env_int("RS_LOG_Q3",       6);
static const int    RS_LOG_SIGMA = rs_env_int("RS_LOG_SIGMA",    8);
static const int    RS_BGBIT     = rs_env_int("RS_BGBIT", 1);
static const int    RS_L         = RS_LOG_Q2 / RS_BGBIT;
static const double RS_SK0_WIDTH = pow(2., RS_LOG_SIGMA);

static const double RS_ALPHA     = pow(2., RS_LOG_SIGMA - RS_LOG_Q2);

static inline Torus64 round_to_q(Torus64 x, int log_q) {
    if (log_q >= 64) return x;
    const int sh        = 64 - log_q;
    const uint64_t half = UINT64_C(1) << (sh - 1);
    const uint64_t mask = ~((UINT64_C(1) << sh) - 1);
    return (Torus64)((((uint64_t)x + half) & mask));
}

static inline void modswitch_poly(Torus64Polynomial* p, int n, int log_q) {
    for (int i = 0; i < n; i++) p->coefs[i] = round_to_q(p->coefs[i], log_q);
}

static void decomp_rs(IntPolynomiala* out, const Torus64Polynomial* in,
                      int n, int bgbit, int l) {
    const uint64_t Bg     = UINT64_C(1) << bgbit;
    const uint64_t mask   = Bg - 1;
    const int64_t  halfBg = (int64_t)(Bg / 2);
    uint64_t offset = 0;
    for (int i = 0; i <= l; i++) offset |= (UINT64_C(1) << (63 - i * bgbit));

    for (int j = 0; j < n; j++) {
        const uint64_t v = (uint64_t)in->coefs[j] + offset;
        for (int p = 0; p < l; p++) {
            const int decal = 64 - (p + 1) * bgbit;
            out[p].coefs[j] = (int)((int64_t)((v >> decal) & mask) - halfBg);
        }
    }
}

static inline void kappa_embed(IntPolynomiala* out, const IntPolynomiala* s0, int N_big, int N0) {
    const int k = N_big / N0;
    for (int i = 0; i < N_big; i++) out->coefs[i] = 0;
    for (int j = 0; j < N0; j++) out->coefs[j * k] = s0->coefs[j];
}

static void tlwePhaseSmall(Torus64* phase, const TLweSample64* ct_small, const IntPolynomiala* s0, int N0) {

    for (int i = 0; i < N0; i++) phase[i] = ct_small->b->coefs[i];

    const Torus64* a = ct_small->a[0].coefs;
    for (int ai = 0; ai < N0; ai++) {
        const Torus64 av = a[ai];
        if (av == 0) continue;
        for (int sj = 0; sj < N0; sj++) {
            const int s = s0->coefs[sj];
            if (!s) continue;
            const Torus64 as = av * (Torus64)s;
            const int t = ai + sj;
            if (t < N0) phase[t] -= as;
            else        phase[t - N0] += as;
        }
    }
}

static TLweSample64** build_ringswitch_key(const Globals* env,
                                          const IntPolynomiala* s_big,
                                          const IntPolynomiala* rep_s0,
                                          int l, int bgbit,
                                          double stdev) {
    const int N = env->N;
    TLweSample64** rswk = new TLweSample64*[l];

    for (int p = 0; p < l; p++) {
        rswk[p] = new TLweSample64(N);

        for (int j = 0; j < N; j++)
            rswk[p]->a[0].coefs[j] = round_to_q((Torus64)random_int64(), RS_LOG_Q2);

        for (int j = 0; j < N; j++) rswk[p]->b->coefs[j] = (Torus64)random_gaussian64(0, stdev);

        torus64PolynomialMultAddKaratsuba_lvl2(rswk[p]->b, rep_s0, &rswk[p]->a[0], env);

        const int shift = 64 - (p + 1) * bgbit;
        const Torus64 Lp = (shift >= 0 && shift < 64) ? (Torus64)(UINT64_C(1) << shift) : 0;
        for (int j = 0; j < N; j++) {
            if (s_big->coefs[j]) rswk[p]->b->coefs[j] += Lp;
        }

        modswitch_poly(rswk[p]->b, N, RS_LOG_Q2);
    }
    return rswk;
}

static void ringSwitch_N_to_N0(TLweSample64* out_small,
                              const TLweSample64* in_big,
                              TLweSample64** rswk,
                              int N0,
                              const Globals* env) {
    const int N = env->N;
    const int k = N / N0;

    Torus64Polynomial a_q2(N);
    Torus64Polynomial b_prime(N);
    for (int i = 0; i < N; i++) a_q2.coefs[i]    = round_to_q(in_big->a[0].coefs[i], RS_LOG_Q2);
    for (int i = 0; i < N; i++) b_prime.coefs[i] = round_to_q(in_big->b->coefs[i],   RS_LOG_Q2);

    IntPolynomiala* adec = new_array1<IntPolynomiala>(RS_L, N);
    decomp_rs(adec, &a_q2, N, RS_BGBIT, RS_L);

    Torus64Polynomial a_prime(N);
    torusPolyClear(&a_prime, N);

    Torus64Polynomial tmp(N);

    for (int p = 0; p < RS_L; p++) {
        torus64PolynomialMultKaratsuba_lvl2(&tmp, &adec[p], &rswk[p]->a[0], env);
        torusPolySubTo(&a_prime, &tmp, N);

        torus64PolynomialMultKaratsuba_lvl2(&tmp, &adec[p], rswk[p]->b, env);
        torusPolySubTo(&b_prime, &tmp, N);
    }

    for (int i = 0; i < N0; i++) {
        out_small->a[0].coefs[i] = round_to_q(a_prime.coefs[i * k], RS_LOG_Q2);
        out_small->b->coefs[i]   = round_to_q(b_prime.coefs[i * k], RS_LOG_Q3);
    }

    delete_array1(adec);
}

const int Globals::k = 1;
const int Globals::N = 2048;
const int Globals::t = 12;
const int Globals::smalln = N;
const int Globals::basebit = 4;

static const int      kDefaultMM1          = 4;
static const uint64_t kDefaultNRecords     = (UINT64_C(1) << 22);
static const uint64_t kDefaultRowPrime     = 0;
static const uint64_t kDefaultColPrime     = 0;
static const bool     kDefaultRandomTarget = true;
static const uint64_t kDefaultTargetIndex  = 0;
static const int      kDefaultBgbit        = 9;
static const int      kDefaultL            = 6;

using namespace std;

Globals::Globals(int bgbit_in, int l_in) : bgbit(bgbit_in), l(l_in) {
    if (bgbit <= 0 || l <= 0 || (int64_t)bgbit * (int64_t)l > 63) {
        std::cerr << "[ERROR] invalid gadget params: bgbit=" << bgbit
                  << ", l=" << l << ". Require bgbit>0, l>0, and bgbit*l <= 63.\n";
        std::exit(1);
    }

    torusDecompOffset = 0;
    for (int i = 0; i <= l; ++i) torusDecompOffset |= (UINT64_C(1) << (63 - i * bgbit));

    torusDecompBuf = new uint64_t[N];

    lwekey = new int[N];
    for (int i = 0; i < N; ++i) lwekey[i] = random_bit();

    tlwekey = new IntPolynomiala(N);
    for (int i = 0; i < N; ++i) tlwekey->coefs[i] = lwekey[i];

    in_key = new int64_t[N + 1];
    for (int i = 0; i < N; ++i) in_key[i] = lwekey[i];
    in_key[N] = -1;
}

static inline size_t bytes_TLweSample64(int k, int N) {
    return (size_t)(k + 1) * (size_t)N * sizeof(int64_t);
}

static inline size_t bytes_TGswSample64(int k, int l, int N) {
    return (size_t)(k + 1) * (size_t)l * bytes_TLweSample64(k, N);
}

static inline size_t bytes_AutoKsKey64(int k, int l, int N) {
    return (size_t)l * bytes_TLweSample64(k, N);
}

static inline size_t bytes_AutoKsKeyFFTa(int l, int N) {
    return (size_t)2 * AKS_LIMBS * l * (size_t)N * sizeof(double);
}

static inline double MB10(size_t bytes) { return (double)bytes / 1e6; }

using Clock = std::chrono::high_resolution_clock;
static inline double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}
static inline void print_ms(const char* label, const Clock::time_point& a, const Clock::time_point& b) {
    const double ms = ms_between(a, b);
    std::cout << "[TIME] " << std::left << std::setw(38) << label
              << " : " << std::right << std::fixed << std::setprecision(3) << ms << " ms"
              << " (" << std::setprecision(3) << (ms / 1000.0) << " s)\n";
}

static inline bool is_power_of_two_u64(uint64_t x) {
    return x && ((x & (x - 1)) == 0);
}
static inline int log2_u64_exact(uint64_t x) {
    int r = 0;
    while (x > 1) { x >>= 1; r++; }
    return r;
}
static inline uint64_t ceil_div_u64(uint64_t a, uint64_t b) {
    return (a + b - 1) / b;
}
static inline uint64_t next_pow2_u64(uint64_t x) {
    if (x <= 1) return 1;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return x + 1;
}
static inline std::string pow2_desc_u64(uint64_t x, const char* unit = nullptr) {
    std::ostringstream oss;
    const bool have_unit = (unit && unit[0] != '\0');
    if (is_power_of_two_u64(x)) {
        oss << "2^" << log2_u64_exact(x);
        if (have_unit) oss << " " << unit;
        oss << " (=" << x;
        if (have_unit) oss << " " << unit;
        oss << ")";
        return oss.str();
    }
    oss << x;
    if (have_unit) oss << " " << unit;
    return oss.str();
}

static inline bool parse_u64_with_pow2(const char* arg, uint64_t* out) {
    if (!arg || !out) return false;
    const std::string s(arg);
    if (s.empty()) return false;

    if (s.size() >= 3 && s[0] == '2' && s[1] == '^') {
        size_t pos = 0;
        const unsigned long long exp = std::stoull(s.substr(2), &pos);
        if (pos != s.size() - 2) return false;
        if (exp >= 64) return false;
        *out = (UINT64_C(1) << exp);
        return true;
    }

    size_t pos = 0;
    const unsigned long long val = std::stoull(s, &pos);
    if (pos != s.size()) return false;
    *out = (uint64_t)val;
    return true;
}

static inline uint64_t parse_u64_or_die(const char* arg, const char* name) {
    uint64_t out = 0;
    try {
        if (parse_u64_with_pow2(arg, &out)) return out;
    } catch (const std::exception&) {
    }
    std::cerr << "[ERROR] invalid " << name << ": '" << (arg ? arg : "")
              << "'. Use a nonnegative integer or 2^k form.\n";
    std::exit(1);
}

static inline bool is_random_target_token(const char* arg) {
    return arg &&
           (std::strcmp(arg, "random") == 0 ||
            std::strcmp(arg, "rand") == 0 ||
            std::strcmp(arg, "r") == 0 ||
            std::strcmp(arg, "-1") == 0);
}

static inline uint64_t choose_random_target_index(uint64_t n, uint64_t* seed_out) {
    const uint64_t now_seed =
        (uint64_t)std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::random_device rd;
    const uint64_t seed = now_seed ^ ((uint64_t)rd() << 32) ^ (uint64_t)rd();
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);
    if (seed_out) *seed_out = seed;
    return dist(rng);
}

static inline uint64_t conservative_safe_row_cap(int MM1, double alpha1, uint64_t N) {
    (void)alpha1;

    if (MM1 == 4)  return std::min<uint64_t>(N, 128);
    if (MM1 == 16) return std::min<uint64_t>(N, 32);
    return std::min<uint64_t>(N, 64);
}

static void torusPolyMulByXNegJ(Torus64* out, const Torus64* in, int N, int j) {
    j %= N;
    if (j < 0) j += N;
    if (j == 0) {
        std::memcpy(out, in, sizeof(Torus64) * N);
        return;
    }
    for (int i = 0; i < N; i++) {
        const int src = i + j;
        if (src < N) out[i] = in[src];
        else         out[i] = -in[src - N];
    }
}

static void tlweMulByXNegJ(TLweSample64* out, const TLweSample64* in, int N, int k, int j) {

    for (int u = 0; u < k; u++) {
        torusPolyMulByXNegJ(out->a[u].coefs, in->a[u].coefs, N, j);
    }

    torusPolyMulByXNegJ(out->b->coefs, in->b->coefs, N, j);

}

static inline void tlweCopy(TLweSample64* dst, const TLweSample64* src, int N, int k) {
    for (int q = 0; q <= k; ++q) {
        std::memcpy(dst->a[q].coefs, src->a[q].coefs, sizeof(Torus64) * (size_t)N);
    }
}

static inline void tlweFFTTotorus(TLweSample64* out, const TLweSampleFFTa* in, int N, int k, const Globals* env) {
    for (int q = 0; q <= k; ++q) {
        TorusPolynomial64_fft_lvl2(&out->a[q], &in->a[q], env);
    }
}

static inline double noise_budget_from_big_phase(const TLweSample64* ct,
                                                 const Torus64Polynomial* expected,
                                                 int MM1,
                                                 const Globals* env,
                                                 int* mismatch_out = nullptr,
                                                 uint64_t* maxerr_out = nullptr) {
    const int N = env->N;
    Torus64Polynomial phase(N);
    tLwe64Phase_lvl2(&phase, ct, env);

    uint64_t maxerr = 0;
    int mismatches = 0;
    for (int i = 0; i < N; ++i) {
        const Torus64 diff = phase.coefs[i] - expected->coefs[i];
        const uint64_t adiff = (uint64_t)((diff < 0) ? -(int64_t)diff : diff);
        if (adiff > maxerr) maxerr = adiff;
        if (torus_to_uint_round(phase.coefs[i], MM1) != torus_to_uint_round(expected->coefs[i], MM1))
            mismatches++;
    }

    if (mismatch_out) *mismatch_out = mismatches;
    if (maxerr_out) *maxerr_out = maxerr;

    if (maxerr == 0) return 64.0 - (double)MM1 - 1.0;
    const double bit_e = std::ceil(std::log2((double)maxerr));
    return 64.0 - (double)MM1 - bit_e - 1.0;
}

static inline double noise_budget_from_small_phase(const Torus64* phase,
                                                   const Torus64* expected,
                                                   int N0,
                                                   int MM1,
                                                   int* mismatch_out = nullptr,
                                                   uint64_t* maxerr_out = nullptr) {
    uint64_t maxerr = 0;
    int mismatches = 0;
    double sumsq = 0.0;
    for (int i = 0; i < N0; ++i) {
        const Torus64 diff = phase[i] - expected[i];
        const uint64_t adiff = (uint64_t)((diff < 0) ? -(int64_t)diff : diff);
        if (adiff > maxerr) maxerr = adiff;
        const double d = (double)(int64_t)diff;
        sumsq += d * d;
        if (torus_to_uint_round(phase[i], MM1) != torus_to_uint_round(expected[i], MM1))
            mismatches++;
    }

    printf("[ERR] final_stddev = %.6e  (log2 = %.3f)\n",
           std::sqrt(sumsq / (double)N0),
           std::log2(std::sqrt(sumsq / (double)N0) + 1e-300));

    if (mismatch_out) *mismatch_out = mismatches;
    if (maxerr_out) *maxerr_out = maxerr;

    if (maxerr == 0) return 64.0 - (double)MM1 - 1.0;
    const double bit_e = std::ceil(std::log2((double)maxerr));
    return 64.0 - (double)MM1 - bit_e - 1.0;
}

static inline void print_nb_unavailable(const char* label, const char* why) {
    std::cout << "[NB] " << label << " = unavailable (" << why << ")\n";
}

static inline void print_nb_diag(const char* label, double nb, int mismatches, uint64_t maxerr) {
    std::cout << "[NB] " << label
              << " = " << std::fixed << std::setprecision(3) << nb
              << "  (mismatch=" << mismatches
              << ", maxerr=" << maxerr;
    if (maxerr > 0) {
        std::cout << ", log2(maxerr)=" << std::ceil(std::log2((double)maxerr));
    }
    std::cout << ")\n";
}

int main(int argc, char** argv) {

    const bool profile = true;
    const auto t_program_start = Clock::now();

    size_t summary_query_col_bytes  = 0;
    size_t summary_query_row_bytes  = 0;
    size_t summary_convk_bytes      = 0;
    size_t summary_aks_bytes        = 0;
    size_t summary_aks_fft_bytes    = 0;
    size_t summary_rswk_bytes       = 0;
    size_t summary_evalkey_bytes    = 0;

    double summary_T_preprocess_ms  = 0.0;
    double summary_T_aks_ms         = 0.0;
    double summary_T_aks_fft_ms     = 0.0;
    double summary_T_rswk_gen_ms    = 0.0;
    double summary_T_convk_gen_ms   = 0.0;
    double summary_T_offline_comp_ms = 0.0;
    double summary_T_row_expand_ms  = 0.0;
    double summary_T_row_query_unpack_ms = 0.0;
    double summary_T_first_dim_ms   = 0.0;
    double summary_T_col_unpack_ms  = 0.0;
    double summary_T_col_query_unpack_ms = 0.0;
    double summary_T_slot_query_prepare_ms = 0.0;
    double summary_T_local_unpack_ms = 0.0;
    double summary_T_second_dim_ms  = 0.0;
    double summary_T_blind_rot_ms   = 0.0;
    double summary_T_local_ext_ms   = 0.0;
    double summary_T_ring_switch_ms = 0.0;
    double summary_T_col_trace_ms    = 0.0;
    double summary_T_query_enc_ms    = 0.0;
    double summary_T_query_unpack_ms = 0.0;
    double summary_T_main_comp_ms   = 0.0;
    double summary_T_extract_time_ms = 0.0;
    double summary_T_total_ms       = 0.0;
    double summary_T_total_no_pre_ms = 0.0;

    double summary_nb_after_row_expand    = 0.0;
    double summary_nb_after_first_dim     = 0.0;
    double summary_nb_after_second_dim    = 0.0;
    double summary_nb_after_local_ext     = 0.0;
    double summary_nb_before_ring_switch  = 0.0;
    double summary_nb_final               = 0.0;

    size_t summary_resp_before_bytes = 0;
    size_t summary_resp_after_bytes  = 0;

    size_t summary_query_local_bytes = 0;

    double summary_output_noise_budget = 0.0;
    bool   summary_decryption_failure  = false;
    bool   summary_record_match        = false;
    int    summary_record_bad_coeffs   = 0;
    int    summary_record_first_bad    = -1;
    uint64_t summary_record_first_exp  = 0;
    uint64_t summary_record_first_got  = 0;

    int      summary_N         = 0;
    int      summary_logN      = 0;
    int      summary_MM1       = 0;
    uint64_t summary_MM        = 0;
    uint64_t summary_n         = 0;
    uint64_t summary_D1        = 0;
    uint64_t summary_M         = 0;
    uint64_t summary_row_prime = 0;
    uint64_t summary_col_prime = 0;

    double alpha1 = pow(2., -55);

    const uint64_t record_bytes_assumed = 256;
    int MM1 = (argc >= 2 ? std::stoi(argv[1]) : kDefaultMM1);

    if (MM1 <= 0 || MM1 >= 63 || (MM1 != 4 && (MM1 % 8) != 0)) {
        std::cerr << "[ERROR] MM1 must be 4 or in {8,16,24,32,40,48,56} (got " << MM1 << ").\n";
        return 1;
    }

    const int bgbit_arg = (argc >= 7 ? (int)parse_u64_or_die(argv[6], "bgbit") : kDefaultBgbit);
    const int l_arg     = (argc >= 8 ? (int)parse_u64_or_die(argv[7], "l") : kDefaultL);
    if (bgbit_arg <= 0 || l_arg <= 0 || (int64_t)bgbit_arg * (int64_t)l_arg > 63) {
        std::cerr << "[ERROR] invalid gadget params: bgbit=" << bgbit_arg
                  << ", l=" << l_arg << ". Require bgbit>0, l>0, and bgbit*l <= 63.\n";
        return 1;
    }

    const int      torus_shift = 64 - MM1;
    const uint64_t plain_mod   = UINT64_C(1) << MM1;
    const uint64_t plain_half  = plain_mod >> 1;
    const uint64_t plain_mask  = plain_mod - 1;
    const uint64_t record_bits = record_bytes_assumed * 8;

    if (record_bits % (uint64_t)MM1 != 0) {
        std::cerr << "[ERROR] record of " << record_bits
                  << " bits is not divisible by MM1=" << MM1 << "\n";
        return 1;
    }

    const uint64_t N1_target = record_bits / (uint64_t)MM1;
    if ((uint64_t)Globals::N % N1_target != 0) {
        std::cerr << "[ERROR] N=" << Globals::N << " not divisible by N1_target=" << N1_target << "\n";
        return 1;
    }
    const uint64_t D1 = (uint64_t)Globals::N / N1_target;

    uint64_t n = (argc >= 3 ? parse_u64_or_die(argv[2], "n") : kDefaultNRecords);
    if (n == 0) {
        std::cerr << "[ERROR] n must be positive.\n";
        return 1;
    }

    const uint64_t M = ceil_div_u64(n, D1);

    uint64_t row_prime = (argc >= 4 ? parse_u64_or_die(argv[3], "row_prime") : kDefaultRowPrime);
    uint64_t col_prime = (argc >= 5 ? parse_u64_or_die(argv[4], "col_prime") : kDefaultColPrime);
    const uint64_t requested_row_prime = row_prime;
    const uint64_t requested_col_prime = col_prime;

    if (row_prime == 0 && col_prime == 0) {
        col_prime = next_pow2_u64((uint64_t)std::ceil(std::sqrt((long double)M)));
        row_prime = next_pow2_u64(ceil_div_u64(M, col_prime));
    } else if (row_prime == 0) {
        if (!is_power_of_two_u64(col_prime)) {
            std::cerr << "[ERROR] col_prime=" << col_prime << " must be a power of 2.\n";
            return 1;
        }
        row_prime = next_pow2_u64(ceil_div_u64(M, col_prime));
    } else if (col_prime == 0) {
        if (!is_power_of_two_u64(row_prime)) {
            std::cerr << "[ERROR] row_prime=" << row_prime << " must be a power of 2.\n";
            return 1;
        }
        col_prime = next_pow2_u64(ceil_div_u64(M, row_prime));
    }

    if (!is_power_of_two_u64(row_prime) || row_prime == 0) {
        std::cerr << "[ERROR] row_prime=" << row_prime << " must be a power of 2.\n";
        return 1;
    }
    if (!is_power_of_two_u64(col_prime) || col_prime == 0) {
        std::cerr << "[ERROR] col_prime=" << col_prime << " must be a power of 2.\n";
        return 1;
    }

    const uint64_t safe_row_cap = conservative_safe_row_cap(MM1, alpha1, (uint64_t)Globals::N);
    bool row_prime_stabilized = false;
    if (row_prime > safe_row_cap) {
        row_prime = safe_row_cap;
        col_prime = next_pow2_u64(ceil_div_u64(M, row_prime));
        row_prime_stabilized = true;
        std::cerr << "[WARN] row' was reduced to conservative safe cap "
                  << pow2_desc_u64(row_prime)
                  << " under fixed alpha1; col' recomputed to "
                  << pow2_desc_u64(col_prime) << ".\n";
    }

    const unsigned __int128 layout_blocks_wide =
        (unsigned __int128)row_prime * (unsigned __int128)col_prime;
    if (layout_blocks_wide < (unsigned __int128)M) {
        std::cerr << "[ERROR] row' * col' = " << (unsigned long long)(row_prime * col_prime)
                  << " is smaller than logical packed blocks M=" << M << "\n";
        return 1;
    }
    const uint64_t layout_blocks = (uint64_t)layout_blocks_wide;
    if (row_prime > (uint64_t)Globals::N) {
        std::cerr << "[ERROR] row_prime=" << row_prime
                  << " exceeds N=" << Globals::N
                  << " but row compression A stores one packed selector over row' coefficients.\n";
        return 1;
    }

    summary_n         = n;
    summary_D1        = D1;
    summary_M         = M;
    summary_row_prime = row_prime;
    summary_col_prime = col_prime;

    uint64_t index = kDefaultTargetIndex;
    bool target_is_random = false;
    uint64_t target_seed = 0;
    if (argc >= 6) {
        if (is_random_target_token(argv[5])) {
            target_is_random = true;
            index = choose_random_target_index(n, &target_seed);
        } else {
            index = parse_u64_or_die(argv[5], "target_index");
        }
    } else if (kDefaultRandomTarget) {
        target_is_random = true;
        index = choose_random_target_index(n, &target_seed);
    }
    if (index >= n) index = 0;

    const uint64_t target_b = index / D1;
    const uint64_t target_s = index % D1;
    const uint64_t target_r = target_b / col_prime;
    const uint64_t target_c = target_b % col_prime;

    std::cout << "[TARGET] mode=" << (target_is_random ? "random" : "fixed");
    if (target_is_random) std::cout << "  seed=" << target_seed;
    std::cout << "  index=" << index
              << "  block=" << target_b
              << "  (r=" << target_r << ", c=" << target_c << ", s=" << target_s << ")\n";
    std::cout << "[LAYOUT] n=" << pow2_desc_u64(n)
              << "  D1=" << pow2_desc_u64(D1)
              << "  M=" << pow2_desc_u64(M)
              << "  row'=" << pow2_desc_u64(row_prime)
              << "  col'=" << pow2_desc_u64(col_prime)
              << "  row'*col'=" << layout_blocks
              << "  padding=" << (layout_blocks - M) << "\n";
    if (row_prime_stabilized) {
        std::cout << "[STABILITY] requested row'="
                  << (requested_row_prime ? pow2_desc_u64(requested_row_prime) : std::string("auto"))
                  << ", requested col'="
                  << (requested_col_prime ? pow2_desc_u64(requested_col_prime) : std::string("auto"))
                  << " -> using row'=" << pow2_desc_u64(row_prime)
                  << ", col'=" << pow2_desc_u64(col_prime)
                  << " to avoid decryption failure under fixed alpha1.\n";
    }

    const auto t_env0 = Clock::now();
    Globals* env = new Globals(bgbit_arg, l_arg);
    const auto t_env1 = Clock::now();
    if (profile) print_ms("Init Globals()", t_env0, t_env1);

    const int32_t k    = env->k;
    const int32_t N    = env->N;
    const int32_t l    = env->l;
    const int     Bgbit = env->bgbit;
    const uint64_t N1  = (uint64_t)N / D1;

    const int      N0  = RS_N0;

    summary_N    = N;
    summary_logN = (int)std::log2((double)N);
    summary_MM1  = MM1;
    summary_MM   = plain_mod;

    const int log2col = log2_u64_exact(col_prime);
    const int log2row = log2_u64_exact(row_prime);
    const int logD1   = log2_u64_exact(D1);

    std::cout << "[DIMS] log2(col')=" << log2col
              << "  log2(row')=" << log2row
              << "  logD1=" << logD1 << "\n"
              << "       CMux depth: col'=" << col_prime
              << " (depth " << log2col << " vs old depth "
              << (int)std::ceil(std::log2((double)(n/D1))) << ")\n";

    TLweSample64*   enc_dataa = new TLweSample64(N);

    {
        const uint64_t bytes_DB_raw      = layout_blocks * (uint64_t)N * sizeof(int);
        const uint64_t bytes_DB_torus    = (uint64_t)N * sizeof(Torus64);
        const uint64_t bytes_PackedDB    = layout_blocks * (uint64_t)N * sizeof(double);
        const uint64_t bytes_row_selects = row_prime * (uint64_t)(k+1) * (uint64_t)N * sizeof(double);
        const uint64_t bytes_enc_temp_m  = col_prime * (uint64_t)(k+1) * (uint64_t)N * sizeof(double);
        const uint64_t bytes_enc_dataa_m = (uint64_t)(k+1) * (uint64_t)N * sizeof(int64_t);
        std::cerr << "\n[MEMORY ESTIMATE (model)]\n";
        std::cerr << "  DB_raw [layout polys, int]         : "
                  << pretty_bytes(bytes_DB_raw) << "\n";
        std::cerr << "  DB_torus [target block only]       : "
                  << pretty_bytes(bytes_DB_torus) << "\n";
        std::cerr << "  PackedDBFFT [row'*col'=" << layout_blocks << " polys, FFT]: "
                  << pretty_bytes(bytes_PackedDB) << "\n";
        std::cerr << "  subtotal (DB preprocess buffers)   : "
                  << pretty_bytes(bytes_DB_raw + bytes_DB_torus + bytes_PackedDB) << "\n";
        std::cerr << "  row selectors FFT [row'=" << row_prime << "]:   "
                  << pretty_bytes(bytes_row_selects) << "\n";
        std::cerr << "  enc_temp (CMux scratch) [col'=" << col_prime << "]: "
                  << pretty_bytes(bytes_enc_temp_m) << "\n";
        std::cerr << "  enc_dataa (coef):                               "
                  << pretty_bytes(bytes_enc_dataa_m) << "\n";
    }

    std::vector<uint8_t> records((size_t)n * (size_t)record_bytes_assumed);
    for (size_t z = 0; z < records.size(); ++z)
        records[z] = (uint8_t)(random_int64() & 0xFF);

    const auto t_pre0 = Clock::now();

    IntPolynomiala*    DB_raw   = new_array1<IntPolynomiala>(layout_blocks, N);
    Torus64Polynomial* DB_torus_tgt = new Torus64Polynomial(N);
    for (int j = 0; j < N; ++j) DB_torus_tgt->coefs[j] = 0;

    const auto t_pre_fill0 = Clock::now();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int64_t bsigned = 0; bsigned < (int64_t)layout_blocks; ++bsigned) {
        const uint64_t b = (uint64_t)bsigned;
        const bool is_target_block = (b == target_b);
        for (int j = 0; j < N; ++j) DB_raw[b].coefs[j] = 0;
        if (b >= M) continue;
        for (uint64_t s = 0; s < D1; ++s) {
            const uint64_t rec = b * D1 + s;
            if (rec >= n) continue;
            const uint8_t* src = &records[(size_t)rec * (size_t)record_bytes_assumed];
            uint64_t acc = 0; int accbits = 0; size_t bp = 0;
            for (uint64_t i = 0; i < N1; ++i) {
                while (accbits < MM1) { acc |= (uint64_t)src[bp++] << accbits; accbits += 8; }
                const uint64_t u_raw = acc & plain_mask;
                acc >>= MM1; accbits -= MM1;
                int64_t msg = (int64_t)u_raw;
                if (u_raw >= plain_half) msg -= (int64_t)plain_mod;
                const uint64_t idx = s + D1 * i;
                DB_raw[b].coefs[idx] = (int)msg;
                if (is_target_block)
                    DB_torus_tgt->coefs[idx] = uint_to_torus_center(u_raw, MM1);
            }
        }
    }
    const auto t_pre_fill1 = Clock::now();
    print_ms("Preprocess: DB reshape", t_pre_fill0, t_pre_fill1);

    LagrangeHalfCPolynomiala* PackedDBFFT = new_array1<LagrangeHalfCPolynomiala>(layout_blocks, N);
    const auto t_pre_fft0 = Clock::now();
#if defined(_OPENMP) && defined(USE_FFT)

#pragma omp parallel
    {
        FFT_Processor_Spqlios local_fft(N);
#pragma omp for schedule(static)
        for (int64_t bsigned = 0; bsigned < (int64_t)layout_blocks; ++bsigned)
            local_fft.execute_reverse_int(PackedDBFFT[bsigned].values,
                                          DB_raw[bsigned].coefs);
    }
#else
    for (uint64_t b = 0; b < layout_blocks; ++b)
        IntPolynomial_ifft_lvl2(&PackedDBFFT[b], &DB_raw[b], env);
#endif
    const auto t_pre_fft1 = Clock::now();
    print_ms("Preprocess: DB FFT", t_pre_fft0, t_pre_fft1);

    delete_array1(DB_raw);

    const auto t_pre1 = Clock::now();
    summary_T_preprocess_ms = ms_between(t_pre0, t_pre1);
    printf("[TIME] preprocessing   = %.3f ms\n", summary_T_preprocess_ms);

    const int n_trace = 1;
    const auto t_aks0 = Clock::now();
    AutoKsKey64** aks_list = AutoKsKeyGenAll(env, n_trace);
    const auto t_aks1 = Clock::now();
    summary_T_aks_ms = ms_between(t_aks0, t_aks1);
    if (profile) print_ms("AutoKsKeyGenAll (aks_list)", t_aks0, t_aks1);

    const auto t_aks_fft0 = Clock::now();
    AutoKsKeyFFTa** aks_fft_list = AutoKsKeyAllToFFT(aks_list, env, n_trace);
    const auto t_aks_fft1 = Clock::now();
    summary_T_aks_fft_ms = ms_between(t_aks_fft0, t_aks_fft1);
    if (profile) print_ms("AutoKsKeyAllToFFT (FFT-AutoKS)", t_aks_fft0, t_aks_fft1);
    std::cout << "[TRACE] RevHomTrace backend = FFT-AutoKS\n";

    IntPolynomiala* minus_sk = new IntPolynomiala(N);
    for (int i = 0; i < N; i++)
        minus_sk->coefs[i] = -env->tlwekey->coefs[i];

    TGswSample64* convk = new TGswSample64(l, N);

    summary_resp_before_bytes = bytes_TLweSample64(k, N);

    summary_resp_after_bytes  = ((size_t)N0 * RS_LOG_Q2 + (size_t)N1 * RS_LOG_Q3 + 7) / 8;

    IntPolynomiala* sk0 = new IntPolynomiala(N0);
    std::mt19937 rng(std::random_device{}());
    if (!RS_ENABLE) {

        for (int i = 0; i < N0; i++) sk0->coefs[i] = env->tlwekey->coefs[i];
    } else {

        std::normal_distribution<double> g(0.0, RS_SK0_WIDTH);
        for (int i = 0; i < N0; i++) sk0->coefs[i] = (int)llround(g(rng));
    }
    IntPolynomiala* rep_sk0 = new IntPolynomiala(N);
    kappa_embed(rep_sk0, sk0, N, N0);

    const double rswk_stdev = RS_ALPHA;
    const auto t_rswk0 = Clock::now();
    TLweSample64** rswk = build_ringswitch_key(env, env->tlwekey, rep_sk0, RS_L, RS_BGBIT, rswk_stdev);
    const auto t_rswk1 = Clock::now();
    summary_T_rswk_gen_ms = ms_between(t_rswk0, t_rswk1);
    if (profile) print_ms("Gen rswk (RingSwitch key)", t_rswk0, t_rswk1);
    summary_rswk_bytes = RS_ENABLE ? (size_t)RS_L * 2 * (size_t)N * sizeof(Torus64) : 0;

    const auto t_convk0 = Clock::now();
    tGsw64Encrypt_poly_2(convk, minus_sk, pow(2., -55), env);
    const auto t_convk1 = Clock::now();
    summary_T_convk_gen_ms = ms_between(t_convk0, t_convk1);
    if (profile) print_ms("Gen convk (tGsw encrypt)", t_convk0, t_convk1);

    const int aks_count = summary_logN;
    summary_convk_bytes   = bytes_TGswSample64(k, l, N);
    summary_aks_bytes     = (size_t)aks_count * bytes_AutoKsKey64(k, l, N);
    summary_aks_fft_bytes = (size_t)aks_count * bytes_AutoKsKeyFFTa(l, N);
    summary_evalkey_bytes = summary_convk_bytes + summary_aks_bytes + summary_rswk_bytes;

    {
        const uint64_t b_dbg = target_r * col_prime + target_c;
        const uint64_t showB = std::min<uint64_t>(16, record_bytes_assumed);
        std::cout << "[DB] expected record (block=" << b_dbg << ", slot=" << target_s
                  << ") first " << showB << " bytes: ";
        const uint64_t mask = (UINT64_C(1) << MM1) - 1;
        uint64_t acc = 0, printed = 0;
        int accbits = 0;
        for (uint64_t i = 0; i < N1 && printed < showB; ++i) {
            const Torus64 tor = DB_torus_tgt->coefs[target_s + D1 * i];
            acc |= (torus_to_uint_round(tor, MM1) & mask) << accbits;
            accbits += MM1;
            while (accbits >= 8 && printed < showB) {
                std::cout << (int)(acc & 0xFF) << " ";
                acc >>= 8; accbits -= 8; ++printed;
            }
        }
        std::cout << "\n";
    }

    const auto t_online_start = Clock::now();

    const uint64_t COL_BASE   = row_prime;
    const uint64_t SLOT_BASE  = COL_BASE + (uint64_t)l * (uint64_t)log2col;
    const uint64_t QUERY_USED = SLOT_BASE + (uint64_t)l * (uint64_t)logD1;
    if (QUERY_USED > (uint64_t)N) {
        std::cerr << "[ERROR] unified query needs " << QUERY_USED
                  << " coefficients but N=" << N
                  << " (row'=" << row_prime << ", l=" << l
                  << ", log2col=" << log2col << ", logD1=" << logD1 << ").\n";
        return 1;
    }
    const uint64_t QUERY_SPAN   = next_pow2_u64(QUERY_USED);
    const int      QUERY_STRIDE = (int)((uint64_t)N / QUERY_SPAN);

    std::vector<int64_t> bit_c(log2col);
    int_to_bin_digit((int)target_c, log2col, bit_c.data());
    std::vector<int64_t> bit_j(logD1, 0);
    int_to_bin_digit((uint32_t)target_s, logD1, bit_j.data());

    const uint64_t query_seed = (uint64_t)random_int64();
    TLweSample64* query_ct = new TLweSample64(N);
    {
        Torus64Polynomial* query_msg = new Torus64Polynomial(N);
        for (int j = 0; j < N; ++j) query_msg->coefs[j] = 0;

        query_msg->coefs[target_r * QUERY_STRIDE] =
            (Torus64)(UINT64_C(1) << torus_shift) / (Torus64)QUERY_SPAN;

        for (int i = 0; i < l; ++i) {
            const int shift_i = 64 - (i + 1) * Bgbit;
            const Torus64 Li = (shift_i >= 0 && shift_i < 64)
                             ? (Torus64)(UINT64_C(1) << shift_i) : 0;
            for (int p = 0; p < log2col; ++p)
                query_msg->coefs[(COL_BASE + (uint64_t)i * log2col + p) * QUERY_STRIDE] =
                    bit_c[log2col - 1 - p] * Li;
            for (int p = 0; p < logD1; ++p)
                query_msg->coefs[(SLOT_BASE + (uint64_t)i * logD1 + p) * QUERY_STRIDE] =
                    bit_j[logD1 - 1 - p] * Li;
        }

        std::mt19937_64 prg(query_seed);
        for (int j = 0; j < N; ++j)
            query_ct->a[0].coefs[j] = (Torus64)prg();
        for (int j = 0; j < N; ++j)
            query_ct->b->coefs[j] = random_gaussian64(0, alpha1);
        torus64PolynomialMultAddKaratsuba_lvl2(query_ct->b, &env->tlwekey[0],
                                               &query_ct->a[0], env);
        for (int j = 0; j < N; ++j)
            query_ct->b->coefs[j] += query_msg->coefs[j];

        delete query_msg;
    }

    summary_query_row_bytes   = SEED_ENABLE
        ? sizeof(query_seed) + (size_t)N * sizeof(Torus64)
        : 2 * (size_t)N * sizeof(Torus64);
    summary_query_col_bytes   = 0;
    summary_query_local_bytes = 0;

    std::cout << "[QUERY] unified: 1 RLWE, " << QUERY_USED << "/" << N
              << " coefficients used; transmitted as "
              << summary_query_row_bytes << " B\n";
    std::cout << "[QUERY] row selectors: coefficient tree over " << QUERY_SPAN << " coefficients, "
              << (QUERY_SPAN - 1) << " key switches\n";

    summary_T_query_enc_ms = ms_between(t_online_start, Clock::now());
    const auto t_row_sel0 = Clock::now();

    TLweSample64** q_ct = new TLweSample64*[row_prime];
    {
        std::vector<GlweCipher64*> leaf(QUERY_SPAN);
        for (auto& p : leaf) p = new GlweCipher64(N);
        SubringCoeffTree_fft(leaf.data(), query_ct, (int)QUERY_SPAN, QUERY_STRIDE, 0, 1, aks_fft_list, env);
        for (uint64_t u = 0; u < QUERY_SPAN; ++u) {
            if (u < row_prime) q_ct[u] = leaf[u];
            else delete leaf[u];
        }
    }

    LagrangeHalfCPolynomiala* q_a_fft = new_array1<LagrangeHalfCPolynomiala>(row_prime, N);
    LagrangeHalfCPolynomiala* q_b_fft = new_array1<LagrangeHalfCPolynomiala>(row_prime, N);
    for (uint64_t u = 0; u < row_prime; ++u) {
        TorusPolynomial64_ifft_lvl2(&q_a_fft[u], &q_ct[u]->a[0], env);
        TorusPolynomial64_ifft_lvl2(&q_b_fft[u], q_ct[u]->b, env);
    }

    const auto t_row_sel1 = Clock::now();
    summary_T_row_expand_ms = ms_between(t_row_sel0, t_row_sel1);
    summary_T_row_query_unpack_ms = summary_T_row_expand_ms;
    printf("[TIME] row_expand      = %.3f ms\n", summary_T_row_expand_ms);

    const auto t_step6_0 = Clock::now();
    TLweSample64*** cipher_prime = new TLweSample64**[l];
    for (int i = 0; i < l; i++) {
        cipher_prime[i] = new TLweSample64*[log2col];
        for (int j = 0; j < log2col; j++) {
            cipher_prime[i][j] = new TLweSample64(N);
            GlweCipher64* ct_in  = new GlweCipher64(N);
            GlweCipher64* ct_out = new GlweCipher64(N);
            tlweMulByXNegJ(ct_in, query_ct, N, k,
                           (int)((COL_BASE + (uint64_t)i * log2col + j) * QUERY_STRIDE));
            RevHomTrace_Alg5_fft(ct_out, ct_in, aks_fft_list, n_trace, env);
            for (int q = 0; q <= k; ++q)
                for (int t = 0; t < N; ++t)
                    cipher_prime[i][j]->a[q].coefs[t] = ct_out->a[q].coefs[t];
            for (int t = 0; t < N; ++t)
                cipher_prime[i][j]->b->coefs[t] = ct_out->b->coefs[t];
            delete ct_in; delete ct_out;
        }
    }
    const auto t_step6_1 = Clock::now();
    summary_T_col_trace_ms = ms_between(t_step6_0, t_step6_1);
    if (profile) print_ms("Col query: RevHomTrace FFT-AutoKS", t_step6_0, t_step6_1);

    const auto t_first_dim0 = Clock::now();
    TLweSampleFFTa* Y_fft = new_array1<TLweSampleFFTa>(col_prime, N);
    for (uint64_t c = 0; c < col_prime; ++c) tlweFFTClear(&Y_fft[c], N, env);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (uint64_t c = 0; c < col_prime; ++c) {
        for (uint64_t u = 0; u < row_prime; ++u) {
            const uint64_t bidx = u * col_prime + c;
            LagrangeHalfCPolynomialAddMul_lvl2(&Y_fft[c].a[0], &q_a_fft[u], &PackedDBFFT[bidx], env);
            LagrangeHalfCPolynomialAddMul_lvl2(Y_fft[c].b, &q_b_fft[u], &PackedDBFFT[bidx], env);
        }
    }
    const auto t_first_dim1 = Clock::now();
    summary_T_first_dim_ms = ms_between(t_first_dim0, t_first_dim1);
    printf("[TIME] first_dim       = %.3f ms\n", summary_T_first_dim_ms);

    {
        Torus64Polynomial expected_row_sel(N);
        for (int j = 0; j < N; ++j) expected_row_sel.coefs[j] = 0;
        expected_row_sel.coefs[0] = (Torus64)(UINT64_C(1) << torus_shift);
        int row_expand_bad = 0;
        uint64_t row_expand_maxerr = 0;
        summary_nb_after_row_expand =
            noise_budget_from_big_phase(q_ct[target_r], &expected_row_sel, MM1, env,
                                        &row_expand_bad, &row_expand_maxerr);
        print_nb_diag("after_row_expand", summary_nb_after_row_expand, row_expand_bad, row_expand_maxerr);
    }

    TLweSample64* first_dim_target = new TLweSample64(N);
    tlweFFTTotorus(first_dim_target, &Y_fft[target_c], N, k, env);
    int first_dim_bad = 0;
    uint64_t first_dim_maxerr = 0;
    summary_nb_after_first_dim =
        noise_budget_from_big_phase(first_dim_target, DB_torus_tgt, MM1, env,
                                    &first_dim_bad, &first_dim_maxerr);
    print_nb_diag("after_first_dim", summary_nb_after_first_dim, first_dim_bad, first_dim_maxerr);

    TGswSample64* extract = new_array1<TGswSample64>(log2col, l, N);
    TGswSampleFFTa* extFFT = new_array1<TGswSampleFFTa>(log2col, l, N);
    int32_t kpl = (k + 1) * l;

    const auto t_unpack0 = Clock::now();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0; j < log2col; ++j) {
        TLweSample64** rlweInputs = new TLweSample64*[l];
        for (int i = 0; i < l; i++) rlweInputs[i] = cipher_prime[i][j];
        unpacking_algorithm4(&extract[j], const_cast<const TLweSample64**>(rlweInputs), convk, env);
        delete[] rlweInputs;
    }
    for (int s = 0; s < log2col; ++s) {
        for (int i = 0; i < kpl; ++i) {
            for (int q = 0; q <= k; ++q) {
                TorusPolynomial64_ifft_lvl2(&extFFT[s].allsamples[i].a[q],
                                            &extract[s].allsamples[i].a[q], env);
            }
        }
    }
    const auto t_unpack1 = Clock::now();
    summary_T_col_unpack_ms = ms_between(t_unpack0, t_unpack1);
    summary_T_col_query_unpack_ms = summary_T_col_unpack_ms;
    printf("[TIME] col_unpack      = %.3f ms\n", summary_T_col_unpack_ms);
    print_nb_unavailable("after_col_unpack", "TGSW selector stage");

    const auto t_second0 = Clock::now();
    TLweSampleFFTa* cur_nodes = new_array1<TLweSampleFFTa>(col_prime, N);
    TLweSampleFFTa* next_nodes = new_array1<TLweSampleFFTa>(col_prime, N);
    for (uint64_t c = 0; c < col_prime; ++c) {
        tlweFFTCopy(&cur_nodes[c], &Y_fft[c], N);
    }

    uint64_t cur_sz = col_prime;
    for (int b = 0; b < log2col; ++b) {
        const uint64_t next_sz = cur_sz >> 1;
        for (uint64_t i = 0; i < next_sz; ++i) {
            CMuxFFTa(&next_nodes[i], &extFFT[b], &cur_nodes[2 * i], &cur_nodes[2 * i + 1], env);
        }
        std::swap(cur_nodes, next_nodes);
        cur_sz = next_sz;
    }

    tlweFFTTotorus(enc_dataa, &cur_nodes[0], N, k, env);
    const auto t_second1 = Clock::now();
    summary_T_second_dim_ms = ms_between(t_second0, t_second1);
    printf("[TIME] second_dim_cmux = %.3f ms\n", summary_T_second_dim_ms);

    int second_dim_bad = 0;
    uint64_t second_dim_maxerr = 0;
    summary_nb_after_second_dim =
        noise_budget_from_big_phase(enc_dataa, DB_torus_tgt, MM1, env,
                                    &second_dim_bad, &second_dim_maxerr);
    print_nb_diag("after_second_dim", summary_nb_after_second_dim, second_dim_bad, second_dim_maxerr);

    const auto t_local0 = Clock::now();
    std::vector<TGswSample64*> j_ctrl(logD1, nullptr);

    TLweSample64*** cipher_prime_j = new TLweSample64**[l];

    const auto t_local_prepare0 = Clock::now();
    for (int i = 0; i < l; i++) {
        cipher_prime_j[i] = new TLweSample64*[logD1];
        for (int t = 0; t < logD1; t++) {
            cipher_prime_j[i][t] = new TLweSample64(N);
            GlweCipher64* ct_in = new GlweCipher64(N);
            GlweCipher64* ct_out = new GlweCipher64(N);
            tlweMulByXNegJ(ct_in, query_ct, N, k,
                           (int)((SLOT_BASE + (uint64_t)i * logD1 + t) * QUERY_STRIDE));
            RevHomTrace_Alg5_fft(ct_out, ct_in, aks_fft_list, n_trace, env);
            for (int u = 0; u < N; u++) {
                cipher_prime_j[i][t]->a[0].coefs[u] = ct_out->a[0].coefs[u];
                cipher_prime_j[i][t]->b->coefs[u]   = ct_out->b->coefs[u];
            }
            delete ct_in;
            delete ct_out;
        }
    }
    const auto t_local_prepare1 = Clock::now();
    summary_T_slot_query_prepare_ms = ms_between(t_local_prepare0, t_local_prepare1);
    printf("[TIME] slot_query_prepare = %.3f ms\n", summary_T_slot_query_prepare_ms);

    const auto t_local_unpack0 = Clock::now();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < logD1; b++) {
        j_ctrl[b] = new TGswSample64(l, N);
        std::vector<const TLweSample64*> tmp_list(l);
        for (int i = 0; i < l; i++) tmp_list[i] = cipher_prime_j[i][b];
        unpacking_algorithm4(j_ctrl[b], tmp_list.data(), convk, env);
    }
    const auto t_local_unpack1 = Clock::now();
    summary_T_local_unpack_ms = ms_between(t_local_unpack0, t_local_unpack1);
    printf("[TIME] slot_query_unpack = %.3f ms\n", summary_T_local_unpack_ms);

    TLweSample64* extracted_block = new TLweSample64(N);
    tlweCopy(extracted_block, enc_dataa, N, k);

    const auto t_blind0 = Clock::now();
    for (int b = 0; b < logD1; ++b) {
        TLweSample64* acc_rot = new TLweSample64(N);
        TLweSample64* acc_next = new TLweSample64(N);
        tlweMulByXNegJ(acc_rot, extracted_block, N, k, (1 << b));
        CMux(acc_next, j_ctrl[b], extracted_block, acc_rot, env);
        delete extracted_block;
        delete acc_rot;
        extracted_block = acc_next;
    }
    const auto t_blind1 = Clock::now();
    summary_T_blind_rot_ms = ms_between(t_blind0, t_blind1);
    printf("[TIME] blind_rotation  = %.3f ms\n", summary_T_blind_rot_ms);
    const auto t_local1 = Clock::now();
    summary_T_local_ext_ms = ms_between(t_local0, t_local1);
    printf("[TIME] local_extract   = %.3f ms\n", summary_T_local_ext_ms);

    Torus64Polynomial expected_extracted(N);
    torusPolyMulByXNegJ(expected_extracted.coefs, DB_torus_tgt->coefs, N, (int)target_s);
    int local_ext_bad = 0;
    uint64_t local_ext_maxerr = 0;
    summary_nb_after_local_ext =
        noise_budget_from_big_phase(extracted_block, &expected_extracted, MM1, env,
                                    &local_ext_bad, &local_ext_maxerr);
    summary_nb_before_ring_switch = summary_nb_after_local_ext;
    print_nb_diag("after_local_extract", summary_nb_after_local_ext, local_ext_bad, local_ext_maxerr);
    print_nb_diag("before_ring_switch", summary_nb_before_ring_switch, local_ext_bad, local_ext_maxerr);

    const auto t_rs0 = Clock::now();
    TLweSample64* outN0 = new TLweSample64(N0);
    if (RS_ENABLE) {
        ringSwitch_N_to_N0(outN0, extracted_block, rswk, N0, env);
    } else {

        for (int i = 0; i < N0; i++) {
            outN0->a[0].coefs[i] = round_to_q(extracted_block->a[0].coefs[i], RS_LOG_Q2);
            outN0->b->coefs[i]   = round_to_q(extracted_block->b->coefs[i],   RS_LOG_Q3);
        }
    }
    const auto t_rs1 = Clock::now();
    summary_T_ring_switch_ms = ms_between(t_rs0, t_rs1);
    printf("[TIME] ring_switch     = %.3f ms\n", summary_T_ring_switch_ms);

    const auto t_online_end = Clock::now();

    const int rs_stride = N / N0;
    Torus64* expected_small = new Torus64[N0];
    for (int i = 0; i < N0; i++)
        expected_small[i] = expected_extracted.coefs[rs_stride * i];

    Torus64* phase_a2 = new Torus64[N0];
    tlwePhaseSmall(phase_a2, outN0, sk0, N0);
    int bad_a2 = 0;
    uint64_t maxerr_a2 = 0;
    summary_nb_final = noise_budget_from_small_phase(phase_a2, expected_small, N0, MM1,
                                                     &bad_a2, &maxerr_a2);
    summary_output_noise_budget = summary_nb_final;
    summary_decryption_failure = (bad_a2 > 0);
    print_nb_diag("final", summary_nb_final, bad_a2, maxerr_a2);

    summary_record_bad_coeffs = 0;
    summary_record_first_bad = -1;
    const int rec_stride = N0 / (int)N1;
    for (int j = 0; j < (int)N1; ++j) {
        const int i = j * rec_stride;
        const uint64_t expected_u = torus_to_uint_round(expected_small[i], MM1);
        const uint64_t got_u = torus_to_uint_round(phase_a2[i], MM1);
        if (expected_u != got_u) {
            if (summary_record_first_bad < 0) {
                summary_record_first_bad = i;
                summary_record_first_exp = expected_u;
                summary_record_first_got = got_u;
            }
            summary_record_bad_coeffs++;
        }
    }
    summary_record_match = (summary_record_bad_coeffs == 0);

    int sent_bad_coeffs = 0;
    {
        TLweSample64 sent(N0);
        tlweCopy(&sent, outN0, N0, k);
        for (int i = 0; i < N0; i++)
            if (i % rec_stride) sent.b->coefs[i] = 0;
        Torus64* phase_sent = new Torus64[N0];
        tlwePhaseSmall(phase_sent, &sent, sk0, N0);
        for (int j = 0; j < (int)N1; ++j) {
            const int i = j * rec_stride;
            if (torus_to_uint_round(phase_sent[i], MM1) != torus_to_uint_round(expected_small[i], MM1))
                sent_bad_coeffs++;
        }
        delete[] phase_sent;
    }

    summary_T_col_query_unpack_ms += summary_T_col_trace_ms;
    summary_T_query_unpack_ms =
        summary_T_row_query_unpack_ms +
        summary_T_col_query_unpack_ms;
    summary_T_main_comp_ms =
        summary_T_first_dim_ms +
        summary_T_second_dim_ms;

    summary_T_extract_time_ms =
        summary_T_local_ext_ms +
        summary_T_ring_switch_ms;
    summary_T_total_ms = ms_between(t_online_start, t_online_end);
    summary_T_total_no_pre_ms = summary_T_total_ms;
    summary_T_offline_comp_ms =
        summary_T_preprocess_ms +
        summary_T_aks_ms +
        summary_T_aks_fft_ms +
        summary_T_rswk_gen_ms +
        summary_T_convk_gen_ms;

    std::cout << "\n==================== SUMMARY ====================\n";
    std::cout << "[PARAMETERS]\n";
    std::cout << "  - N (poly degree)         : " << summary_N << " (logN=" << summary_logN << ")\n";
    std::cout << "  - plaintext modulus p     : 2^" << summary_MM1 << " (p=" << summary_MM << ")\n";
    std::cout << "  - gadget (bgbit, l)       : (" << Bgbit << ", " << l << ")\n";
    std::cout << "  - logical records n       : " << pow2_desc_u64(summary_n) << "\n";
    std::cout << "  - data element bytes      : " << pow2_desc_u64(record_bytes_assumed, "bytes") << "\n";
    std::cout << "  - coeff payload bits      : " << summary_MM1 << " bits\n";
    std::cout << "  - packed block bytes      : "
              << pow2_desc_u64((uint64_t)summary_N * (uint64_t)summary_MM1 / 8, "bytes") << "\n";
    std::cout << "  - records / packed block  : " << summary_D1 << "\n";
    std::cout << "  - packed blocks (M)       : " << pow2_desc_u64(summary_M) << "\n";
    std::cout << "  - layout blocks           : " << pow2_desc_u64(layout_blocks) << "\n";
    std::cout << "  - padding blocks          : " << (layout_blocks - summary_M) << "\n";
    std::cout << "  - safe row' cap           : " << pow2_desc_u64(safe_row_cap) << "\n";
    std::cout << "  - row'                    : " << pow2_desc_u64(summary_row_prime) << "\n";
    std::cout << "  - col'                    : " << pow2_desc_u64(summary_col_prime) << "\n";
    std::cout << "  - target (r,c,s)          : (" << target_r << ", " << target_c << ", " << target_s << ")\n";

    const size_t summary_query_total_bytes =
        summary_query_row_bytes + summary_query_col_bytes + summary_query_local_bytes;
    std::cout << "\n[OFFLINE / COMM]\n";
    std::cout << "  - Query(row compressed A) : " << summary_query_row_bytes << " bytes (" << MB10(summary_query_row_bytes) << " MB)\n";
    std::cout << "  - Query(col compressed)   : " << summary_query_col_bytes << " bytes (" << MB10(summary_query_col_bytes) << " MB)\n";
    std::cout << "  - Query(local extract)    : " << summary_query_local_bytes << " bytes (" << MB10(summary_query_local_bytes) << " MB)\n";
    std::cout << "  - Query(total)            : " << summary_query_total_bytes << " bytes (" << MB10(summary_query_total_bytes) << " MB)\n";
    std::cout << "  - EvalKey size (coef)     : " << summary_evalkey_bytes << " bytes (" << MB10(summary_evalkey_bytes) << " MB)"
              << "  sent by the client\n";
    std::cout << "  - AutoKS FFT view         : " << summary_aks_fft_bytes << " bytes (" << MB10(summary_aks_fft_bytes) << " MB)"
              << "  derived by the server, not sent\n";
    std::cout << "  - EvalKey stored total    : " << (summary_evalkey_bytes + summary_aks_fft_bytes)
              << " bytes (" << MB10(summary_evalkey_bytes + summary_aks_fft_bytes) << " MB)  server storage\n";

    std::cout << "\n[ONLINE / COMM]\n";
    std::cout << "  - Response before RS      : " << summary_resp_before_bytes << " bytes (" << MB10(summary_resp_before_bytes) << " MB)\n";
    std::cout << "  - Response after  RS      : " << summary_resp_after_bytes  << " bytes (" << MB10(summary_resp_after_bytes)  << " MB)\n";

    std::cout << "\n[OFFLINE / COMP]\n";
    std::cout << "  - preprocessing_ms        : " << summary_T_preprocess_ms << "\n";
    std::cout << "  - aks_keygen_ms           : " << summary_T_aks_ms << "\n";
    std::cout << "  - aks_fft_precompute_ms   : " << summary_T_aks_fft_ms << "\n";
    std::cout << "  - rswk_keygen_ms          : " << summary_T_rswk_gen_ms << "\n";
    std::cout << "  - convk_keygen_ms         : " << summary_T_convk_gen_ms << "\n";
    std::cout << "  - offline_comp_ms         : " << summary_T_offline_comp_ms << "\n";

    std::cout << "\n[ONLINE / COMP]\n";
    std::cout << "  - query_encrypt_ms        : " << summary_T_query_enc_ms << "\n";
    std::cout << "  - col_trace_ms            : " << summary_T_col_trace_ms << "\n";
    std::cout << "  - row_query_unpack_ms     : " << summary_T_row_query_unpack_ms << "\n";
    std::cout << "  - col_query_unpack_ms     : " << summary_T_col_query_unpack_ms << "\n";
    std::cout << "  - query_unpack_ms         : " << summary_T_query_unpack_ms << "\n";
    std::cout << "  - first_dim_ms            : " << summary_T_first_dim_ms << "\n";
    std::cout << "  - second_dim_ms           : " << summary_T_second_dim_ms << "\n";
    std::cout << "  - slot_query_prepare_ms   : " << summary_T_slot_query_prepare_ms << "\n";
    std::cout << "  - slot_query_unpack_ms    : " << summary_T_local_unpack_ms << "\n";
    std::cout << "  - blind_rotation_ms       : " << summary_T_blind_rot_ms << "\n";
    std::cout << "  - local_extract_ms        : " << summary_T_local_ext_ms << "\n";
    std::cout << "  - ring_switch_ms          : " << summary_T_ring_switch_ms << "\n";
    std::cout << "  - main_comp_ms            : " << summary_T_main_comp_ms << "\n";
    std::cout << "  - extract_time_ms         : " << summary_T_extract_time_ms << "\n";
    std::cout << "  - total_ms                : " << summary_T_total_ms << "\n";

    std::cout << "\n[NB]\n";
    std::cout << "  - after_row_expand        : " << summary_nb_after_row_expand << "\n";
    std::cout << "  - after_first_dim         : " << summary_nb_after_first_dim << "\n";
    std::cout << "  - after_col_unpack        : unavailable (TGSW selector stage)\n";
    std::cout << "  - after_second_dim        : " << summary_nb_after_second_dim << "\n";
    std::cout << "  - after_local_extract     : " << summary_nb_after_local_ext << "\n";
    std::cout << "  - before_ring_switch      : " << summary_nb_before_ring_switch << "\n";
    std::cout << "  - final                   : " << summary_nb_final << "\n";

    std::cout << "\n[OUTPUT]\n";
    std::cout << "  - output noise budget     : " << summary_output_noise_budget << "\n";
    std::cout << "  - target index/block/slot : " << index << " / " << target_b << " / " << target_s << "\n";
    print_record_bytes_from_torus("  - expected record bytes   : ", expected_small, N0,
                                  N0 / (int)N1, MM1,
                                  std::min<uint64_t>(16, record_bytes_assumed));
    print_record_bytes_from_torus("  - decrypted record bytes  : ", phase_a2, N0,
                                  N0 / (int)N1, MM1,
                                  std::min<uint64_t>(16, record_bytes_assumed));
    std::cout << "  - record match?           : " << (summary_record_match ? "Yes" : "No")
              << " (bad_coeffs=" << summary_record_bad_coeffs << "/" << N1;
    if (!summary_record_match) {
        std::cout << ", first_bad_coeff=" << summary_record_first_bad
                  << ", expected=" << summary_record_first_exp
                  << ", got=" << summary_record_first_got;
    }
    std::cout << ")\n";
    std::cout << "  - record match (sent b)?  : " << (sent_bad_coeffs == 0 ? "Yes" : "No")
              << " (b at " << N1 << " of " << N0 << " coefficients, bad_coeffs=" << sent_bad_coeffs << ")\n";
    std::cout << "  - decryption failure?     : " << (summary_decryption_failure ? "Yes" : "No") << "\n";
    std::cout << "=================================================\n";

    DeleteAutoKsKeyFFTAll(aks_fft_list, env, n_trace);
    DeleteAutoKsKeyAll(aks_list, env, n_trace);

    const auto t_end_no_pre = t_online_end;
    if (profile) {
        print_ms("TOTAL PROGRAM (with preprocess)", t_program_start, t_end_no_pre);
        print_ms("TOTAL PROGRAM (no preprocess)", t_online_start, t_end_no_pre);
        print_ms("Setup total (synthetic DB)", t_program_start, t_online_start);
        if (summary_T_total_no_pre_ms > 0.0) {
            const double online_sec = summary_T_total_no_pre_ms / 1000.0;
            const double throughput_MBps =
                ((double)summary_n * 256.0 / online_sec) / 1e6;
            printf("[THROUGHPUT] online_total = %.3f MB/s\n", throughput_MBps);
        }
    }

return 0;
}
