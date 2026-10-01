// Closed wormlike-chain Monte Carlo for the solution conformation of a supercoiled or relaxed
// plasmid. The chain is a closed polygon of n equal segments; its energy is
//
//   E/kT = sum_k (P/l0)(1 - cos theta_k)            bending, persistence length P
//        + (2 pi^2 C / L)(dLk - Wr)^2                twist, C = torsional rigidity / kT
//
// with a hard-core excluded volume of diameter d_eff between segments farther apart than
// the core along the contour. dLk is fixed by the superhelical density, Wr is the writhe of the
// polygon (Klenin and Langowski 2000). This is the model of Vologodskii et al. 1992 and
// Hammermann et al. 1998; the moves are crankshaft rotations of a run of vertices about the
// line through the two fixed vertices that bound it.
//
// Build:  c++ -O2 -std=c++17 -o plasmid_mc plasmid_mc.cpp
// Output: one line "x y z" (nm) per vertex to --out, and a one-line statistics record per
//         sampling interval on stdout.
//
// This file is not part of the TOPAS build, which compiles *.cc only.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

struct V3 { double x, y, z; };
static inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline V3 operator*(double s, V3 a) { return {s * a.x, s * a.y, s * a.z}; }
static inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline V3 cross(V3 a, V3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
static inline double norm(V3 a) { return std::sqrt(dot(a, a)); }
static inline V3 unit(V3 a) { double n = norm(a); return n < 1e-14 ? V3{0, 0, 0} : (1.0 / n) * a; }

static const double PI = 3.14159265358979323846;

// Gauss linking integral of two straight segments (p1,p2) and (p3,p4), as a signed solid-angle
// sum; half of this over 2 pi summed on the pair list gives the writhe (see writhe_pair).
static double omega(V3 p1, V3 p2, V3 p3, V3 p4)
{
    V3 r12 = p2 - p1, r34 = p4 - p3, r13 = p3 - p1, r14 = p4 - p1, r23 = p3 - p2, r24 = p4 - p2;
    V3 n1 = unit(cross(r13, r14)), n2 = unit(cross(r14, r24));
    V3 n3 = unit(cross(r24, r23)), n4 = unit(cross(r23, r13));
    auto as = [](double v) { return std::asin(std::max(-1.0, std::min(1.0, v))); };
    double s = as(dot(n1, n2)) + as(dot(n2, n3)) + as(dot(n3, n4)) + as(dot(n4, n1));
    double sg = dot(cross(r34, r12), r13);
    return sg > 0 ? s : -s;
}

// Minimum distance between segments (a0,a1) and (b0,b1).
static double seg_dist(V3 a0, V3 a1, V3 b0, V3 b1)
{
    V3 d1 = a1 - a0, d2 = b1 - b0, r = a0 - b0;
    double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    double s, t;
    double c = dot(d1, r);
    double b = dot(d1, d2), den = a * e - b * b;
    s = den > 1e-12 ? std::max(0.0, std::min(1.0, (b * f - c * e) / den)) : 0.0;
    t = (b * s + f) / e;
    if (t < 0) { t = 0; s = std::max(0.0, std::min(1.0, -c / a)); }
    else if (t > 1) { t = 1; s = std::max(0.0, std::min(1.0, (b - c) / a)); }
    V3 cp = (a0 + s * d1) - (b0 + t * d2);
    return norm(cp);
}

struct Chain
{
    int n;
    double l0, P, Ctw, L, dLk, deff, h;
    int kmin;                      // segments closer than this along the contour skip the core
    std::vector<V3> q;             // vertices
    std::vector<double> W;         // pair writhe contributions Omega/(2 pi), a<b, n*n
    double wr;

    int cyc(int a, int b) const { int d = std::abs(a - b); return std::min(d, n - d); }
    V3 v(int k) const { return q[((k % n) + n) % n]; }

    double pair_w(int a, int b) const
    {
        return omega(q[a], q[(a + 1) % n], q[b], q[(b + 1) % n]) / (2 * PI);
    }
    double total_writhe()
    {
        W.assign((size_t)n * n, 0.0);
        double s = 0;
        for (int a = 0; a < n; ++a)
            for (int b = a + 2; b < n; ++b)
            {
                if (a == 0 && b == n - 1) continue;
                double w = pair_w(a, b);
                W[(size_t)a * n + b] = w;
                s += w;
            }
        return s;
    }
    double bend_at(int k) const
    {
        V3 u1 = unit(v(k) - v(k - 1)), u2 = unit(v(k + 1) - v(k));
        return (P / l0) * (1.0 - dot(u1, u2));
    }
    double bend_total() const
    {
        double e = 0;
        for (int k = 0; k < n; ++k) e += bend_at(k);
        return e;
    }
    double twist_e(double w) const { double d = dLk - w; return 2 * PI * PI * Ctw / L * d * d; }
    bool core_ok_full() const
    {
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b)
            {
                if (cyc(a, b) < kmin) continue;
                if (seg_dist(q[a], q[(a + 1) % n], q[b], q[(b + 1) % n]) < deff) return false;
            }
        return true;
    }
    double rg() const
    {
        V3 c{0, 0, 0};
        for (auto& p : q) c = c + p;
        c = (1.0 / n) * c;
        double s = 0;
        for (auto& p : q) s += dot(p - c, p - c);
        return std::sqrt(s / n);
    }
};

static void rotate_about(V3 o, V3 ax, double th, V3& p)
{
    V3 r = p - o;
    double c = std::cos(th), s = std::sin(th);
    V3 k = ax;
    p = o + (c * r + s * cross(k, r) + ((1 - c) * dot(k, r)) * k);
}

// Right-handed interwound template (the form negatively supercoiled DNA takes, Wr < 0): two strands of one superhelix joined by a half-circle loop
// at each end, arc-length resampled to n equal segments.
static void plectoneme_template(Chain& ch, double r, double alpha_deg)
{
    double P = 2 * PI * r * std::tan(alpha_deg * PI / 180.0);   // pitch
    double k = std::sqrt(1.0 + std::pow(2 * PI * r / P, 2));
    double H = (ch.L - 2 * PI * r) / (2 * k);
    int M = 40000;
    std::vector<V3> d;
    // strand A rises, right-handed: (r cos t, r sin t, z)
    for (int i = 0; i < M; ++i)
    {
        double z = H * i / M, t = 2 * PI * z / P;
        d.push_back({r * std::cos(t), r * std::sin(t), z});
    }
    double tA = 2 * PI * H / P;
    V3 a_end{r * std::cos(tA), r * std::sin(tA), H};
    // top half-circle loop from A(H) to the opposite point, in the plane of the axis and A(H)
    for (int i = 0; i < M / 20; ++i)
    {
        double f = PI * (double)i / (M / 20);
        V3 cdir = unit(V3{a_end.x, a_end.y, 0});
        d.push_back({r * std::cos(f) * cdir.x, r * std::cos(f) * cdir.y, H + r * std::sin(f)});
    }
    // strand B descends: opposite phase, same handedness
    for (int i = 0; i < M; ++i)
    {
        double z = H - H * i / M, t = 2 * PI * z / P + PI;
        d.push_back({r * std::cos(t), r * std::sin(t), z});
    }
    double tB = PI;
    V3 b_end{r * std::cos(tB), r * std::sin(tB), 0.0};
    for (int i = 0; i < M / 20; ++i)
    {
        double f = PI * (double)i / (M / 20);
        V3 cdir = unit(V3{b_end.x, b_end.y, 0});
        d.push_back({r * std::cos(f) * cdir.x, r * std::cos(f) * cdir.y, -r * std::sin(f)});
    }
    // resample closed polyline at equal arc length
    int m = (int)d.size();
    std::vector<double> s(m + 1, 0.0);
    for (int i = 0; i < m; ++i) s[i + 1] = s[i] + norm(d[(i + 1) % m] - d[i]);
    double Lt = s[m];
    ch.q.clear();
    int j = 0;
    for (int i = 0; i < ch.n; ++i)
    {
        double target = Lt * i / ch.n;
        while (j < m - 1 && s[j + 1] < target) ++j;
        double f = (target - s[j]) / std::max(1e-12, s[j + 1] - s[j]);
        ch.q.push_back(d[j] + f * (d[(j + 1) % m] - d[j]));
    }
    // equalize: rescale so the polygon perimeter is exactly L
    double per = 0;
    for (int i = 0; i < ch.n; ++i) per += norm(ch.q[(i + 1) % ch.n] - ch.q[i]);
    V3 c{0, 0, 0};
    for (auto& p : ch.q) c = c + p;
    c = (1.0 / ch.n) * c;
    for (auto& p : ch.q) p = (ch.L / per) * (p - c);
}

static void circle_template(Chain& ch)
{
    double R = ch.L / (2 * PI);
    ch.q.clear();
    for (int i = 0; i < ch.n; ++i)
    {
        double t = 2 * PI * i / ch.n;
        ch.q.push_back({R * std::cos(t), R * std::sin(t), 0});
    }
    // polygon perimeter is 2 n R sin(pi/n); rescale to L
    double per = 0;
    for (int i = 0; i < ch.n; ++i) per += norm(ch.q[(i + 1) % ch.n] - ch.q[i]);
    for (auto& p : ch.q) p = (ch.L / per) * p;
}

int main(int argc, char** argv)
{
    int bp = 4363;
    double sigma = -0.06, h = 10.5, segnm = 5.0, P = 50.0, Ctw = 49.0, deff = 11.5;
    double r0 = 7.0, alpha0 = 55.0;
    long long moves = 4000000, sample_every = 200000;
    int seed = 1;
    std::string start = "plectoneme", out = "conformation.xyz";
    double rise = 0.34, amax = 1.5;
    for (int i = 1; i < argc; ++i)
    {
        auto is = [&](const char* k) { return !std::strcmp(argv[i], k) && i + 1 < argc; };
        if (is("--bp")) bp = std::atoi(argv[++i]);
        else if (is("--sigma")) sigma = std::atof(argv[++i]);
        else if (is("--helical-repeat")) h = std::atof(argv[++i]);
        else if (is("--segment")) segnm = std::atof(argv[++i]);
        else if (is("--persistence")) P = std::atof(argv[++i]);
        else if (is("--torsion")) Ctw = std::atof(argv[++i]);
        else if (is("--core-diameter")) deff = std::atof(argv[++i]);
        else if (is("--start")) start = argv[++i];
        else if (is("--start-radius")) r0 = std::atof(argv[++i]);
        else if (is("--start-angle")) alpha0 = std::atof(argv[++i]);
        else if (is("--moves")) moves = std::atoll(argv[++i]);
        else if (is("--sample-every")) sample_every = std::atoll(argv[++i]);
        else if (is("--max-angle")) amax = std::atof(argv[++i]);
        else if (is("--seed")) seed = std::atoi(argv[++i]);
        else if (is("--out")) out = argv[++i];
        else { std::fprintf(stderr, "unknown or incomplete option %s\n", argv[i]); return 2; }
    }

    Chain ch;
    ch.L = rise * bp;
    ch.n = (int)std::lround(ch.L / segnm);
    ch.l0 = ch.L / ch.n;
    ch.P = P; ch.Ctw = Ctw; ch.deff = deff; ch.h = h;
    ch.dLk = sigma * (double)bp / h;
    ch.kmin = (int)std::ceil(deff / ch.l0) + 1;
    if (start == "circle") circle_template(ch);
    else plectoneme_template(ch, r0, alpha0);

    if (!ch.core_ok_full())
    {
        std::fprintf(stderr, "starting conformation violates the hard core; adjust --start-radius\n");
        return 3;
    }
    ch.wr = ch.total_writhe();
    double eb = ch.bend_total();
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::printf("# n=%d l0=%.3f nm L=%.1f nm dLk=%.3f kmin=%d start Wr=%.3f Rg=%.1f\n",
                ch.n, ch.l0, ch.L, ch.dLk, ch.kmin, ch.wr, ch.rg());

    long long acc = 0, tried = 0;
    std::vector<V3> saved;
    struct Upd { size_t idx; double val; };
    std::vector<Upd> upd;
    for (long long it = 1; it <= moves; ++it)
    {
        int i = (int)(U(rng) * ch.n);
        // run length m, log-uniform between 1 and n/2 - 2
        double lm = U(rng) * std::log((double)(ch.n / 2 - 2));
        int m = std::max(1, (int)std::exp(lm));
        int j = i + m + 1;
        V3 o = ch.v(i), ax = unit(ch.v(j) - ch.v(i));
        double th = (2 * U(rng) - 1) * amax / std::sqrt((double)m);
        tried++;
        saved.assign(m, V3{0, 0, 0});
        for (int k = 0; k < m; ++k)
        {
            int idx = (i + 1 + k) % ch.n;
            saved[k] = ch.q[idx];
            rotate_about(o, ax, th, ch.q[idx]);
        }
        // bending change at the two bounding vertices
        double de = 0;
        auto ang = [&](int k) { return ch.bend_at(((k % ch.n) + ch.n) % ch.n); };
        double eb_new = ang(i) + ang(j);
        // old values: compute by restoring temporarily
        std::vector<V3> cur(m);
        for (int k = 0; k < m; ++k) { int idx = (i + 1 + k) % ch.n; cur[k] = ch.q[idx]; ch.q[idx] = saved[k]; }
        double eb_old = ang(i) + ang(j);
        for (int k = 0; k < m; ++k) { int idx = (i + 1 + k) % ch.n; ch.q[idx] = cur[k]; }
        de += eb_new - eb_old;

        // affected segments: i .. j-1 (cyclic); rigid interior: i+1 .. j-2
        bool ok = true;
        auto inA = [&](int s) { int d = ((s - i) % ch.n + ch.n) % ch.n; return d <= m; };           // i .. j-1
        auto inR = [&](int s) { int d = ((s - i) % ch.n + ch.n) % ch.n; return d >= 1 && d <= m - 1; };
        upd.clear();
        double dw = 0;
        for (int da = 0; da <= m && ok; ++da)
        {
            int a = (i + da) % ch.n;
            for (int b = 0; b < ch.n; ++b)
            {
                int c = ch.cyc(a, b);
                if (c <= 1) continue;
                bool bA = inA(b);
                if (bA && b < a) continue;
                if (bA && inR(a) && inR(b)) continue;
                if (c >= ch.kmin &&
                    seg_dist(ch.q[a], ch.q[(a + 1) % ch.n], ch.q[b], ch.q[(b + 1) % ch.n]) < ch.deff)
                { ok = false; break; }
            }
        }
        if (ok)
        {
            for (int da = 0; da <= m; ++da)
            {
                int a = (i + da) % ch.n;
                for (int b = 0; b < ch.n; ++b)
                {
                    int c = ch.cyc(a, b);
                    if (c <= 1) continue;
                    bool bA = inA(b);
                    if (bA && b < a) continue;
                    if (bA && inR(a) && inR(b)) continue;
                    int lo = std::min(a, b), hi = std::max(a, b);
                    double w = ch.pair_w(lo, hi);
                    size_t idx = (size_t)lo * ch.n + hi;
                    dw += w - ch.W[idx];
                    upd.push_back({idx, w});
                }
            }
            double wr_new = ch.wr + dw;
            de += ch.twist_e(wr_new) - ch.twist_e(ch.wr);
            if (de <= 0 || U(rng) < std::exp(-de))
            {
                for (auto& u : upd) ch.W[u.idx] = u.val;
                ch.wr = wr_new;
                acc++;
            }
            else ok = false;
        }
        if (!ok)
            for (int k = 0; k < m; ++k) ch.q[(i + 1 + k) % ch.n] = saved[k];

        if (it % sample_every == 0)
        {
            double wcheck = ch.total_writhe();       // also refreshes the cache against drift
            ch.wr = wcheck;
            std::printf("%lld acc=%.3f Wr=%.3f dLk=%.2f Rg=%.2f E_bend=%.1f E_tw=%.2f\n", it,
                        (double)acc / tried, ch.wr, ch.dLk, ch.rg(), ch.bend_total(),
                        ch.twist_e(ch.wr));
            std::fflush(stdout);
            acc = tried = 0;
        }
    }
    if (!ch.core_ok_full()) { std::fprintf(stderr, "hard core violated at end: bug\n"); return 4; }
    FILE* f = std::fopen(out.c_str(), "w");
    for (auto& p : ch.q) std::fprintf(f, "%.5f %.5f %.5f\n", p.x, p.y, p.z);
    std::fclose(f);
    return 0;
}
