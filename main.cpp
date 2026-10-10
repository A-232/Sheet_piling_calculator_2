// ============================================================
//  Шпунтовая стенка — предельное равновесие с перебором h
//  Схемы: console / strut / frame / ring
//
//  Сборка:
//    g++ -O2 -static -static-libgcc -static-libstdc++ -std=c++17 main.cpp -o main.exe
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;

const double PI = 3.14159265358979;

string trim(string s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
string toLower(string s) {
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}
map<string,string> readKV(const string& fname) {
    ifstream in(fname);
    map<string,string> kv;
    string line;
    while (getline(in, line)) {
        size_t c = line.find('#');
        if (c != string::npos) line = line.substr(0, c);
        size_t eq = line.find('=');
        if (eq == string::npos) continue;
        kv[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }
    return kv;
}
double num(const map<string,string>& kv, const string& key, double def) {
    auto it = kv.find(key);
    return it == kv.end() ? def : stod(it->second);
}
string strk(const map<string,string>& kv, const string& key, string def) {
    auto it = kv.find(key);
    return it == kv.end() ? def : it->second;
}

struct Params {
    string scheme;
    double h, support_z;
    double gamma, phi, q, c;
    double m, gamma_c;
    double E, I, W, Ry, y_lim;
    double h_start, h_step, h_max;
    double strut_s, strut_A, strut_Rc;
    double ring_R, ring_A;
    double frame_B, frame_W, frame_Ry;
    double N_gamma, N_c, N_q;
};
struct Coef { double A, B; };

double Ka(double phi_deg) {
    double phi = phi_deg * PI / 180.0;
    return pow(tan(PI/4 - phi/2), 2);
}
double Kp(double phi_deg) {
    double phi = phi_deg * PI / 180.0;
    return pow(tan(PI/4 + phi/2), 2);
}

// Активное: σa(z) = γ·z·Ka + q·Ka − 2c·√Ka  (z от верха стенки)
Coef activeCoef(const Params& P, double Ka_v) {
    Coef c;
    c.A = P.gamma * Ka_v;
    c.B = P.q * Ka_v - 2.0 * P.c * sqrt(Ka_v);
    return c;
}
// Пассивное: σp(z') = γ·z'·Kp + 2c·√Kp  (z' от дна котлована, БЕЗ γ·H·Kp!)
Coef passiveCoef(const Params& P, double Kp_v) {
    Coef c;
    c.A = P.gamma * Kp_v;
    c.B = 2.0 * P.c * sqrt(Kp_v);
    return c;
}
void forceCentroid(const Coef& c, double L, double& F, double& zc) {
    F = c.A * L * L / 2.0 + c.B * L;
    double M0 = c.A * L * L * L / 3.0 + c.B * L * L / 2.0;
    zc = (fabs(F) > 1e-12) ? M0 / F : 0.0;
}

// Mmax численно по всей длине стенки.
// Na — усилие распорки (>0 — вклад учитывается; точка приложения — support_z)
double MmaxNumerical(double h, const Params& P, const Coef& ca, const Coef& cp,
                     double Na, double& z_out)
{
    double L = P.h + h;
    double Mmax = 0.0; z_out = 0.0;
    int N = 1000;
    for (int i = 0; i <= N; ++i) {
        double z = L * i / (double)N;
        double Ma = ca.A * z*z*z / 6.0 + ca.B * z*z / 2.0;
        double M = -Ma;
        if (Na > 0 && z > P.support_z)
            M += Na * (z - P.support_z);
        if (z > P.h) {
            double zp = z - P.h;
            double Mp = cp.A * zp*zp*zp / 6.0 + cp.B * zp*zp / 2.0;
            M += Mp;
        }
        if (fabs(M) > fabs(Mmax)) { Mmax = M; z_out = z; }
    }
    return Mmax;
}

// Проверка грунта: (Kp−Ka)·γ·z' ≤ p_пред
struct SoilCheck { double maxRatio, z_crit, p_act, p_lim; bool ok; };
SoilCheck checkSoil(const Params& P, const Coef& cp, double h, double Ka_v) {
    SoilCheck r; r.maxRatio = 0; r.z_crit = 0; r.p_act = 0; r.p_lim = 0; r.ok = true;
    int N = 100;
    for (int i = 1; i <= N; ++i) {
        double z = h * i / (double)N;
        double p_act = P.gamma * z * (Kp(P.phi) - Ka_v);  // чистое сопротивление
        double p_lim = P.gamma * z * P.N_gamma + P.c * P.N_c + P.q * P.N_q;
        double ratio = (p_lim > 1e-9) ? p_act / p_lim : 999.0;
        if (ratio > r.maxRatio) {
            r.maxRatio = ratio; r.z_crit = z;
            r.p_act = p_act; r.p_lim = p_lim;
        }
    }
    r.ok = (r.maxRatio <= 1.0);
    return r;
}

int main(int argc, char* argv[]) {
    string inName  = argc > 1 ? argv[1] : "parameters.txt";
    string outName = argc > 2 ? argv[2] : "result.txt";

    ifstream fin(inName);
    if (!fin) { cout << "Не найден " << inName << "\n"; return 1; }
    fin.close();
    auto kv = readKV(inName);
    if (kv.empty()) { cout << "Файл " << inName << " пуст.\n"; return 1; }

    Params P;
    P.scheme   = toLower(trim(strk(kv, "scheme", "console")));
    P.h        = num(kv, "h", 3.0);
    P.gamma    = num(kv, "gamma", 18.0);
    P.phi      = num(kv, "phi", 30.0);
    P.q        = num(kv, "q", 0.0);
    P.c        = num(kv, "c", 0.0);
    P.m        = num(kv, "m", 0.95);
    P.gamma_c  = num(kv, "gamma_c", 1.0);
    P.E        = num(kv, "E", 2.06e8);
    P.I        = num(kv, "I", 3.0e-4);
    P.W        = num(kv, "W", 1.67e-3);
    P.Ry       = num(kv, "Ry", 240000.0);
    P.y_lim    = num(kv, "y_lim", 0.03);
    P.h_start  = num(kv, "h_start", 0.10);
    P.h_step   = num(kv, "h_step", 0.10);
    P.h_max    = num(kv, "h_max", 12.0);
    P.support_z= num(kv, "support_z", 0.0);
    P.strut_s  = num(kv, "strut_s", 4.0);
    P.strut_A  = num(kv, "strut_A", 56e-4);
    P.strut_Rc = num(kv, "strut_Rc", 200000.0);
    P.ring_R   = num(kv, "ring_R", 8.0);
    P.ring_A   = num(kv, "ring_A", 61.2e-4);
    P.frame_B  = num(kv, "frame_B", 6.0);
    P.frame_W  = num(kv, "frame_W", 242e-6);
    P.frame_Ry = num(kv, "frame_Ry", 240000.0);
    P.N_gamma  = num(kv, "N_gamma", 15.0);
    P.N_c      = num(kv, "N_c", 30.0);
    P.N_q      = num(kv, "N_q", 18.0);

    if (P.scheme != "console" && P.scheme != "strut" &&
        P.scheme != "frame"   && P.scheme != "ring") {
        cout << "Неизвестная схема '" << P.scheme << "'\n"; return 1;
    }

    double Ka_v = Ka(P.phi), Kp_v = Kp(P.phi);
    Coef ca = activeCoef(P, Ka_v);
    Coef cp = passiveCoef(P, Kp_v);

    ofstream out(outName);
    if (!out) { cout << "Не открыть " << outName << "\n"; return 1; }
    out << fixed << setprecision(4);

    out << "============================================================\n";
    out << "  РАСЧЁТ ШПУНТОВОЙ СТЕНКИ (предельное равновесие)\n";
    out << "============================================================\n";
    out << "Схема: " << P.scheme << "\n";
    out << "H = " << P.h << " м; γ = " << P.gamma << " кН/м³; φ = "
        << P.phi << "°; c = " << P.c << " кПа; q = " << P.q << " кПа\n";
    out << "Ka = " << Ka_v << "; Kp = " << Kp_v
        << "; m = " << P.m << "; γc = " << P.gamma_c << "\n";
    out << "Профиль: I = " << P.I << " м⁴/м; W = " << P.W
        << " м³/м; Ry = " << P.Ry << " кПа\n";
    out << "Перебор h: от " << P.h_start << " до " << P.h_max
        << " м с шагом " << P.h_step << " м\n";
    out << "------------------------------------------------------------\n\n";
    out << ">>> Итерационный подбор h\n\n";

    bool isStrutScheme = (P.scheme != "console");
    bool found = false;
    double h_found = -1;
    double Ea_f=0, za_f=0, Ep_f=0, zp_f=0, Na_f=0;
    double Mmax_f=0, zM_f=0, sigma_f=0, y_f=0, ratio_f=0;
    double N_strut_f=0, N_cap_strut_f=0;
    double M_frame_f=0, M_cap_frame_f=0;
    double N_ring_f=0,  N_cap_ring_f=0;

    int iter = 0;
    for (double h = P.h_start; h <= P.h_max + 1e-9; h += P.h_step) {
        if (++iter > 5000) break;
        double L = P.h + h;

        double Ea, za, Ep, zp;
        forceCentroid(ca, L, Ea, za);
        forceCentroid(cp, h, Ep, zp);

        // --- равновесие моментов ---
        double M_opr, M_ud;
        if (P.scheme == "console") {
            M_opr = Ea * (L - za);
            M_ud  = Ep * (h - zp);
        } else {
            M_opr = Ea * (za - P.support_z);
            M_ud  = Ep * ((P.h + zp) - P.support_z);
        }
        bool eq_ok = (M_opr <= P.m * M_ud + 1e-6);

        // --- усилие распорки ---
        double Na = (P.scheme == "console") ? 0.0 : (Ea - Ep);
        bool na_ok = !isStrutScheme || (Na > 1e-9);   // распорка должна работать на сжатие

        // --- Mmax и прочность ---
        double Mmax, zM;
        Mmax = MmaxNumerical(h, P, ca, cp, Na, zM);
        double sigma = fabs(Mmax) / P.W;
        bool str_ok = (sigma <= P.Ry * P.gamma_c + 1e-6);

        // --- грунт ---
        SoilCheck sc = checkSoil(P, cp, h, Ka_v);
        bool soil_ok = sc.ok;

        // --- прогиб (оценка) ---
        double qa_H = ca.A * P.h + ca.B;
        double q_avg = 0.5 * qa_H;
        double l_span = (P.scheme == "console") ? P.h : (P.h - P.support_z);
        if (l_span <= 0) l_span = P.h;
        double y_max = 5.0 * q_avg * pow(l_span, 4) / (384.0 * P.E * P.I);
        bool defl_ok = (y_max <= P.y_lim + 1e-9);

        // --- распорка / обвязка ---
        bool sup_ok = true;
        double N_strut = 0, N_cap_strut = 0;
        double M_frame = 0, M_cap_frame = 0;
        double N_ring  = 0, N_cap_ring  = 0;
        if (P.scheme == "strut") {
            N_strut = fabs(Na) * P.strut_s;
            N_cap_strut = P.strut_Rc * P.strut_A;
            sup_ok = (N_strut <= N_cap_strut + 1e-6);
        } else if (P.scheme == "frame") {
            M_frame = fabs(Na) * P.frame_B * P.frame_B / 8.0;
            M_cap_frame = P.frame_Ry * P.frame_W;
            sup_ok = (M_frame <= M_cap_frame + 1e-6);
        } else if (P.scheme == "ring") {
            N_ring = fabs(Na) * P.ring_R;
            N_cap_ring = P.Ry * P.ring_A;
            sup_ok = (N_ring <= N_cap_ring + 1e-6);
        }

        // --- лог строки ---
        out << "  h = " << setw(6) << h << " м:  "
            << "равн." << (eq_ok   ? "[+]" : "[-]") << "  "
            << "расп." << (na_ok   ? "[+]" : "[-]") << "  "
            << "прочн."<< (str_ok  ? "[+]" : "[-]") << "  "
            << "грунт" << (soil_ok ? "[+]" : "[-]") << "  "
            << "прогиб"<< (defl_ok ? "[+]" : "[-]") << "  "
            << "подп." << (sup_ok  ? "[+]" : "[-]");
        if (eq_ok && na_ok && str_ok && soil_ok && defl_ok && sup_ok) {
            out << "   <<< ОК\n";
            found = true;
            h_found = h;
            Ea_f=Ea; za_f=za; Ep_f=Ep; zp_f=zp; Na_f=Na;
            Mmax_f=Mmax; zM_f=zM; sigma_f=sigma; y_f=y_max;
            ratio_f = sc.maxRatio;
            N_strut_f=N_strut; N_cap_strut_f=N_cap_strut;
            M_frame_f=M_frame; M_cap_frame_f=M_cap_frame;
            N_ring_f=N_ring;   N_cap_ring_f=N_cap_ring;
            break;
        }
        out << "   (";
        if (!eq_ok)   out << "момент, ";
        if (!na_ok)   out << "распорка в растяжении, ";
        if (!str_ok)  out << "σ>" << P.Ry*P.gamma_c << ", ";
        if (!soil_ok) out << "грунт, ";
        if (!defl_ok) out << "прогиб, ";
        if (!sup_ok)  out << "подпор, ";
        out << ")\n";
    }

    if (!found) {
        out << "\n!!! Решение не найдено в диапазоне до " << P.h_max << " м\n";
        out << "Увеличьте h_max, измените профиль или схему.\n";
        out.close();
        cout << "Решение не найдено. См. " << outName << "\n";
        return 1;
    }

    out << "\n============================================================\n";
    out << "         РЕЗУЛЬТАТ\n";
    out << "============================================================\n";
    out << "Минимальная глубина заделки: h = " << h_found << " м\n";
    out << "Полная длина: L = H + h = " << (P.h + h_found) << " м\n";
    out << "Ea = " << Ea_f << " кН/м  (центр za = " << za_f << " м)\n";
    out << "Ep = " << Ep_f << " кН/м  (центр zp = " << zp_f
        << " м от дна котлована)\n";
    if (P.scheme != "console")
        out << "Na = " << Na_f << " кН/м  (распорка, сжатие)\n";
    out << "Mmax = " << fabs(Mmax_f) << " кН·м при z = " << zM_f << " м\n";
    out << "σ = " << sigma_f << " кПа при Ry·γc = "
        << (P.Ry*P.gamma_c) << " кПа\n";
    out << "y_max ≈ " << (y_f*1000.0) << " мм при [y] = "
        << (P.y_lim*1000.0) << " мм\n";
    out << "Грунт: max(p/p_пред) = " << ratio_f << "\n";
    if (P.scheme == "strut")
        out << "Распорка: N = " << N_strut_f
            << " кН; Rc·A = " << N_cap_strut_f << " кН\n";
    else if (P.scheme == "frame")
        out << "Обвязка: M = " << M_frame_f
            << " кН·м; Ry·W = " << M_cap_frame_f << " кН·м\n";
    else if (P.scheme == "ring")
        out << "Кольцо: N = " << N_ring_f
            << " кН; Ry·A = " << N_cap_ring_f << " кН\n";
    out << "============================================================\n";

    out.close();
    cout << "Готово. Найдено h = " << h_found << " м. См. " << outName << "\n";
    return 0;
}