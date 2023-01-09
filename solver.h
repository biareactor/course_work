#ifndef SOLVER_H
#define SOLVER_H

#include <cstdlib>
#include <vector>
#include <functional>
#include <utility>
#include <string>
#include <random>

#include <QString>

class Solver
{
    using vec = std::vector<double>;
    using vecvec = std::vector<std::vector<double>>;

public:
    struct Params
    {
        double x0;
        double u1;
        double u2;
        double h0;
        double border;
        double eps_border;
        size_t n;
        double a;
        double alpha;
        double A;
        double omega;
        bool step_control;
        QString method;
        double D;

        double delta;
        double eps;
        size_t p;
    };

    struct Result
    {
        vec x;
        vecvec v; // vector{v1, v2}
        vec step;

        Result() { v.resize(2); }

        vec get_values(size_t idx) const
        {
            vec res{x[idx]};
            for (size_t i = 0; i < v.size(); ++i)
                res.push_back(v[i][idx]);
            res.push_back(step[idx]);

            return res;
        }

        QString get_fields_csv() const
        {
            std::string s = "x,v1,v2,step";

            return s.c_str();
        }

        QString get_values_csv(size_t idx) const
        {
            std::string res;
            for (const auto& value : get_values(idx))
                res += std::to_string(value) + ",";
            res.pop_back();

            return res.c_str();
        }
    };

    Solver(const Params&p): params(p) {}

    template <typename func1, typename func2>
    Result eq_solve(func1, func2);

private:
    Params params;

private:
    template <typename func1, typename func2>
    std::pair<double, double> Euler(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    template <typename func1, typename func2>
    std::pair<double, double> Heun(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    template <typename func1, typename func2>
    std::pair<double, double> RK4(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    template <typename func1, typename func2>
    std::pair<double, double> RK4_double(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    template <typename func1, typename func2>
    std::pair<double, double> Heun_double(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    template <typename func1, typename func2>
    std::pair<double, double> get_v(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    template <typename func1, typename func2>
    std::pair<double, double> get_v2(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h);

    double get_new_step(double s, double step, double x, bool& last_step);
    double get_s(const std::pair<double, double>& v, const std::pair<double, double>& v2);

    void make_step(Solver::Result& res, double x, const std::pair<double, double>& v, double step);
};

template <typename func1, typename func2>
Solver::Result Solver::eq_solve(func1 f1, func2 f2)
{
    Result res;
    if (params.x0 >= params.border - params.eps_border)
        return res;

    make_step(res, params.x0, {params.u1, params.u2}, 0);

    double x = params.x0;
    double step = params.h0;

    size_t i = 0;

    std::pair<double, double> v = {0,0};
    std::pair<double, double> v2 = {0,0};
    double s = 0;
    bool last_step = false;

    while(i < params.n)
    {
        v = get_v(x, v, f1, f2, step);

        if (params.step_control)
        {
            v2 = get_v2(x, v, f1, f2, step);
            s = get_s(v, v2);
            step = get_new_step(s, step, x, last_step);
        }

        if ((x + step) > params.border)
        {
            step = params.border - x - params.eps_border / 2.0;
            last_step = true;
        }

        if (!params.step_control || s < params.eps + params.delta)
        {
            make_step(res, x, v, step);
            x += step;
            ++i;

            if (last_step)
                break;
        }
    }

    return res;
}

template <typename func1, typename func2>
std::pair<double, double> Solver::Euler(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    std::pair<double, double> new_v{0, 0};

    std::random_device rd;
    std::mt19937 RNG(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    double r = sqrt(2.0 * params.D * h);

    new_v.first = v.first + r * dist(RNG) + h * f1(xn, v.first, v.second);
    new_v.second = v.second + r * dist(RNG) + h * f2(xn, v.first, v.second);

    return new_v;
}

template <typename func1, typename func2>
std::pair<double, double> Solver::Heun(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    std::pair<double, double> new_v{0, 0};

    std::random_device rd;
    std::mt19937 RNG(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    double r = sqrt(2.0 * params.D * h);

    double k1_1 = f1(xn, v.first, v.second);
    double k1_2 = f2(xn, v.first, v.second);

    double k2_1 = f1(xn + h, v.first + h * k1_1 + r * dist(RNG), v.second + h * k1_2 + r * dist(RNG));
    double k2_2 = f2(xn + h, v.first + h * k1_1 + r * dist(RNG), v.second + h * k1_2 + r * dist(RNG));

    new_v.first = v.first + h * (k1_1 + k2_1) / 2.0 + r * dist(RNG);
    new_v.second = v.second + h * (k1_2 + k2_2) / 2.0 + r * dist(RNG);

    return new_v;
}

template <typename func1, typename func2>
std::pair<double, double> Solver::RK4(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    std::pair<double, double> new_v{0, 0};

    double k1, k2, k3, k4;

    double k1_1 = f1(xn, v.first, v.second);
    double k1_2 = f2(xn, v.first, v.second);

    double k2_1 = f1(xn + h / 2.0, v.first + h / 2.0 * k1_1, v.second + h / 2.0 * k1_2);
    double k2_2 = f2(xn + h / 2.0, v.first + h / 2.0 * k1_1, v.second + h / 2.0 * k1_2);

    double k3_1 = f1(xn + h / 2.0, v.first + h / 2.0 * k2_1, v.second + h / 2.0 * k2_2);
    double k3_2 = f2(xn + h / 2.0, v.first + h / 2.0 * k2_1, v.second + h / 2.0 * k2_2);

    double k4_1 = f1(xn + h, v.first + h * k3_1, v.second + h * k3_2);
    double k4_2 = f2(xn + h, v.first + h * k3_1, v.second + h * k3_2);

    new_v.first = v.first + h / 6.0 * (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1);
    new_v.second = v.second + h / 6.0 * (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2);

    return new_v;
}

template <typename func1, typename func2>
std::pair<double, double> Solver::RK4_double(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    auto v_half = RK4(xn, v, f1, f2, h/2.0);
    return RK4(xn+h/2.0, v_half, f1, f2, h/2.0);
}

template <typename func1, typename func2>
std::pair<double, double> Solver::Heun_double(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    auto v_half = Heun(xn, v, f1, f2, h/2.0);
    return Heun(xn+h/2.0, v_half, f1, f2, h/2.0);
}

template <typename func1, typename func2>
std::pair<double, double> Solver::get_v(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    std::pair<double, double> v_{0, 0};

    if (params.method == "Хюна")
        v_ = Heun(xn, v, f1, f2, h);
    else if (params.method == "РК 4")
        v_ = RK4(xn, v, f1, f2, h);

    return v_;
}

template <typename func1, typename func2>
std::pair<double, double> Solver::get_v2(double xn, const std::pair<double, double>& v, func1 f1, func2 f2, double h)
{
    std::pair<double, double> v2{0, 0};

    if (params.method == "Хюна")
        v2 = Heun_double(xn, v, f1, f2, h);
    else if (params.method == "РК 4")
        v2 = RK4_double(xn, v, f1, f2, h);

    return v2;
}

#endif // SOLVER_H
