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
    };

    struct Result
    {
        vec x;
        vecvec v; // vector{v1, v2}
        vec step;

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
};

template <typename func1, typename func2>
std::pair<double, double> EulerScheme(double xn, double v1, double v2, func1 f1, func2 f2, double h, double r)
{
    std::random_device rd;
    std::mt19937 RNG(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    v1 += r * dist(RNG) + h * f1(xn, v1, v2);
    v2 += r * dist(RNG) + h * f2(xn, v1, v2);
    return std::make_pair(v1, v2);
}

template <typename func1, typename func2>
std::pair<double, double> HeunScheme(double xn, double v1, double v2, func1 f1, func2 f2, double h, double r)
{
    std::random_device rd;
    std::mt19937 RNG(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    double k1_1 = f1(xn, v1, v2);
    double k1_2 = f2(xn, v1, v2);

    double k2_1 = f1(xn + h, v1 + h * k1_1 + r * dist(RNG), v2 + h * k1_2 + r * dist(RNG));
                //f1(xn + h, v1 + h * k1_1 + r * dist(RNG), v2 + h * k1_1 + r * dist(RNG));
    double k2_2 = f2(xn + h, v1 + h * k1_1 + r * dist(RNG), v2 + h * k1_2 + r * dist(RNG));
                //f2(xn + h, v1 + h * k1_2 + r * dist(RNG), v2 + h * k1_2 + r * dist(RNG));

    v1 += h * (k1_1 + k2_1) / 2.0 + r * dist(RNG);
    v2 += h * (k1_2 + k2_2) / 2.0 + r * dist(RNG);

    return std::make_pair(v1, v2);
}

template <typename func1, typename func2>
std::pair<double, double> RK4(double xn, double v1, double v2, func1 f1, func2 f2, double h)
{
    double k1, k2, k3, k4;

    double k1_1 = f1(xn, v1, v2);
    double k1_2 = f2(xn, v1, v2);

    double k2_1 = //f1(xn + h / 2.0, v1 + h / 2.0 * k1_1, v2 + h / 2.0 * k1_2);
                f1(xn + h / 2.0, v1 + h / 2.0 * k1_1, v2 + h / 2.0 * k1_1);
    double k2_2 = //f2(xn + h / 2.0, v1 + h / 2.0 * k1_1, v2 + h / 2.0 * k1_2);
                f2(xn + h / 2.0, v1 + h / 2.0 * k1_2, v2 + h / 2.0 * k1_2);

    double k3_1 = f1(xn + h / 2.0, v1 + h / 2.0 * k2_1, v2 + h / 2.0 * k2_2);
    double k3_2 = f2(xn + h / 2.0, v1 + h / 2.0 * k2_1, v2 + h / 2.0 * k2_2);

    double k4_1 = f1(xn + h, v1 + h * k3_1, v2 + h * k3_2);
    double k4_2 = f2(xn + h, v1 + h * k3_1, v2 + h * k3_2);

    v1 += h / 6.0 * (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1);
    v2 += h / 6.0 * (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2);

    return std::make_pair(v1, v2);
}

template <typename func1, typename func2>
std::pair<double, double> RK4_double(double xn, double v1, double v2, func1 f1, func2 f2, double h)
{
    auto v = RK4(xn, v1, v2, f1, f2, h/2.0);
    return RK4(xn+h/2.0, v.first, v.second, f1, f2, h/2.0);
}

template <typename func1, typename func2>
std::pair<double, double> HeunScheme_double(double xn, double v1, double v2, func1 f1, func2 f2, double h)
{
    auto v = HeunScheme(xn, v1, v2, f1, f2, h/2.0, 0);
    return HeunScheme(xn+h/2.0, v.first, v.second, f1, f2, h/2.0, 0);
}

static double modulo(double a)
{
    double r = 2 * M_PI;
    return a - floor (a / r) * r;
}

static double delta = std::pow(10, -10);
static double eps = std::pow(10, -3);
static size_t p = 4;
//size_t p = 2;

static void step_control(double s, size_t& inc_step, size_t& dec_step)
{
    double eps_min = eps / pow(2, p + 1);

    if (s < eps_min - delta)
        ++inc_step;
    else if (s > eps + delta)
        ++dec_step;
}

template <typename func1, typename func2>
Solver::Result Solver::eq_solve(func1 f1, func2 f2)
{
    Result res;
    if (params.x0 >= params.border - params.eps_border)
        return res;

    res.x.push_back(params.x0);
    res.v.resize(2);
    res.v[0].push_back(params.u1);
    res.v[1].push_back(params.u2);
    res.step.push_back(0);

    double x = params.x0;
    double step = params.h0;

    size_t i = 0;
    size_t inc_step = 0;
    size_t dec_step = 0;
    bool last_step = false;

    std::pair<double, double> v = {0,0};
    std::pair<double, double> v2 = {0,0};
    double s = 0;

    while(i < params.n || last_step)
    {
        if (params.method == "Хюна")
        {
            double r = sqrt(2.0 * params.D * step);
            v = HeunScheme(x, v.first, v.second, f1, f2, step, r);

            if (params.step_control)
                v2 = HeunScheme_double(x, v.first, v.second, f1, f2, step);
        }
        else if (params.method == "РК 4")
        {
            v = RK4(x, v.first, v.second, f1, f2, step);

            if (params.step_control)
                v2 = RK4_double(x, v.first, v.second, f1, f2, step);
        }

        if (params.step_control)
        {
            double s1 = (v2.first - v.first) / (pow(2, p) - 1);
            double s2 = (v2.second - v.second) / (pow(2, p) - 1);
            s = std::max(std::abs(s1), std::abs(s2));
            //double local_err = pow(2, p)*s;

            step_control(s, inc_step, dec_step);
            step /= std::pow(2, dec_step);
        }

        if (!params.step_control || s < eps + delta)
        {
            res.x.push_back(x);
            res.v[0].push_back(modulo(v.first));
     //        res.v[0].push_back(v.first);
            res.v[1].push_back(v.second);
            res.step.push_back(step);
            x += step;
            ++i;

            if (last_step)
                break;

            if (params.step_control)
                step *= std::pow(2, inc_step);

            if ((x + step) > params.border)
            {
                step = params.border - x - params.eps_border / 2.0;
                last_step = true;
            }

            inc_step = 0;
            dec_step = 0;
        }
    }

    return res;
}

//template <typename func1, typename func2>
//Solver::Result Solver::eq_solve(func1 f1, func2 f2)
//{
//    Result res;
//    res.x.push_back(params.x0);
//    res.v.resize(2);
//    res.v[0].push_back(params.u1);
//    res.v[1].push_back(params.u2);
//    res.step.push_back(0);

//    double x = params.x0;
//    double step = params.h0;

//    size_t i = 0;
//    while(i < params.n && x < params.border - params.eps_border)
//    {
//        while ((x + step) > params.border)
//            step = params.border - x - params.eps_border / 2.0;

//        x += step;

//        double r = sqrt(2.0 * params.D * step);
// //        auto v = HeunScheme(x, res.v[0].back(), res.v[1].back(), f1, f2, step, r);
//        auto v = RK4(x, res.v[0].back(), res.v[1].back(), f1, f2, step);

//        res.x.push_back(x);
//        res.v[0].push_back(modulo(v.first));
// //        res.v[0].push_back(v.first);
//        res.v[1].push_back(v.second);
//        res.step.push_back(step);

//        ++i;
//    }

//    return res;
//}

#endif // SOLVER_H
