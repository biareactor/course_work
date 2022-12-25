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
    double k2_2 = f2(xn + h, v1 + h * k1_1 + r * dist(RNG), v2 + h * k1_2 + r * dist(RNG));

    v1 += h * (k1_1 + k2_1) / 2.0 + r * dist(RNG);
    v2 += h * (k1_2 + k2_2) / 2.0 + r * dist(RNG);

    return std::make_pair(v1, v2);
}

static double modulo(double a)
{
    double r = 2 * M_PI;
    return a - floor (a / r) * r;
}

template <typename func1, typename func2>
Solver::Result Solver::eq_solve(func1 f1, func2 f2)
{
    Result res;
    res.x.push_back(params.x0);
    res.v.resize(2);
    res.v[0].push_back(params.u1);
    res.v[1].push_back(params.u2);
    res.step.push_back(0);

    double x = params.x0;
    double step = params.h0;

    size_t i = 0;
    while(i < params.n && x < params.border - params.eps_border)
    {
        while ((x + step) > params.border)
            step = params.border - x - params.eps_border / 2.0;

        x += step;

        double r = sqrt(2.0 * params.D * step);
        auto v = HeunScheme(x, res.v[0].back(), res.v[1].back(), f1, f2, step, r);

        res.x.push_back(x);
        res.v[0].push_back(modulo(v.first));
//        res.v[0].push_back(v.first);
        res.v[1].push_back(v.second);
        res.step.push_back(step);

        ++i;
    }

    return res;
}

#endif // SOLVER_H
