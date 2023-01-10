#include "solver.h"
#include <cmath>

double Solver::get_new_step(double s, double step, double x, bool& last_step)
{
    double eps_min = params.eps / pow(2, params.p + 1);

    if (s < eps_min - params.delta)
        step *= 2;
    else if (s > params.eps + params.delta)
        step /= 2;

    return step;
}

double Solver::get_s(const std::pair<double, double>& v, const std::pair<double, double>& v2)
{
    double s1 = (v2.first - v.first) / (pow(2, params.p) - 1);
    double s2 = (v2.second - v.second) / (pow(2, params.p) - 1);
    return std::max(std::abs(s1), std::abs(s2));
}



static void push_val_to_v2(Solver::Result& res, const std::pair<double, double>& v)
{
    double prev_val = !res.v1.empty() ? res.v1.back() : 0;

    if (res.v2.empty() || (prev_val < (2*M_PI) * std::ceil(prev_val/(2*M_PI)) && v.first > (2*M_PI) * std::ceil(prev_val/(2*M_PI))))
        res.v2.push_back({});

    res.v2.back().push_back(v.second);
}

void Solver::make_step(Solver::Result& res, double x, const std::pair<double, double>& v, double step)
{
    res.x.push_back(x);
    push_val_to_v2(res, v);
    res.v1.push_back(v.first);
    res.step.push_back(step);
}
