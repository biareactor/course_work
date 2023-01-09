#include "solver.h"

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

static double modulo(double a)
{
    double r = 2 * M_PI;
    return a - floor (a / r) * r;
}

void Solver::make_step(Solver::Result& res, double x, const std::pair<double, double>& v, double step)
{
    res.x.push_back(x);
    res.v[0].push_back(modulo(v.first));
    res.v[1].push_back(v.second);
    res.step.push_back(step);
}
