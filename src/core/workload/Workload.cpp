#include "workload/Workload.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mf {

std::string ConstantWorkload::describe() const {
    std::ostringstream s;
    s << "constant u=" << u_;
    return s.str();
}

std::string StepWorkload::describe() const {
    std::ostringstream s;
    s << "step u=" << a_ << "->" << b_ << " at t=" << t_ << "s";
    return s.str();
}

double SquareWorkload::u(double t) const {
    double x = std::fmod(t + phase_, period_);
    if (x < 0) x += period_;
    return x < duty_ * period_ ? hi_ : lo_;
}

std::string SquareWorkload::describe() const {
    std::ostringstream s;
    s << "square u=" << lo_ << "/" << hi_ << " period=" << period_ << "s duty=" << duty_;
    return s.str();
}

TraceWorkload::TraceWorkload(const std::string& path) : path_(path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open workload trace " + path);
    std::string line;
    while (std::getline(in, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        double t, u;
        if (ss >> t >> u) { t_.push_back(t); u_.push_back(std::clamp(u, 0.0, 1.0)); }
    }
    if (t_.empty()) throw std::runtime_error("workload trace " + path + " has no data rows");
}

double TraceWorkload::u(double t) const {
    if (t <= t_.front()) return u_.front();
    if (t >= t_.back()) return u_.back();
    const auto it = std::upper_bound(t_.begin(), t_.end(), t);
    const size_t i = static_cast<size_t>(it - t_.begin());
    const double f = (t - t_[i - 1]) / (t_[i] - t_[i - 1]);
    return u_[i - 1] + f * (u_[i] - u_[i - 1]);
}

} // namespace mf
