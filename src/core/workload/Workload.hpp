#pragma once
#include <memory>
#include <string>
#include <vector>

namespace mf {

/// Computational activity u(t) in [0, 1] for one server. The CNN inference workload will be
/// another implementation of this interface (requests -> activity) in a later milestone.
class Workload {
public:
    virtual ~Workload() = default;
    virtual double u(double t) const = 0;
    virtual std::string describe() const = 0;
};

class ConstantWorkload final : public Workload {
public:
    explicit ConstantWorkload(double u) : u_(u) {}
    double u(double) const override { return u_; }
    void set(double u) { u_ = u; }
    std::string describe() const override;
private:
    double u_;
};

class StepWorkload final : public Workload {
public:
    StepWorkload(double before, double after, double tStep) : a_(before), b_(after), t_(tStep) {}
    double u(double t) const override { return t < t_ ? a_ : b_; }
    std::string describe() const override;
private:
    double a_, b_, t_;
};

class SquareWorkload final : public Workload {
public:
    SquareWorkload(double low, double high, double period, double duty, double phase = 0.0)
        : lo_(low), hi_(high), period_(period), duty_(duty), phase_(phase) {}
    double u(double t) const override;
    std::string describe() const override;
private:
    double lo_, hi_, period_, duty_, phase_;
};

/// Piecewise-linear trace read from a CSV with columns "t,u" (header line optional).
class TraceWorkload final : public Workload {
public:
    explicit TraceWorkload(const std::string& csvPath);
    double u(double t) const override;
    std::string describe() const override { return "trace " + path_; }
private:
    std::string path_;
    std::vector<double> t_, u_;
};

} // namespace mf
