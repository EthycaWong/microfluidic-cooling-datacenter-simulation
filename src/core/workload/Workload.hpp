#pragma once
#include <chrono>
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

/// One sample of a recorded compute job (e.g. the PyTorch CNN training run).
struct TracePoint {
    double t = 0, u = 0;
    std::string phase;
    double trainSamples = 0, inferSamples = 0;
};

/// Replays a job trace written by ml/activity.py (columns t_s,u,phase,train_samples,infer_samples).
/// Simulated time t maps to trace time (t - startAt) / timeScale: before the job starts and after it
/// ends the servers are idle (u = 0), so the job arrives "on demand". With live = true the file is
/// re-read as it grows, so a running training job drives the simulation in real time.
class JobTraceWorkload : public Workload {
public:
    JobTraceWorkload(std::string path, double startAt = 0.0, double timeScale = 1.0, bool live = false);
    double u(double t) const override;
    std::string describe() const override;
    /// Sample in effect at simulated time t (phase "waiting" before the job, "done" after it).
    TracePoint at(double t) const;
    double durationSim() const;           ///< simulated seconds the recorded job lasts
    double startAt() const { return start_; }
    void setStartAt(double t) { start_ = t; }
    double timeScale() const { return scale_; }
    bool live() const { return live_; }
    size_t size() const { return pts_.size(); }
    const std::vector<TracePoint>& points() const { return pts_; }
    /// Live mode: read any new rows (cheap when nothing changed). Called automatically by u().
    void refresh() const;

private:
    std::string path_;
    double start_, scale_;
    bool live_;
    mutable std::vector<TracePoint> pts_;
    mutable long long offset_ = 0;
    mutable std::string partial_;
    mutable std::chrono::steady_clock::time_point lastRead_{};
};

} // namespace mf
