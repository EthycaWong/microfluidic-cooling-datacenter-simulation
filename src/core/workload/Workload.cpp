#include "workload/Workload.hpp"

#include <algorithm>
#include <chrono>
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

// ---- JobTraceWorkload ------------------------------------------------------------------------------
JobTraceWorkload::JobTraceWorkload(std::string path, double startAt, double timeScale, bool live)
    : path_(std::move(path)), start_(startAt), scale_(timeScale > 0 ? timeScale : 1.0), live_(live) {
    refresh();
    if (!live_ && pts_.empty()) throw std::runtime_error("job trace " + path_ + " has no data rows");
}

void JobTraceWorkload::refresh() const {
    std::ifstream in(path_, std::ios::binary);
    if (!in) return; // live: the job has not started writing yet
    in.seekg(0, std::ios::end);
    const long long size = static_cast<long long>(in.tellg());
    if (size < offset_) { pts_.clear(); offset_ = 0; partial_.clear(); } // file was restarted
    if (size == offset_) return;
    in.seekg(offset_);
    std::string chunk(static_cast<size_t>(size - offset_), '\0');
    in.read(chunk.data(), static_cast<std::streamsize>(chunk.size()));
    offset_ = size;
    partial_ += chunk;
    size_t pos;
    while ((pos = partial_.find('\n')) != std::string::npos) {
        std::string line = partial_.substr(0, pos);
        partial_.erase(0, pos + 1);
        if (line.empty() || line[0] == 't') continue; // header
        std::stringstream ss(line);
        std::string f[5];
        int n = 0;
        while (n < 5 && std::getline(ss, f[n], ',')) ++n;
        if (n < 2) continue;
        try {
            TracePoint p;
            p.t = std::stod(f[0]);
            p.u = std::clamp(std::stod(f[1]), 0.0, 1.0);
            if (n > 2) p.phase = f[2];
            if (n > 3) p.trainSamples = std::stod(f[3]);
            if (n > 4) p.inferSamples = std::stod(f[4]);
            pts_.push_back(p);
        } catch (...) { /* skip malformed line */ }
    }
}

TracePoint JobTraceWorkload::at(double t) const {
    TracePoint none;
    if (live_) {
        // Live: follow the newest sample of the running job (re-read the file at most every 50 ms of
        // wall time, however many servers ask), independent of simulation speed.
        const auto now = std::chrono::steady_clock::now();
        if (now - lastRead_ > std::chrono::milliseconds(50)) { refresh(); lastRead_ = now; }
        if (pts_.empty()) { none.phase = "waiting"; return none; }
        return pts_.back();
    }
    if (pts_.empty()) { none.phase = "waiting"; return none; }
    const double tt = (t - start_) / scale_;
    if (tt < 0) { none.phase = "waiting"; return none; }
    if (tt > pts_.back().t + 1e-9) {
        TracePoint last = pts_.back();
        last.u = 0.0; // the job is over: servers idle until the next one
        last.phase = "done";
        return last;
    }
    // first sample whose window ends at or after tt (each row covers the window ending at its t)
    const auto it = std::lower_bound(pts_.begin(), pts_.end(), tt,
                                     [](const TracePoint& p, double v) { return p.t < v; });
    return *it;
}

double JobTraceWorkload::u(double t) const { return at(t).u; }

double JobTraceWorkload::durationSim() const { return pts_.empty() ? 0.0 : pts_.back().t * scale_; }

std::string JobTraceWorkload::describe() const {
    std::ostringstream s;
    s << (live_ ? "live job trace " : "job trace ") << path_ << " (" << pts_.size() << " samples, x" << scale_ << " time)";
    return s.str();
}

} // namespace mf
