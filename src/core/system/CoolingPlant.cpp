#include "system/CoolingPlant.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>

namespace mf {

// ---- epsilon-NTU ---------------------------------------------------------------------------------
double hx::effectivenessCounterflow(double ntu, double cr) {
    if (ntu <= 0) return 0.0;
    if (std::abs(1.0 - cr) < 1e-9) return ntu / (1.0 + ntu);
    const double e = std::exp(-ntu * (1.0 - cr));
    return (1.0 - e) / (1.0 - cr * e);
}

double hx::effectivenessCrossflow(double ntu, double cr, bool cminIsMixed) {
    if (ntu <= 0) return 0.0;
    if (cr < 1e-9) return 1.0 - std::exp(-ntu);
    if (!cminIsMixed) // Cmax mixed, Cmin unmixed
        return (1.0 / cr) * (1.0 - std::exp(-cr * (1.0 - std::exp(-ntu))));
    // Cmin mixed, Cmax unmixed
    return 1.0 - std::exp(-(1.0 / cr) * (1.0 - std::exp(-cr * ntu)));
}

double WeatherParams::at(double t) const {
    if (amplitudeK == 0.0) return meanC;
    return meanC - amplitudeK * std::cos(2.0 * M_PI * (t - phaseS) / periodS);
}

namespace {
// Design-point reference properties: plain water (never the scaled coolant), so a more viscous
// equipment coolant is penalised relative to the design.
const FluidProps& refWater(double tC) {
    static WaterIAPWS w;
    static const FluidProps p35 = w.props(35.0), p30 = w.props(30.0);
    return tC > 32.5 ? p35 : p30;
}
double sideHA(double hARef, double m, double mRef, const FluidProps& p, const FluidProps& ref, double n) {
    const double x = std::max(m / mRef, 1e-6);
    return hARef * std::pow(x, n) * std::pow(p.mu / ref.mu, 1.0 / 3.0 - n) * std::pow(p.k / ref.k, 2.0 / 3.0) *
           std::pow(p.cp / ref.cp, 1.0 / 3.0);
}
double implicitMix(double tOld, double cap, double dt, double flowC, double tIn) {
    const double cdt = cap / dt;
    return (cdt * tOld + flowC * tIn) / (cdt + flowC);
}
} // namespace

CoolingPlant::CoolingPlant(const PlantParams& p, std::shared_ptr<const Fluid> equipFluid,
                           std::shared_ptr<const Fluid> facilityFluid)
    : p_(p), fe_(std::move(equipFluid)), ff_(std::move(facilityFluid)) {
    tES_ = tER_ = p_.supplySetpointC;
    tFS_ = tFR_ = p_.supplySetpointC - 4.0;
    s_.equipSupplyT = tES_;
}

double CoolingPlant::facilityMdot() const {
    return ff_->props(0.5 * (tFS_ + tFR_)).rho * lpmToM3s(p_.facility.flowLpm);
}

void CoolingPlant::setFanSpeed(double s) { fanSpeed_ = std::clamp(s, p_.dryCooler.minFanSpeed, 1.0); }

void CoolingPlant::cduPerformance(double mE, double tE, double tF, double& ua, double& eps) const {
    const CduParams& c = p_.cdu;
    const FluidProps pe = fe_->props(tE), pf = ff_->props(tF);
    const FluidProps& re = refWater(35.0);
    const FluidProps& rf = refWater(30.0);
    const double mF = facilityMdot();
    const double hAe = sideHA(c.hAEquipRef, mE, re.rho * lpmToM3s(c.equipRefLpm), pe, re, c.exponent);
    const double hAf = sideHA(c.hAFacilityRef, mF, rf.rho * lpmToM3s(c.facilityRefLpm), pf, rf, c.exponent);
    ua = 1.0 / (1.0 / hAe + 1.0 / hAf + c.wallR);
    const double ce = mE * pe.cp, cf = mF * pf.cp;
    const double cmin = std::min(ce, cf), cmax = std::max(ce, cf);
    eps = cmin > 0 ? hx::effectivenessCounterflow(ua / cmin, cmin / cmax) : 0.0;
}

void CoolingPlant::dryCoolerPerformance(double s, double tW, double, double& ua, double& eps) const {
    const DryCoolerParams& d = p_.dryCooler;
    const double hAairRef = d.uaRef / d.airShare, hAwRef = d.uaRef / (1.0 - d.airShare);
    const FluidProps pw = ff_->props(tW);
    const double mF = facilityMdot();
    const double mRef = refWater(30.0).rho * lpmToM3s(p_.facility.flowLpm);
    const double hAair = hAairRef * std::pow(std::max(s, 1e-3), d.airExponent);
    const double hAw = hAwRef * std::pow(mF / mRef, d.waterExponent);
    ua = 1.0 / (1.0 / hAair + 1.0 / hAw);
    const double cAir = d.rhoAir * d.cpAir * d.maxAirflow * std::max(s, 1e-3);
    const double cW = mF * pw.cp;
    const double cmin = std::min(cAir, cW), cmax = std::max(cAir, cW);
    eps = hx::effectivenessCrossflow(ua / cmin, cmin / cmax, /*cminIsMixed=*/cAir < cW);
}

double CoolingPlant::chillerCop(double tChilled, double tOut) const {
    const ChillerParams& c = p_.chiller;
    const double te = toKelvin(tChilled - c.evapApproachK), tc = toKelvin(tOut + c.condApproachK);
    if (tc <= te + 1.0) return c.maxCop;
    return std::min(c.maxCop, c.carnotFraction * te / (tc - te));
}

double CoolingPlant::facilityPumpPower() const {
    const FacilityParams& f = p_.facility;
    const double dp = p_.cdu.facilityDpRefPa + p_.dryCooler.coilDpRefPa + f.pipingDpRefPa; // at design flow
    return dp * lpmToM3s(f.flowLpm) / (f.pumpEta * f.motorEta);
}

void CoolingPlant::step(double dt, double time, double mE, double rackReturnT) {
    if (!(dt > 0.0)) return;
    const double tOut = p_.weather.at(time);
    const FluidProps pe = fe_->props(0.5 * (tES_ + tER_));
    const FluidProps pf = ff_->props(0.5 * (tFS_ + tFR_));
    const double ce = mE * pe.cp, mF = facilityMdot(), cf = mF * pf.cp;
    const double capES = pe.rho * p_.equipSupplyVolumeL * 1e-3 * pe.cp;
    const double capER = pe.rho * p_.equipReturnVolumeL * 1e-3 * pe.cp;
    const double capF = pf.rho * 0.5 * p_.facilityVolumeL * 1e-3 * pf.cp;
    const double old[4] = {tES_, tER_, tFS_, tFR_};

    // 1. coolant returning from the rack fills the return volume
    tER_ = implicitMix(tER_, capER, dt, ce, rackReturnT);
    // 2. CDU
    double ua, eps;
    cduPerformance(mE, 0.5 * (tER_ + tES_), 0.5 * (tFS_ + tFR_), ua, eps);
    const double cmin = std::min(ce, cf);
    const double qCdu = eps * cmin * (tER_ - tFS_);
    const double hotOut = ce > 0 ? tER_ - qCdu / ce : tER_;
    const double coldOut = tFS_ + qCdu / cf;
    // 3. supply volume, 4. facility return volume
    tES_ = implicitMix(tES_, capES, dt, ce, hotOut);
    tFR_ = implicitMix(tFR_, capF, dt, cf, coldOut);
    // 5. dry cooler
    double uaDc, epsDc;
    dryCoolerPerformance(fanSpeed_, tFR_, tOut, uaDc, epsDc);
    const double cAir = p_.dryCooler.rhoAir * p_.dryCooler.cpAir * p_.dryCooler.maxAirflow * fanSpeed_;
    // Bypass: when outdoor air is no colder than the water, running the coil would heat the loop,
    // so the water bypasses the dry cooler and its fans stop.
    const double qDcFans = epsDc * std::min(cAir, cf) * (tFR_ - tOut);
    // Economizer rule: while the chiller is carrying load with the fans at full speed, keep the
    // dry cooler only if the chiller electricity it saves exceeds its fan power.
    bool econBypass = false;
    if (p_.chiller.enabled && fanSpeed_ > 0.98 && tES_ > p_.supplySetpointC && qDcFans > 0) {
        const double fanFull = p_.dryCooler.maxFanPower * std::pow(fanSpeed_, 3);
        econBypass = qDcFans / chillerCop(tFS_, tOut) < fanFull;
    }
    const bool bypass = tFR_ <= tOut || econBypass;
    const double qDc = bypass ? 0.0 : qDcFans;
    const double dcOut = tFR_ - qDc / cf;
    // 6. fan PI on the equipment supply temperature
    const double e = tES_ - p_.supplySetpointC;
    fanSpeed_ = std::clamp(fanSpeed_ + p_.fanKp * (e - lastErr_) + p_.fanKi * e * dt, p_.dryCooler.minFanSpeed, 1.0);
    lastErr_ = e;
    // 7. trim chiller: only with fans at (almost) full speed; aims 0.1 K above the setpoint so the
    //    fan loop stays saturated instead of chattering with the chiller.
    double qCh = 0.0;
    if (p_.chiller.enabled && fanSpeed_ > 0.98 && tES_ > p_.supplySetpointC && eps > 0) {
        const double target = p_.supplySetpointC + 0.1;
        const double tReq = tER_ - (tER_ - target) * ce / (eps * cmin);
        qCh = std::clamp(cf * (dcOut - tReq), 0.0, p_.chiller.capacityW);
    }
    const double chOut = dcOut - qCh / cf;
    // 8. facility supply volume
    tFS_ = implicitMix(tFS_, capF, dt, cf, chOut);

    // ---- report ----
    const double stor = (capES * (tES_ - old[0]) + capER * (tER_ - old[1]) + capF * (tFS_ - old[2]) + capF * (tFR_ - old[3])) / dt;
    s_.equipSupplyT = tES_;
    s_.equipReturnT = rackReturnT;
    s_.cduHotInT = tER_;
    s_.cduHotOutT = hotOut;
    s_.facilitySupplyT = tFS_;
    s_.facilityReturnT = tFR_;
    s_.dryCoolerOutT = dcOut;
    s_.outdoorT = tOut;
    s_.cduHeat = qCdu;
    s_.dryCoolerHeat = qDc;
    s_.chillerHeat = qCh;
    s_.cduUA = ua;
    s_.cduEffectiveness = eps;
    s_.dryCoolerUA = uaDc;
    s_.dryCoolerEffectiveness = epsDc;
    s_.fanSpeed = bypass ? 0.0 : fanSpeed_;
    s_.fanPower = bypass ? 0.0 : p_.dryCooler.maxFanPower * std::pow(fanSpeed_, 3);
    s_.chillerCop = qCh > 0 ? chillerCop(chOut, tOut) : 0.0;
    s_.chillerPower = qCh > 0 ? qCh / s_.chillerCop : 0.0;
    s_.facilityPumpPower = facilityPumpPower();
    s_.residual = ce * (rackReturnT - tES_) - qDc - qCh - stor;
}

void CoolingPlant::solveSteady(double time, double mE, double q) {
    const double tOut = p_.weather.at(time);
    const double set = p_.supplySetpointC;

    struct Result { double tES, tER, tFS, tFR, dcOut, qDc, qCh, ua, eps, uaDc, epsDc; };
    // Steady loop temperatures at fan speed s; with `chill`, the chiller holds the setpoint.
    auto solveAt = [&](double s, bool chill) {
        Result r{};
        r.tES = tES_; r.tER = tER_; r.tFS = tFS_; r.tFR = tFR_;
        for (int it = 0; it < 8; ++it) { // properties and effectiveness depend weakly on temperature
            const double ce = mE * fe_->props(0.5 * (r.tES + r.tER)).cp;
            const double cf = facilityMdot() * ff_->props(0.5 * (r.tFS + r.tFR)).cp;
            cduPerformance(mE, 0.5 * (r.tES + r.tER), 0.5 * (r.tFS + r.tFR), r.ua, r.eps);
            dryCoolerPerformance(s, r.tFR, tOut, r.uaDc, r.epsDc);
            const double cAir = p_.dryCooler.rhoAir * p_.dryCooler.cpAir * p_.dryCooler.maxAirflow * s;
            const double gDc = r.epsDc * std::min(cAir, cf); // dry cooler heat per K (water in - air)
            const double gCdu = r.eps * std::min(ce, cf);
            if (!chill) {
                r.qDc = q; r.qCh = 0.0;
                r.tFR = tOut + q / gDc;
                r.tFS = r.tFR - q / cf;
                r.dcOut = r.tFS;
            } else {
                r.tES = set;
                r.tER = set + q / ce;
                const double tFSreq = r.tER - q / gCdu;
                const double tFR = tFSreq + q / cf;
                const double qDc = std::max(0.0, gDc * (tFR - tOut)); // bypass if air is warmer
                double qCh = q - qDc;
                if (qCh > p_.chiller.capacityW) qCh = p_.chiller.capacityW;
                if (qCh < 0) qCh = 0;
                r.qCh = qCh;
                r.qDc = q - qCh;
                r.tFR = r.qDc > 0 ? tOut + r.qDc / gDc : tFSreq + q / cf;
                r.dcOut = r.tFR - r.qDc / cf;
                r.tFS = r.dcOut - r.qCh / cf;
            }
            r.tER = r.tFS + q / gCdu;
            r.tES = r.tER - q / ce;
        }
        return r;
    };

    const double smin = p_.dryCooler.minFanSpeed;
    double s;
    Result r = solveAt(smin, false);
    if (r.tES <= set) {
        s = smin; // cold weather: minimum fans already undershoot the setpoint
    } else {
        r = solveAt(1.0, false);
        if (r.tES >= set) {
            s = 1.0;
            if (p_.chiller.enabled) {
                r = solveAt(1.0, true);
                // Economizer: compare with bypassing the dry cooler entirely (fans off, chiller only).
                const Result b = solveAt(1e-6, true);
                auto elec = [&](const Result& x, double fanW) {
                    return fanW + (x.qCh > 0 ? x.qCh / chillerCop(x.tFS, tOut) : 0.0);
                };
                if (elec(b, 0.0) < elec(r, p_.dryCooler.maxFanPower)) { r = b; r.qDc = 0.0; }
            }
        } else {
            double lo = smin, hi = 1.0;
            for (int it = 0; it < 60; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (solveAt(mid, false).tES > set) lo = mid; else hi = mid;
            }
            s = hi;
            r = solveAt(s, false);
        }
    }
    fanSpeed_ = s;
    lastErr_ = r.tES - set;
    tES_ = r.tES; tER_ = r.tER; tFS_ = r.tFS; tFR_ = r.tFR;

    s_.equipSupplyT = r.tES;
    s_.equipReturnT = r.tER;
    s_.cduHotInT = r.tER;
    s_.cduHotOutT = r.tES;
    s_.facilitySupplyT = r.tFS;
    s_.facilityReturnT = r.tFR;
    s_.dryCoolerOutT = r.dcOut;
    s_.outdoorT = tOut;
    s_.cduHeat = q;
    s_.dryCoolerHeat = r.qDc;
    s_.chillerHeat = r.qCh;
    s_.cduUA = r.ua;
    s_.cduEffectiveness = r.eps;
    s_.dryCoolerUA = r.uaDc;
    s_.dryCoolerEffectiveness = r.epsDc;
    const bool bypass = r.qDc <= 0.0;
    s_.fanSpeed = bypass ? 0.0 : s;
    s_.fanPower = bypass ? 0.0 : p_.dryCooler.maxFanPower * s * s * s;
    s_.chillerCop = r.qCh > 0 ? chillerCop(r.tFS, tOut) : 0.0;
    s_.chillerPower = r.qCh > 0 ? r.qCh / s_.chillerCop : 0.0;
    s_.facilityPumpPower = facilityPumpPower();
    s_.residual = q - r.qDc - r.qCh;
}

} // namespace mf
