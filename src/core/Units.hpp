#pragma once
// Unit helpers. The physics core works in SI units throughout, except that temperatures are
// passed around in degrees Celsius for readability (converted to kelvin where a formula needs it).
namespace mf {

constexpr double kZeroCelsius = 273.15;     // K
constexpr double kGravity     = 9.80665;    // m/s^2
constexpr double kRackUnit    = 0.04445;    // m (1U = 1.75 in)

constexpr double toKelvin(double tC) { return tC + kZeroCelsius; }
constexpr double lpmToM3s(double lpm) { return lpm / 60000.0; }
constexpr double m3sToLpm(double m3s) { return m3s * 60000.0; }

} // namespace mf
