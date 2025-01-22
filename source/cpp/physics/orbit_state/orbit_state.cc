#include "lupnt/physics/orbit_state.h"

#include <magic_enum/magic_enum.hpp>
namespace lupnt {

  // ****************************************************************************
  // OrbitState
  // ****************************************************************************

  OrbitState::OrbitState(const Vec6 &x, Frame coord, OrbitStateRepres repres,
                         const std::array<const char *, 6> &names,
                         const std::array<const char *, 6> &units)
      : x_(x), frame_(coord), repres_(repres), names_(names), units_(units) {}

  // Overrides
  int OrbitState::GetSize() const { return kOrbitStateSize; }
  VecX OrbitState::GetVec() const { return x_; }
  void OrbitState::SetVec(const VecX &x) { x_ = x; }
  Real OrbitState::GetValue(int i) const { return x_(i); }
  void OrbitState::SetValue(int idx, Real val) { x_(idx) = val; }
  StateType OrbitState::GetStateType() const { return static_cast<StateType>(repres_); }

  Vec6 OrbitState::GetVec6() const { return x_; }
  Frame OrbitState::GetFrame() const { return frame_; }
  std::array<const char *, 6> OrbitState::GetUnits() const { return units_; }
  std::array<const char *, 6> OrbitState::GetNames() const { return names_; }
  OrbitStateRepres OrbitState::GetOrbitStateRepres() const { return repres_; }
  OrbitState OrbitState::CreateCopyWithValue(const Vec6 &x) const {
    return OrbitState(x, frame_, repres_, names_, units_);
  }

  void OrbitState::SetOrbitStateRepres(const OrbitStateRepres rep) { repres_ = rep; }
  void OrbitState::SetCoordSystem(Frame frame) { frame_ = frame; }

  Real OrbitState::operator()(int idx) const { return x_(idx); }
  std::ostream &OrbitState::operator<<(std::ostream &os) const {
    os << "<OrbitState(" << x_.transpose() << ", " << frame_ << ", "
       << magic_enum::enum_name(repres_) << ")>";
    return os;
  }

  // ****************************************************************************
  // CartesianOrbitState
  // ****************************************************************************

  void CheckOrbitStateRepres(const OrbitState &state, OrbitStateRepres repres) {
    if (state.GetOrbitStateRepres() != repres)
      throw std::runtime_error("OrbitState type must be "
                               + std::string(magic_enum::enum_name(repres)));
  }

}  // namespace lupnt
