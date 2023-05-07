#pragma once

#include <functional>
#include <unordered_map>

#include "properties/properties_manager.h"
#include "properties/properties_wrapper.h"
#include "utils/computation_grid.h"

enum class SolverType {
  SemiImplicit,
  Explicit
};

const std::string kSemiImplicitSolverName = "semi-implicit";
const std::string kExplicitSolverName = "explicit";

const std::unordered_map<std::string, SolverType> kSolverNameToType = {
  {kSemiImplicitSolverName, SolverType::SemiImplicit},
  {kExplicitSolverName, SolverType::Explicit}
};

class SolverBase : public PropertiesWrapper {
 public:
  using Callback = std::function<void(const Matrix&)>;

  SolverBase(SolverType type, int p_rank, int p_size, Properties* properties, Callback callback)
      : PropertiesWrapper(properties),
        type_(type),
        p_rank_(p_rank),
        p_size_(p_size),
        nx_(properties->GetGridWidth()),
        nz_(properties->GetGridHeight()),
        properties_(properties),
        on_layer_ready_(std::move(callback)) {
    Initialize();
  }

  virtual void Solve() = 0;

 private:
  // Initialize fields with values from Properties...
  void Initialize() {
    // current_temp_ = properties_->InitComputationGrid(type_ == SolverType::SemiImplicit);
    current_temp_ = properties_->InitComputationGrid(false);
    previous_temp_ = current_temp_;
    on_layer_ready_(current_temp_);
  }

 protected:
  const int kSchedulerRank = 0;

  SolverType type_;

  int p_rank_;
  int p_size_;

  int nx_;
  int nz_;

  ComputationGrid current_temp_;
  ComputationGrid previous_temp_;

  Properties* properties_;
  Callback on_layer_ready_;
};
