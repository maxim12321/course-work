#include "semi_implicit_solver.h"

#include <iostream>

#include <chrono>
#include <utility>
#include <cassert>

#include "../utils/math.h"

SemiImplicitSolver::SemiImplicitSolver(int p_rank,
                                       int p_size,
                                       Properties* properties,
                                       Callback callback)
    : SolverBase(SolverType::SemiImplicit, p_rank, p_size, properties, std::move(callback)),
      x_derivates_ids_grid_(nx_, std::vector<std::pair<int, int>>(nz_)),
      z_derivates_ids_grid_(nx_, std::vector<std::pair<int, int>>(nz_)),
      semi_next_temp_(nx_, nz_, z_plate_finish, x_tool_start, x_tool_finish),
      next_temp_(nx_, nz_, z_plate_finish, x_tool_start, x_tool_finish),
      row_system_(nx_),
      column_system_(nz_)
      {
  InitDerivativeApproxs();
}

void SemiImplicitSolver::InitDerivativeApproxs() {
  // Default approximations for forward derivates
  int default_forward_x = x_derivates_.size();
  x_derivates_.push_back([&](int x, int z, double from_temp, double to_temp) -> DerivateCoeffs {
    // lambda{x + 0.5, z} / dxb{x + 1} * (to_temp - from_temp)
    double l = lambda(x + 0.5, z) / dxb(x + 1, z);
    double c = static_cast<int>(is_zero(to_temp) || is_zero(from_temp));
    return {-l * c, l * c, l * (to_temp - from_temp)};
  });
  int default_forward_z = z_derivates_.size();
  z_derivates_.push_back([&](int x, int z, double from_temp, double to_temp) -> DerivateCoeffs{
    // lambda{x, z + 0.5} / dzb{z + 1} * (to_temp - from_temp)
    double l = lambda(x, z + 0.5) / dzb(x, z + 1);
    double c = static_cast<int>(is_zero(to_temp) || is_zero(from_temp));
    return {-l * c, l * c, l * (to_temp - from_temp)};
  });
  int default_backward_x = x_derivates_.size();
  x_derivates_.push_back([&](int x, int z, double from_temp, double to_temp) -> DerivateCoeffs {
    // lambda{x - 0.5, z} / dxb{x} * (to_temp - from_temp)
    double l = lambda(x - 0.5, z) / dxb(x, z);
    double c = static_cast<int>(is_zero(to_temp) || is_zero(from_temp));
    return {-l * c, l * c, l * (to_temp - from_temp)};
  });
  int default_backward_z = z_derivates_.size();
  z_derivates_.push_back([&](int x, int z, double from_temp, double to_temp) -> DerivateCoeffs{
    // lambda{x, z - 0.5} / dzb{z} * (to_temp - from_temp)
    double l = lambda(x, z - 0.5) / dzb(x, z);
    double c = static_cast<int>(is_zero(to_temp) || is_zero(from_temp));
    return {-l * c, l * c, l * (to_temp - from_temp)};
  });

  // Boundary condition for x backward derivate on Gamma_1
  int gamma1_x = x_derivates_.size();
  x_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    // alpha{1} * (T{x, z} - T{out})
    return {0.l, alpha1, -alpha1 * t_out} ;
  });
  // Boundary condition for x forward derivate on Gamma_2
  int gamma2_x = x_derivates_.size();
  x_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    // -alpha{2} * (T{x, z} - T{out})
    return std::make_tuple(-alpha2, 0.l, alpha2 * t_out);
  });
  // Boundary condition for z backward derivate on Gamma_3
  int gamma3_z = z_derivates_.size();
  z_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    // alpha{3} * (T{x, z} - T{out})
    return {0.l, alpha3, -alpha3 * t_out};
  });
  // Boundary condition for z forward derivate on Gamma_4
  int gamma4_z = z_derivates_.size();
  z_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    // -alpha{4} * (T{x, z} - T{out})
    return {-alpha4, 0.l, alpha4 * t_out};
  });
  // Boundary condition for z forward derivate on Gamma_5
  int gamma5_z = z_derivates_.size();
  z_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    return {0.l, 0.l, 0.l};
  });
  // Boundary condition for x backward derivate on Gamma_6
  int gamma6_x = x_derivates_.size();
  x_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    // alpha{6} * (T{x, z} - T{out})
    return {0.l, alpha6, -alpha6 * t_out};
  });
  // Boundary condition for x forward derivate on Gamma_7
  int gamma7_x = x_derivates_.size();
  x_derivates_.push_back([&](int x, int z, double, double) -> DerivateCoeffs{
    // -alpha{7} * (T{x, z} - T{out})
    return {-alpha7, 0.l, alpha7 * t_out};
  });

  // Boundary condition for x backward derivate on Gamma4 + Gamma6
  int gamma46_backward_x = x_derivates_.size();
  x_derivates_.push_back([&, default_backward_x](int x, int z, double from_temp, double to_temp) -> DerivateCoeffs{
    // magic_split * alpha6 * (T{x, z} - T{out}) +
    // + (1 - magic_split) * x_derivates_[default_x](x, z, from, to)
    // + (1 - magic_split) * w_fr_x
    auto coefs = x_derivates_[default_backward_x](x, z, from_temp, to_temp);
    std::get<0>(coefs) *= (1 - magic_split);
    std::get<1>(coefs) *= (1 - magic_split);
    std::get<2>(coefs) *= (1 - magic_split);

    std::get<1>(coefs) += magic_split * alpha6;

    std::get<2>(coefs) -= magic_split * alpha6 * t_out;

    return coefs;
  });
  // Boundary condition for x forward derivate on Gamma4 + Gamma7
  int gamma47_forward_x = x_derivates_.size();
  x_derivates_.push_back([&, default_forward_x](int x, int z, double from_temp, double to_temp) -> DerivateCoeffs{
    // -magic_split * alpha{7} * (T{x, z} - T{out}) + 
    // + (1 - magic_split) * x_derivates_[default_x](x, z, from, to)
    // + (1 - magic_split) * w_fr_x
    auto coefs = x_derivates_[default_forward_x](x, z, from_temp, to_temp);
    std::get<0>(coefs) *= (1 - magic_split);
    std::get<1>(coefs) *= (1 - magic_split);
    std::get<2>(coefs) *= (1 - magic_split);

    std::get<0>(coefs) -= magic_split * alpha7;

    std::get<2>(coefs) += magic_split * alpha7 * t_out;

    return coefs;
  });

  for (int i = 0; i < nx_; ++i) {
    for (int k = 0; k < nz_; ++k) {
      // set invalid ids for nodes outside calculation area
      if (k >= z_plate_finish && (i < x_tool_start || i >= x_tool_finish)) {
        x_derivates_ids_grid_[i][k] = {-1, -1};
        z_derivates_ids_grid_[i][k] = {-1, -1};
        continue;
      }

      // process x derivates
      if (i == 0) {
        // gamma1
        x_derivates_ids_grid_[i][k] = {default_forward_x, gamma1_x};
      } else if (i == x_tool_start && k == z_plate_finish - 1) {
        // gamma4 + gamma6
        x_derivates_ids_grid_[i][k] = {default_forward_x, gamma46_backward_x};
      } else if (i == x_tool_start && k >= z_plate_finish) {
        // gamma6
        x_derivates_ids_grid_[i][k] = {default_forward_x, gamma6_x};
      } else if (i == x_tool_finish - 1 && k == z_plate_finish - 1) {
        // gamma4 + gamma7
        x_derivates_ids_grid_[i][k] = {gamma47_forward_x, default_backward_x};
      } else if (i == x_tool_finish - 1 && k >= z_plate_finish) {
        // gamma7
        x_derivates_ids_grid_[i][k] = {gamma7_x, default_backward_x};
      } else if (i == nx_ - 1) {
        // gamma2
        x_derivates_ids_grid_[i][k] = {gamma2_x, default_backward_x};
      } else {
        // Default
        x_derivates_ids_grid_[i][k] = {default_forward_x, default_backward_x};
      }

      // process z derivates
      if (k == 0) {
        // gamma3
        z_derivates_ids_grid_[i][k] = {default_forward_z, gamma3_z};
      } else if (k == z_plate_finish - 1 && (i < x_tool_start || i >= x_tool_finish)) {
        // gamma4
        z_derivates_ids_grid_[i][k] = {gamma4_z, default_backward_z};
      } else if (k == nz_ - 1 && (i >= x_tool_start && i < x_tool_finish)) {
        // gamma5
        z_derivates_ids_grid_[i][k] = {gamma5_z, default_backward_z};
      } else {
        z_derivates_ids_grid_[i][k] = {default_forward_z, default_backward_z};
      }
    }
  }
}

void SemiImplicitSolver::Solve() {
  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < properties_->GetTimeLayers(); ++i) {
    CalculateNextLayer();
    if (i % 100 == 0) {
      // Print some debug info for long calculations
    }
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
}

void SemiImplicitSolver::CalculateNextLayer() {
  int iteration = 0;
  while (iteration < properties_->GetMaxIterations()) {
  // for (int iteration = 0; iteration < properties_->GetMaxIterations();
      //  ++iteration) {
    RunRowIters(previous_temp_, current_temp_);
    RunColumnIters(current_temp_, semi_next_temp_);

    if (HasConverged(current_temp_, next_temp_)) {
      std::swap(current_temp_, next_temp_);
      // current_temp_ = next_temp_;
      break;
    }
    // TODO assign more optimal
    std::swap(current_temp_, next_temp_);
    // current_temp_ = next_temp_;
    ++iteration;
  }
    std::cout << "Iterations " << iteration << std::endl;

  previous_temp_ = current_temp_;
  on_layer_ready_(current_temp_);
}

bool SemiImplicitSolver::HasConverged(const ComputationGrid& current, const ComputationGrid& next) {
  for (int i = 0; i < current_temp_.GetRowCount(); ++i) {
    for (int j = 0; j < current_temp_.GetColumnCount(); ++j) {
      double delta = std::fabs(current(i, j) - next(i, j));

      if (delta > properties_->GetEpsilon1() * current(i, j) +
                      properties_->GetEpsilon2()) {
        return false;
      }
    }
  }
  return true;
}

void SemiImplicitSolver::RunRowIters(const ComputationGrid &prev_layer, const ComputationGrid &prev_iter) {
  UpdateTemperatureLambda(prev_layer);

  for (size_t i = 0; i < semi_next_temp_.GetColumnCount(); ++i) {
    auto column = semi_next_temp_.GetColumn(i);
    BuildRowIterSystem(prev_iter, prev_layer, column);
    assert(TridiagonalAlgorithm(row_system_, column));
  }
}

void SemiImplicitSolver::BuildRowIterSystem(
    const ComputationGrid &prev_layer,
    const ComputationGrid &prev_iter,
    Vector &row) {
  row_system_.SetN(row.GetSize());

  int j = 0;
  for (auto it = row.begin(); it != row.end(); ++it) {
    auto [i, k] = it.Index();

    auto forward_x = XForwardDerivateCoeffs(i, k);
    auto backward_x = XBackwardDerivateCoeffs(i, k);

    double z_to_temp = k < nz_ - 1 ? prev_iter(i, k + 1) : 0.l;
    double z_from_temp = k > 0 ? prev_iter(i, k - 1) : 0.l;
    auto forward_z = ZForwardDerivateCoeffs(i, k, prev_iter(i, k), z_to_temp);
    auto backward_z = ZBackwardDerivateCoeffs(i, k, z_from_temp, prev_iter(i, k));

    auto heat_x_pred = static_cast<int>(
      (i == x_tool_start - 1 || i == x_tool_start || i == x_tool_finish - 1 || i == x_tool_finish) && 
      // (i == x_tool_start - 1 || i == x_tool_finish) && 
      (k >= z_tool_start && k < z_plate_finish)
    );
    auto heat_x_sign = 2 * static_cast<int>(i == x_tool_start - 1 || i == x_tool_finish - 1) - 1; 
    heat_x_sign = 1; 
    auto auxilary_coef = -magic_split * static_cast<int>((i == x_tool_start || i == x_tool_finish - 1) && k == z_plate_finish - 1) + 1;
    double Qx = dt * (w_fr_x / dz(i, k)) / dx(i, k);
    Qx *= heat_x_pred * auxilary_coef * heat_x_sign;
    
    auto heat_z_pred = static_cast<int>(
      (k == z_tool_start - 1 || k == z_tool_start) && (i >= x_tool_start && i < x_tool_finish)
      // (k == z_tool_start - 1) && (i >= x_tool_start && i < x_tool_finish)
    );
    auto heat_z_sign = 2 * static_cast<int>(k == z_tool_start - 1) - 1;
    heat_z_sign = 1;
    double Qz = dt * (w_fr_z(i) / dx(i, k)) / dz(i, k);
    Qz *= heat_z_pred * heat_z_sign;
    
    // T{i - 1, k}
    row_system_[j][0] = dt / dx(i, k) * std::get<0>(backward_x);

    // T{i, k}
    row_system_[j][1] = c(i, k) * rho(i, k);
    row_system_[j][1] -= dt / dx(i, k) * (std::get<0>(forward_x) - std::get<1>(backward_x));
    row_system_[j][1] -= dt / dz(i, k) * (std::get<0>(forward_z) - std::get<1>(backward_z));

    // T{i + 1, k}
    row_system_[j][2] = -dt / dx(i, k) * std::get<1>(forward_x);

    // C
    *it = c(i, k) * rho(i, k) * prev_layer(i, k);
    *it += dt / dx(i, k) * (std::get<2>(forward_x) - std::get<2>(backward_x));
    *it += dt / dz(i, k) * (std::get<2>(forward_z) - std::get<2>(backward_z));
    *it += Qx + Qz;

    ++j;
  }
  row_system_.LeftShiftRow(0);
}


void SemiImplicitSolver::RunColumnIters(const ComputationGrid &prev_iter, const ComputationGrid &semi_next_iter) {
  UpdateTemperatureLambda(prev_iter);

  for (size_t i = 0; i < next_temp_.GetRowCount(); ++i) {
    auto row = next_temp_.GetRow(i);
    BuildColumnIterSystem(prev_iter, semi_next_iter, row);
    assert(TridiagonalAlgorithm(column_system_, row));
  }
}

void SemiImplicitSolver::BuildColumnIterSystem(
    const ComputationGrid &prev_iter,
    const ComputationGrid &semi_next_iter,
    Vector &column) {
  column_system_.SetN(column.GetSize());

  int j = 0;
  for (auto it = column.begin(); it != column.end(); ++it) {
    auto [i, k] = it.Index();

    double x_to_temp = i < nx_ - 1 ? semi_next_iter(i + 1, k) : 0.l;
    double x_from_temp = i > 0 ? semi_next_iter(i - 1, k) : 0.l;
    auto forward_x_semi = XForwardDerivateCoeffs(i, k, semi_next_iter(i, k), x_to_temp);
    auto backward_x_semi = XBackwardDerivateCoeffs(i, k, x_from_temp, semi_next_iter(i, k));

    auto forward_z_next = ZForwardDerivateCoeffs(i, k);
    auto backward_z_next = ZBackwardDerivateCoeffs(i, k);

    double z_to_temp = k < nz_ - 1 ? prev_iter(i, k + 1) : 0.l;
    double z_from_temp = k > 0 ? prev_iter(i, k - 1) : 0.l;
    auto forward_z_prev = ZForwardDerivateCoeffs(i, k, prev_iter(i, k), z_to_temp);
    auto backward_z_prev = ZBackwardDerivateCoeffs(i, k, z_from_temp, prev_iter(i, k));

    auto dt_dx = dt / dx(i, k);
    auto dt_dz = dt / dz(i, k);

    // T{i, k - 1}
    column_system_[j][0] = 0.l;
    column_system_[j][0] += dt_dz * std::get<0>(backward_z_next);

    // T{i, k}
    column_system_[j][1] = 0.l;
    column_system_[j][1] += c(i, k) * rho(i, k);
    column_system_[j][1] -= dt_dx * (std::get<0>(forward_x_semi) - std::get<1>(backward_x_semi));
    column_system_[j][1] -= dt_dz * (std::get<0>(forward_z_next) - std::get<1>(backward_z_next));

    // T{i, k + 1}
    column_system_[j][2] = 0.l;
    column_system_[j][2] -= dt_dz * std::get<1>(forward_z_next);
    
    // C
    *it = 0.l;
    *it += c(i, k) * rho(i, k) * semi_next_iter(i, k);
    *it -= dt_dx * semi_next_iter(i, k) * (std::get<0>(forward_x_semi) - std::get<1>(backward_x_semi));
    *it -= dt_dz * semi_next_iter(i, k) * (std::get<0>(forward_z_prev) - std::get<1>(backward_z_prev));
    // *it += dt_dz * (std::get<2>(forward_z_next) - std::get<2>(backward_z_next));
    auto f_z_prev_cond = static_cast<int>(is_zero(std::get<0>(forward_z_prev)) && is_zero(std::get<1>(forward_z_prev)));
    auto b_z_prev_cond = static_cast<int>(is_zero(std::get<0>(backward_z_prev)) && is_zero(std::get<1>(backward_z_prev)));
    *it -= dt_dz * (f_z_prev_cond * std::get<2>(forward_z_prev) - b_z_prev_cond * std::get<2>(backward_z_prev));

    ++j;
  }

  column_system_.LeftShiftRow(0);
}


DerivateCoeffs SemiImplicitSolver::XForwardDerivateCoeffs(int x, int z, double from_temp, double to_temp) {
  auto id = x_derivates_ids_grid_[x][z].first;
  return x_derivates_[id](x, z, from_temp, to_temp);
}

DerivateCoeffs SemiImplicitSolver::XBackwardDerivateCoeffs(int x, int z, double from_temp, double to_temp) {
  auto id = x_derivates_ids_grid_[x][z].second;
  return x_derivates_[id](x, z, from_temp, to_temp);
}

DerivateCoeffs SemiImplicitSolver::ZForwardDerivateCoeffs(int x, int z, double from_temp, double to_temp) {
  auto id = z_derivates_ids_grid_[x][z].first;
  return z_derivates_[id](x, z, from_temp, to_temp);
}

DerivateCoeffs SemiImplicitSolver::ZBackwardDerivateCoeffs(int x, int z, double from_temp, double to_temp) {
  auto id = z_derivates_ids_grid_[x][z].second;
  return z_derivates_[id](x, z, from_temp, to_temp);
}

void SemiImplicitSolver::SchedulerRoutine() {
}

void SemiImplicitSolver::WorkerRoutine() {
}
