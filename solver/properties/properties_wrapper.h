#pragma once

#include "properties_manager.h"
#include <functional>

class PropertiesWrapper {
public:
  PropertiesWrapper(Properties *properties);
  
  void UpdateTemperatureLambda(const Matrix& temp);

  double c(int x, int z);
  double rho(int x, int z);
  double lambda(double x, double y);

  double dx(int i, int k);
  double dxb(int i, int k);
  double dz(int i, int k);
  double dzb(int i, int k);

  const double &dt;
  const double &alpha1;
  const double &alpha2;
  const double &alpha3;
  const double &alpha4;
  const double &alpha6;
  const double &alpha7;
  const double &t_out;

  double w_fr_x;
  double w_fr_z(int i);

  // nodes in [z_plate_start, z_plate_finish) belongs to plate zone
  const int &z_plate_finish;
  const int &z_plate_start;

  const int &z_tool_start;
  
  // nodes in [x_tool_start_, x_tool_finish_) belongs to instrument zone
  const int &x_tool_start;
  const int &x_tool_finish;

  const double &magic_split;

  Properties *properties_;

private:
  double tool_center_pos_x_;
  double w_fr_z_;

  std::function<double(int x, int z)> get_temp_;
};
