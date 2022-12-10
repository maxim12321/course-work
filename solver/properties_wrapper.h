#pragma once

#include "properties_manager.h"
#include <functional>

class PropertiesWrapper {
public:
  PropertiesWrapper(PropertiesManager *manager);
  
  void UpdateTemperatureLambda(const Matrix& temp);

  double c(int x, int z);
  double rho(int x, int z);
  double lambda(double x, double y);

  double dx(int i);
  double dxb(int i);
  double dz(int i);
  double dzb(int i);

  double dt;
  double alpha1;
  double alpha2;
  double alpha3;
  double alpha4;
  double t_out;

  PropertiesManager *manager_;

private:  
  std::function<double(int x, int z)> get_temp_;
};
