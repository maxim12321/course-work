#include "properties_wrapper.h"

PropertiesWrapper::PropertiesWrapper(PropertiesManager *manager)
    : dt(manager->GetDeltaT()),
      alpha1(manager->GetAlpha1()),
      alpha2(manager->GetAlpha2()),
      alpha3(manager->GetAlpha3()),
      alpha4(manager->GetAlpha4()),
      t_out(manager->GetOutTemperature()),
      manager_(manager) {}

void PropertiesWrapper::UpdateTemperatureLambda(const Matrix& temp) {
    get_temp_ = [temp](int x, int z) -> double {
        return temp(x, z);
    };
}

double PropertiesWrapper::c(int x, int z) {
    auto temp = get_temp_(x, z);
    return manager_->GetHeatCapacity(x, z, temp);
}

double PropertiesWrapper::rho(int x, int z) {
    auto temp = get_temp_(x, z);
    return manager_->GetDensity(x, z, temp);
}

double PropertiesWrapper::lambda(double x, double z) {
    int x1 = std::floor(x);
    int x2 = std::ceil(x);

    int z1 = std::floor(z);
    int z2 = std::ceil(z);

    double first_temp = get_temp_(x1, z1);
    double first_value = manager_->GetThermalConductivity(x1, z1, first_temp);

    double second_temp = get_temp_(x2, z2);
    double second_value = manager_->GetThermalConductivity(x2, z2, second_temp);

    return 2 * first_value * second_value / (first_value + second_value);
}

double PropertiesWrapper::dx(int i) {
    return manager_->GetDeltaX(i);
}

double PropertiesWrapper::dxb(int i) {
    return manager_->GetDeltaBackX(i);
}

double PropertiesWrapper::dz(int i) {
    return manager_->GetDeltaZ(i);
}

double PropertiesWrapper::dzb(int i) {
    return manager_->GetDeltaBackZ(i);
}
