#include "properties_wrapper.h"

PropertiesWrapper::PropertiesWrapper(Properties *properties)
    : dt(properties->delta_t_),
      alpha1(properties->alpha1_),
      alpha2(properties->alpha2_),
      alpha3(properties->alpha3_),
      alpha4(properties->alpha4_),
      alpha6(properties->alpha6_),
      alpha7(properties->alpha7_),
      t_out(properties->out_temp_),
      z_plate_start(properties->z_plate_start_),
      z_plate_finish(properties->z_plate_finish_),
      z_tool_start(properties->z_tool_start_),
      x_tool_start(properties->x_tool_start_),
      x_tool_finish(properties->x_tool_finish_),
      magic_split(properties->magic_split_),
      w_fr_x(properties->friction_coef_ * properties->tool_angular_velo_ * properties->tool_radius_ * properties->f_x_),
      tool_center_pos_x_(std::max(properties->plate_lenght_, properties->backing_length_) / 2),
      w_fr_z_(properties->friction_coef_ * properties->tool_angular_velo_ * properties->f_z_ / 2),
      properties_(properties) {}

void PropertiesWrapper::UpdateTemperatureLambda(const Matrix& temp) {
    get_temp_ = [temp](int x, int z) -> double {
        return temp(x, z);
    };
}

double PropertiesWrapper::c(int x, int z) {
    auto temp = get_temp_(x, z);
    return properties_->GetHeatCapacity(x, z, temp);
}

double PropertiesWrapper::rho(int x, int z) {
    auto temp = get_temp_(x, z);
    return properties_->GetDensity(x, z, temp);
}

double PropertiesWrapper::lambda(double x, double z) {
    int x1 = std::floor(x);
    int x2 = std::ceil(x);

    int z1 = std::floor(z);
    int z2 = std::ceil(z);

    double first_temp = get_temp_(x1, z1);
    double first_value = properties_->GetThermalConductivity(x1, z1, first_temp);

    double second_temp = get_temp_(x2, z2);
    double second_value = properties_->GetThermalConductivity(x2, z2, second_temp);

    return 2 * first_value * second_value / (first_value + second_value);
}

double PropertiesWrapper::dx(int i, int k) {
    return properties_->GetDeltaX(i, k);
}

double PropertiesWrapper::dxb(int i, int k) {
    return (properties_->GetDeltaX(i, k) + properties_->GetDeltaX(i - 1, k)) / 2;
}

double PropertiesWrapper::dz(int i, int k) {
    return properties_->GetDeltaZ(i, k);
}

double PropertiesWrapper::dzb(int i, int k) {
    return (properties_->GetDeltaZ(i, k) + properties_->GetDeltaZ(i, k - 1)) / 2;
}

double PropertiesWrapper::w_fr_z(int i) {
    double radius = std::fabs(tool_center_pos_x_ - properties_->x_position_[i]);
    return w_fr_z_ * radius;
}