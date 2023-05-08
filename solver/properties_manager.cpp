#include "properties_manager.h"
#include "utils/strings.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>
#include <utility>

const std::map<std::string, PropertiesManagerFabric::PropertiesType>
    PropertiesManagerFabric::name_to_type_ = {
        {"plate", PropertiesType::Plate},
        {"tool", PropertiesType::Tool},
        {"backing", PropertiesType::Backing},
        {"method", PropertiesType::Method},
        {"heat_exchange", PropertiesType::HeatExchange},
};

PropertiesManager
PropertiesManagerFabric::Create(const std::string &properties_filepath) {
  std::ifstream input;
  input.open(properties_filepath, std::ios::in);
  if (input.fail()) {
    std::cerr << "An error occurred while opening the input file" << std::endl;
    assert(false);
  }

  PropertiesManager manager;

  LoadProperties(input, &manager);
  NormalizePropertiesMeasurmentUnits(&manager);

  return manager;
}

std::map<std::string, double *>
PropertiesManagerFabric::CreateFieldsMap(PropertiesManager *manager) {
  std::map<std::string, double *> name_to_field;
  name_to_field["Hplast"] = &manager->plate_height_;
  name_to_field["Lplast"] = &manager->plate_lenght_;
  name_to_field["TePlast_init"] = &manager->plate_init_temp_;
  name_to_field["mat_plast"] = nullptr;

  name_to_field["HSub"] = &manager->backing_height_;
  name_to_field["LSub"] = &manager->backing_length_;
  name_to_field["TSub_init"] = &manager->backing_init_temp_;
  name_to_field["mat_sub"] = nullptr;

  name_to_field["Ftool,x"] = &manager->f_x_;
  name_to_field["Ftool,z"] = &manager->f_z_;
  name_to_field["Htool"] = &manager->tool_height_;
  name_to_field["Hza_tool"] = &manager->tool_penetration_depth_;
  name_to_field["Rtool"] = &manager->tool_radius_;
  name_to_field["TeTool_init"] = &manager->tool_init_temp_;
  name_to_field["coef_fr"] = &manager->friction_coef_;
  name_to_field["mat_tool"] = nullptr;
  name_to_field["obTool"] = &manager->tool_angular_velo_;

  name_to_field["H"] = &manager->nz_;
  name_to_field["W"] = &manager->nx_;
  name_to_field["eps1"] = &manager->eps1_;
  name_to_field["eps2"] = &manager->eps2_;
  name_to_field["ki_max"] = &manager->max_iter_;
  name_to_field["ksloi_fin"] = &manager->time_layers_;
  name_to_field["nzTool"] = &manager->tool_words_;
  name_to_field["tau"] = &manager->delta_t_;

  name_to_field["TeOut"] = &manager->out_temp_;
  name_to_field["alf_toG1"] = &manager->alpha1_;
  name_to_field["alf_toG2"] = &manager->alpha2_;
  name_to_field["alf_toG3"] = &manager->alpha3_;
  name_to_field["alf_toG4g"] = &manager->alpha4_;
  name_to_field["alf_toG4gi"] = &manager->alpha4_tool_;

  return std::move(name_to_field);
}

void PropertiesManagerFabric::LoadProperties(std::ifstream &file,
                                             PropertiesManager *manager) {
  auto name_to_field = CreateFieldsMap(manager);
  std::string type_name;
  int material_id;
  while (std::getline(file, type_name)) {
    auto type = name_to_type_.at(type_name);
    switch (type) {
    case PropertiesType::Tool:
      material_id = manager->materials_.size();
      manager->materials_.emplace_back(file, material_id);
      manager->tool_material_i_ = material_id;
      LoadNamedProperties(file, name_to_field, tool_properties_count_);
      break;
    case PropertiesType::Plate:
      material_id = manager->materials_.size();
      manager->materials_.emplace_back(file, material_id);
      manager->plate_material_i_ = material_id;
      LoadNamedProperties(file, name_to_field, plate_properties_count_);
      break;
    case PropertiesType::Backing:
      material_id = manager->materials_.size();
      manager->materials_.emplace_back(file, material_id);
      manager->backing_material_i_ = material_id;
      LoadNamedProperties(file, name_to_field, backing_properties_count_);
      break;
    case PropertiesType::Method:
      LoadNamedProperties(file, name_to_field, method_properties_count_);
      break;
    case PropertiesType::HeatExchange:
      LoadNamedProperties(file, name_to_field, heat_exchange_properties_count_);
      break;
    default:
      assert(false);
      break;
    }
  }
}

void PropertiesManagerFabric::LoadNamedProperties(
    std::ifstream &file,
    const std::map<std::string, double *> &property_name_to_value,
    size_t count) {
  for (int i = 0; i < count; i++) {
    std::string line = ReadLine(file);
    auto key_value = SplitString(line, "|");
    assert(key_value.size() == 2);
    auto key = key_value[0];
    auto value = std::stod(key_value[1]);

    if (property_name_to_value.count(key) == 0) {
      std::cout << "Key [" << key << "] is missed from name to value map"
                << std::endl;
      assert(property_name_to_value.count(key) > 0);
    }
    auto ptr = property_name_to_value.at(key);
    if (ptr != nullptr) {
      *ptr = value;
    }
  }
}

void PropertiesManagerFabric::NormalizePropertiesMeasurmentUnits(
    PropertiesManager *manager) {
  // Convert mm to m
  manager->tool_radius_ *= 1e-3;
  manager->tool_height_ *= 1e-3;
  manager->tool_penetration_depth_ *= 1e-3;
  manager->plate_lenght_ *= 1e-3;
  manager->plate_height_ *= 1e-3;
  manager->backing_length_ *= 1e-3;
  manager->backing_height_ *= 1e-3;
  // Convert rotation per minute to padian per second
  manager->tool_angular_velo_ = 2 * M_PI * manager->tool_angular_velo_ / 60;
}

std::pair<int, int> PropertiesManager::ComputeDeltas(int n, long double dx,
                                                     Vector &delta,
                                                     Vector borders) {
  long double remainder = 0;
  long double curr_pos = 0;

  int j = 0;
  for (int i = 1; i <= n; ++i) {
    delta[i] = std::min(dx + remainder, borders[j] - curr_pos);
    if (delta[i] < 1e-9) {
      borders[j] = i - 1;
      ++j;
      if (j == borders.GetSize()) {
        std::cout << curr_pos << remainder << '\n';
        assert(i == n && curr_pos == borders[2] && remainder < 1e-9);
      }
      --i;
      continue;
    }
    remainder = remainder > 1e-9 ? 0 : dx - delta[i];
    curr_pos += delta[i];
  }
  return std::make_pair(static_cast<int>(borders[0]),
                        static_cast<int>(borders[1]));
}

Matrix PropertiesManager::InitializeGrids(int nx, int nz) {
  total_length_ = std::max(plate_lenght_, backing_length_);
  total_height_ = plate_height_ + backing_height_;
  height_without_penetration_ = total_height_ - tool_penetration_depth_;

  tool_start_ = total_length_ / 2 - tool_radius_;
  tool_finish_ = total_length_ / 2 + tool_radius_;

  // + 2 according to enumeration in doc
  delta_x_ = Vector(nx + 2);
  auto hor_indxs = ComputeDeltas(nx, total_length_ / nx, delta_x_,
                                 {tool_start_, tool_finish_, total_length_});
  i_tool_start_ = hor_indxs.first;
  i_tool_finish_ = hor_indxs.second;

  double current_x_position = 0;
  x_position_ = Vector(nx + 2);
  for (size_t i = 0; i < nx + 2; ++i) {
    x_position_[i] = current_x_position + delta_x_[i] / 2;
    current_x_position += delta_x_[i];
  }

  delta_z_ = Vector(nz + 2);
  auto vert_indxs = ComputeDeltas(
      nz, total_height_ / nz, delta_z_,
      {backing_height_, height_without_penetration_, total_height_});
  i_plate_start_ = vert_indxs.first;
  i_tool_bottom_start_ = vert_indxs.second;

  tool_wave_height_ =
      tool_height_ - tool_penetration_depth_ + 0.25 * delta_z_[nz];

  // + 2 because there are nodes on borders;
  Matrix init_temp(nx + 2, nz + 2);
  materials_grid_.assign(nx + 2, std::vector<int>(nz + 2));

  for (size_t x = 0; x < nx + 2; ++x) {
    for (size_t z = 0; z < nz + 2; ++z) {
      if (z <= i_plate_start_) {
        init_temp[x][z] = backing_init_temp_;
        materials_grid_[x][z] = backing_material_i_;
      } else if (z <= i_tool_bottom_start_) {
        init_temp[x][z] = plate_init_temp_;
        materials_grid_[x][z] = plate_material_i_;
      } else if (x <= i_tool_start_ || x > i_tool_finish_) {
        init_temp[x][z] = plate_init_temp_;
        materials_grid_[x][z] = plate_material_i_;
      } else {
        init_temp[x][z] = tool_init_temp_;
        materials_grid_[x][z] = tool_material_i_;
      }
    }
  }

  return init_temp;
}

void PropertiesManager::PrintAllProperties() {
  std::cout << "--- Method ---" << std::endl;
  std::cout << "Delta t = " << delta_t_ << std::endl;
  std::cout << "Eps1 = " << eps1_ << std::endl;
  std::cout << "Eps2 = " << eps2_ << std::endl;
  std::cout << "Max iter = " << max_iter_ << std::endl;
  std::cout << "Time layers = " << time_layers_ << std::endl;
  std::cout << "Tool layers = " << tool_words_ << std::endl;
  std::cout << std::endl;
  std::cout << "--- Plate ---" << std::endl;
  std::cout << "Height = " << plate_height_ << std::endl;
  std::cout << "Length = " << plate_lenght_ << std::endl;
  std::cout << "Init temp = " << plate_init_temp_ << std::endl;
  std::cout << "Density = "
            << materials_[plate_material_i_].GetDensity().ApproximateAt(
                   plate_init_temp_)
            << std::endl;
  std::cout << "Heat capacity = "
            << materials_[plate_material_i_].GetHeatCapacity().ApproximateAt(
                   plate_init_temp_)
            << std::endl;
  std::cout
      << "Thermal cond = "
      << materials_[plate_material_i_].GetThermalConductivity().ApproximateAt(
             plate_init_temp_)
      << std::endl;
  std::cout << std::endl;
  std::cout << "--- Backing ---" << std::endl;
  std::cout << "Height = " << backing_height_ << std::endl;
  std::cout << "Length = " << backing_length_ << std::endl;
  std::cout << "Init temp = " << backing_init_temp_ << std::endl;
  std::cout << "Density = "
            << materials_[backing_material_i_].GetDensity().ApproximateAt(
                   backing_init_temp_)
            << std::endl;
  std::cout << "Heat capacity = "
            << materials_[backing_material_i_].GetHeatCapacity().ApproximateAt(
                   backing_init_temp_)
            << std::endl;
  std::cout
      << "Thermal cond = "
      << materials_[backing_material_i_].GetThermalConductivity().ApproximateAt(
             backing_init_temp_)
      << std::endl;
  std::cout << std::endl;
  std::cout << "--- Tool ---" << std::endl;
  std::cout << "Height = " << tool_height_ << std::endl;
  std::cout << "Length = " << tool_radius_ << std::endl;
  std::cout << "Penetration = " << tool_penetration_depth_ << std::endl;
  std::cout << "Angular Velo = " << tool_angular_velo_ << std::endl;
  std::cout << "Friction = " << friction_coef_ << std::endl;
  std::cout << "Fx = " << f_x_ << std::endl;
  std::cout << "Fz = " << f_z_ << std::endl;
  std::cout << "Init temp = " << tool_init_temp_ << std::endl;
  std::cout << "Density = "
            << materials_[tool_material_i_].GetDensity().ApproximateAt(
                   tool_init_temp_)
            << std::endl;
  std::cout << "Heat capacity = "
            << materials_[tool_material_i_].GetHeatCapacity().ApproximateAt(
                   tool_init_temp_)
            << std::endl;
  std::cout
      << "Thermal cond = "
      << materials_[tool_material_i_].GetThermalConductivity().ApproximateAt(
             tool_init_temp_)
      << std::endl;
  std::cout << std::endl;
  std::cout << "--- Heat Exchange ---" << std::endl;
  std::cout << "alpha1 = " << alpha1_ << std::endl;
  std::cout << "alpha2 = " << alpha2_ << std::endl;
  std::cout << "alpha3 = " << alpha3_ << std::endl;
  std::cout << "alpha4 = " << alpha4_ << std::endl;
  std::cout << "alpha4_tool = " << alpha4_tool_ << std::endl;
  std::cout << "Out temp = " << out_temp_ << std::endl;
}

// Alpha getters
double PropertiesManager::GetAlpha1() { return alpha1_; }

double PropertiesManager::GetAlpha2() { return alpha2_; }

double PropertiesManager::GetAlpha3() { return alpha3_; }

double PropertiesManager::GetAlpha4() { return alpha4_; }

double PropertiesManager::GetAlpha4Tool() { return alpha4_tool_; }

double PropertiesManager::GetOutTemperature() { return out_temp_; }

// Physical properties getters
double PropertiesManager::GetDensity(int x, int z, double temp) {
  auto id = materials_grid_[x][z];
  return materials_[id].GetDensity().ApproximateAt(temp);
}

double PropertiesManager::GetHeatCapacity(int x, int z, double temp) {
  auto id = materials_grid_[x][z];
  return materials_[id].GetHeatCapacity().ApproximateAt(temp);
}

double PropertiesManager::GetThermalConductivity(int x, int z, double temp) {
  auto id = materials_grid_[x][z];
  return materials_[id].GetThermalConductivity().ApproximateAt(temp);
}

// Delta getters
double PropertiesManager::GetDeltaT() { return delta_t_; }

double PropertiesManager::GetDeltaX(int i) {
  assert(i > 0 && i + 1 < delta_x_.GetSize());
  return delta_x_[i];
}

double PropertiesManager::GetDeltaBackX(int i) {
  assert(i > 0 && i < delta_x_.GetSize());
  return (delta_x_[i] + delta_x_[i - 1]) / 2;
}

double PropertiesManager::GetDeltaZ(int i) {
  assert(i > 0 && i + 1 < delta_z_.GetSize());
  return delta_z_[i];
}

double PropertiesManager::GetDeltaBackZ(int i) {
  assert(i > 0 && i < delta_z_.GetSize());
  return (delta_z_[i] + delta_z_[i - 1]) / 2;
}

// Tool getters
double PropertiesManager::GetToolHeight() { return tool_height_; }

double PropertiesManager::GetToolPenetration() {
  return tool_penetration_depth_;
}

double PropertiesManager::GetToolWaveHeight() { return tool_wave_height_; }

double PropertiesManager::GetToolInitTemperature() { return tool_init_temp_; }

int PropertiesManager::GetToolStartI() { return i_tool_start_; }

int PropertiesManager::GetToolFinishI() { return i_tool_finish_; }

// Heat getters
double PropertiesManager::GetHeatOutputX() {
  // TODO: check values
  double pressure = f_x_ / delta_z_[i_tool_bottom_start_];
  double velocity = tool_angular_velo_ * tool_radius_;
  return friction_coef_ * velocity * pressure;
}

double PropertiesManager::GetHeatOutputZ(int x) {
  // TODO: check values
  double tool_center = tool_start_ + tool_radius_;
  double radius = std::fabs(tool_center - x_position_[x]);

  double pressure = f_z_ / (delta_x_[x] * 2);
  double velocity = tool_angular_velo_ * radius;
  return friction_coef_ * velocity * pressure;
}

double PropertiesManager::GetHeatX(int x, int z) {
  // Ignoring nodes below tool
  if (z < i_tool_bottom_start_) {
    return 0;
  }

  if (x == i_tool_start_ || x == i_tool_finish_ + 1) {
    return GetHeatOutputX() / delta_x_[x];
  }
  return 0;
}

double PropertiesManager::GetHeatZ(int x, int z) {
  // Ignoring node in tool
  if (z != i_tool_bottom_start_) {
    return 0;
  }
  if (x < i_tool_start_ || x > i_tool_finish_ + 2) {
    return 0;
  }

  return GetHeatOutputZ(x) / delta_z_[z];
}

// Method getters
int PropertiesManager::GetGridHeight() { return nz_; }

int PropertiesManager::GetGridWidth() { return nx_; }

double PropertiesManager::GetEpsilon1() { return eps1_; }

double PropertiesManager::GetEpsilon2() { return eps2_; }

int PropertiesManager::GetMaxIterations() { return max_iter_; }

int PropertiesManager::GetTimeLayers() { return time_layers_; }
