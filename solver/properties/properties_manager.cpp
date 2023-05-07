#include "properties_manager.h"
#include "../utils/strings.h"
#include "../utils/math.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <cmath>
#include <utility>

const std::map<std::string, PropertiesLoader::PropertiesType>
    PropertiesLoader::name_to_type_ = {
        {"plate", PropertiesType::Plate},
        {"tool", PropertiesType::Tool},
        {"backing", PropertiesType::Backing},
        {"method", PropertiesType::Method},
        {"heat_exchange", PropertiesType::HeatExchange},
};

Properties
PropertiesLoader::Load(const std::string &properties_filepath) {
  std::ifstream input;
  input.open(properties_filepath, std::ios::in);
  if (input.fail()) {
    std::cerr << "An error occurred while opening the input file" << std::endl;
    assert(false);
  }

  Properties properties;

  LoadProperties(input, &properties);
  NormalizePropertiesMeasurmentUnits(&properties);

  return properties;
}

std::map<std::string, double *>
PropertiesLoader::CreateFieldsMap(Properties *properties) {
  std::map<std::string, double *> name_to_field;
  name_to_field["Hplast"] = &properties->plate_height_;
  name_to_field["Lplast"] = &properties->plate_lenght_;
  name_to_field["TePlast_init"] = &properties->plate_init_temp_;
  name_to_field["mat_plast"] = nullptr;

  name_to_field["HSub"] = &properties->backing_height_;
  name_to_field["LSub"] = &properties->backing_length_;
  name_to_field["TSub_init"] = &properties->backing_init_temp_;
  name_to_field["mat_sub"] = nullptr;

  name_to_field["Ftool,x"] = &properties->f_x_;
  name_to_field["Ftool,z"] = &properties->f_z_;
  name_to_field["Htool"] = &properties->tool_height_;
  name_to_field["Hza_tool"] = &properties->tool_penetration_depth_;
  name_to_field["Rtool"] = &properties->tool_radius_;
  name_to_field["TeTool_init"] = &properties->tool_init_temp_;
  name_to_field["coef_fr"] = &properties->friction_coef_;
  name_to_field["mat_tool"] = nullptr;
  name_to_field["obTool"] = &properties->tool_angular_velo_;

  name_to_field["H"] = &properties->nz_;
  name_to_field["W"] = &properties->nx_;
  name_to_field["eps1"] = &properties->eps1_;
  name_to_field["eps2"] = &properties->eps2_;
  name_to_field["ki_max"] = &properties->max_iter_;
  name_to_field["ksloi_fin"] = &properties->time_layers_;
  name_to_field["tau"] = &properties->delta_t_;

  name_to_field["TeOut"] = &properties->out_temp_;
  name_to_field["alf_toG1"] = &properties->alpha1_;
  name_to_field["alf_toG2"] = &properties->alpha2_;
  name_to_field["alf_toG3"] = &properties->alpha3_;
  name_to_field["alf_toG4"] = &properties->alpha4_;
  name_to_field["alf_toG6"] = &properties->alpha6_;
  name_to_field["alf_toG7"] = &properties->alpha7_;

  return name_to_field;
}

void PropertiesLoader::LoadProperties(std::ifstream &file,
                                             Properties *properties) {
  auto name_to_field = CreateFieldsMap(properties);
  std::string type_name;
  int material_id;
  while (std::getline(file, type_name)) {
    if (name_to_type_.count(type_name) == 0) {
      std::cerr << "Unknown properties type " << type_name << std::endl;
      abort();
    }
    auto type = name_to_type_.at(type_name);
    switch (type) {
    case PropertiesType::Tool:
      material_id = properties->materials_.size();
      properties->materials_.emplace_back(file, material_id);
      properties->tool_material_i_ = material_id;
      LoadNamedProperties(file, name_to_field, tool_properties_count_);
      break;
    case PropertiesType::Plate:
      material_id = properties->materials_.size();
      properties->materials_.emplace_back(file, material_id);
      properties->plate_material_i_ = material_id;
      LoadNamedProperties(file, name_to_field, plate_properties_count_);
      break;
    case PropertiesType::Backing:
      material_id = properties->materials_.size();
      properties->materials_.emplace_back(file, material_id);
      properties->backing_material_i_ = material_id;
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

void PropertiesLoader::LoadNamedProperties(
    std::ifstream &file,
    const std::map<std::string, double *> &property_name_to_value,
    size_t count) {
  for (size_t i = 0; i < count; i++) {
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

void PropertiesLoader::NormalizePropertiesMeasurmentUnits(
    Properties *properties) {
  // Convert mm to m
  properties->tool_radius_ *= 1e-3;
  properties->tool_height_ *= 1e-3;
  properties->tool_penetration_depth_ *= 1e-3;
  properties->plate_lenght_ *= 1e-3;
  properties->plate_height_ *= 1e-3;
  properties->backing_length_ *= 1e-3;
  properties->backing_height_ *= 1e-3;
  // Convert rotation per minute to padian per second
  properties->tool_angular_velo_ = 2 * M_PI * properties->tool_angular_velo_ / 60;
}

std::vector<int> Properties::DistributeNodes(const std::vector<double> &lens, int total_nodes_count) {
  double total_len = 0;
  for (const auto &len : lens) {
    total_len += len;
  }
  std::vector<std::pair<double, size_t>> parts(lens.size());
  for (size_t i = 0; i < lens.size(); ++i) {
    parts[i] = {lens[i] / total_len, i};
  }
  std::sort(parts.begin(), parts.end());
  std::vector<int> nodes(parts.size());
  int used = 0;
  for (const auto &info : parts) {
    int n = static_cast<int>(std::ceil(info.first * total_nodes_count));
    n = std::min(n, total_nodes_count - used);
    used += n;
    nodes[info.second] = n;
  }

  return nodes;
}

void Properties::ComputeDeltas(bool log) {
  auto nx_parts = std::vector<double>{
    std::max(plate_lenght_, backing_length_) / 2 - tool_radius_,
    2 * tool_radius_,
    std::max(plate_lenght_, backing_length_) / 2 - tool_radius_,
  };
  auto nx_nodes = DistributeNodes(nx_parts, nx_ - 2 /*disregard boundary nodes*/);
  x_tool_start_ = nx_nodes[0] + 1 /*one more node for border Gamma1*/;
  x_tool_finish_ = x_tool_start_ + nx_nodes[1];

  for (size_t i = 0; i < nx_parts.size(); ++i) {
    auto deltas = log ? LogSegmentSplit(nx_parts[i], nx_nodes[i]) : UniformSegmentSplit(nx_parts[i], nx_nodes[i]);
    for (const auto &dx : deltas) {
      if (new_delta_x_.empty()) {
        // split first delta for border node
        new_delta_x_.push_back({0.25 * dx});
        
        new_delta_x_.push_back({0.75 * dx});
        continue;
      }
      new_delta_x_.push_back({dx});
    }
    if (i + 1 == nx_parts.size()) {
      // split last delta for border node
      double dx = new_delta_x_.back().back();
      new_delta_x_.back().back() = 0.75 * dx;
      new_delta_x_.push_back({0.25 * dx});
    }
  }
  assert(new_delta_x_.size() == nx_);

  size_t count_x = 0;
  double current_x_pos = 0.l;
  for (const auto &dxs : new_delta_x_) {
    current_x_pos += dxs.back() / 2;
    x_position_.push_back(current_x_pos);
    current_x_pos += dxs.back() / 2;
    count_x += dxs.size();
  }
  assert(count_x == nx_);
  assert(is_zero(std::fabs(std::max(plate_lenght_, backing_length_) - current_x_pos)));
  assert(x_position_.size() == nx_);

  auto nz_parts = std::vector<double>{
    backing_height_,
    plate_height_ - tool_penetration_depth_,
    tool_penetration_depth_,
    tool_height_ - tool_penetration_depth_,
  };
  auto nz_nodes = DistributeNodes(nz_parts, nz_ - 2 /*disregard boundary nodes*/);
  z_plate_start_ = nz_nodes[0] + 1 /*one more node for border Gamma3*/;
  z_tool_start_ = z_plate_start_ + nz_nodes[1];
  z_plate_finish_ = z_tool_start_ + nz_nodes[2] + 1 /*one more node for border Gamma4*/;

  bool split = false;
  for (size_t i = 0; i < nz_parts.size(); ++i) {
    auto deltas = log ? LogSegmentSplit(nz_parts[i], nz_nodes[i]) : UniformSegmentSplit(nz_parts[i], nz_nodes[i]);
    for (const auto &dz : deltas) {
      // split first delta for border node
      if (new_delta_z_.empty()) {
        new_delta_z_.push_back({0.25 * dz});
        new_delta_z_.push_back({0.75 * dz});
        continue;;
      }

      // split first row above plate
      if (split) {
        for (int x = x_tool_start_; x < x_tool_finish_; ++x) {
          new_delta_z_.back()[x] += 0.25 * dz;
          magic_split_ = 0.25 * dz / new_delta_z_.back()[x];
        }
        new_delta_z_.push_back({0.75 * dz});
        split = false;
        continue;
      }
      
      new_delta_z_.push_back({dz});
    }

    if (i + 2 == nz_parts.size()) {
      double dz = new_delta_z_.back().back();
      new_delta_z_.back().back() = dz * 0.75;
      new_delta_z_.push_back(std::vector<double>(nx_, 0.25 * dz));
      split = true;
    }
  }
  assert(new_delta_z_.size() == nz_);
  size_t count_z = 0;
  for (const auto &dzs : new_delta_z_) {
    count_z += dzs.size();
  }
  assert(count_z == nz_ + nx_ - 1);
  double first = 0.l;
  double second = 0.l;
  double third = 0.l;
  for (int z = 0; z < nz_; ++z) {
    second += GetDeltaZ(x_tool_start_, z);
    if (z < z_plate_finish_) {
      first += GetDeltaZ(0, z);
      third += GetDeltaZ(x_tool_finish_, z);
    }
  }
  assert(is_zero(std::fabs(first - backing_height_ - plate_height_)));
  assert(is_zero(std::fabs(first - third)));
  assert(is_zero(std::fabs(second - backing_height_ - plate_height_ - tool_height_ + tool_penetration_depth_)));
}

std::vector<double> Properties::UniformSegmentSplit(double segment_length, int parts) {
  double average = segment_length / parts;
  return std::vector<double>(parts, average);
}

std::vector<double> Properties::LogSegmentSplit(double segment_length, int parts) {
  std::vector<double> result(parts);

  if (parts == 1) {
    result[0] = segment_length;
    return result;
  }

  if (parts % 2 != 0) {
    double avg = segment_length / parts;
    segment_length -= avg;
    result[parts / 2] = avg;
  }

  double b = segment_length / 2 * (kLogFactor - 1) / (pow(kLogFactor, parts / 2) - 1);
  for (int i = 0; i < parts / 2; ++i) {
    result[i] = b;
    result[parts - 1 - i] = b;
    b *= kLogFactor;
  } 

  return result;
}


ComputationGrid Properties::InitComputationGrid(bool log) {
  ComputeDeltas(log);
  ComputationGrid init_temp(nx_, nz_, z_plate_finish_, x_tool_start_, x_tool_finish_);
  materials_grid_.assign(nx_, std::vector<int>(nz_));

  for (int x = 0; x < nx_; ++x) {
    for (int z = 0; z < nz_; ++z) {
      if (z < z_plate_start_) {
        init_temp(x, z) = backing_init_temp_;
        materials_grid_[x][z] = backing_material_i_;
      } else if (z < z_tool_start_) {
        init_temp(x, z) = plate_init_temp_;
        materials_grid_[x][z] = plate_material_i_;
      } else if (x >= x_tool_start_ && x < x_tool_finish_) {
        init_temp(x, z) = tool_init_temp_;
        materials_grid_[x][z] = tool_material_i_;
      } else if (z < z_plate_finish_) {
        init_temp(x, z) = plate_init_temp_;
        materials_grid_[x][z] = plate_material_i_;
      } else {
        // air zone
        init_temp(x, z) = 0;
        materials_grid_[x][z] = -1;
      }
    }
  }

  return init_temp;  

}

// Physical properties getters
double Properties::GetDensity(int x, int z, double temp) {
  auto id = materials_grid_[x][z];
  return materials_[id].GetDensity().ApproximateAt(temp);
}

double Properties::GetHeatCapacity(int x, int z, double temp) {
  auto id = materials_grid_[x][z];
  return materials_[id].GetHeatCapacity().ApproximateAt(temp);
}

double Properties::GetThermalConductivity(int x, int z, double temp) {
  auto id = materials_grid_[x][z];
  return materials_[id].GetThermalConductivity().ApproximateAt(temp);
}

double Properties::GetDeltaX(int i, int k) {
  auto &dxs = new_delta_x_[i];
  return dxs[k % dxs.size()];
}

double Properties::GetDeltaZ(int i, int k) {
  auto &dzs = new_delta_z_[k];
  return dzs[i % dzs.size()];
}

// Tool getters
double Properties::GetToolHeight() { return tool_height_; }

double Properties::GetToolPenetration() {
  return tool_penetration_depth_;
}

double Properties::GetToolWaveHeight() { return tool_wave_height_; }

double Properties::GetToolInitTemperature() { return tool_init_temp_; }

// Method getters
int Properties::GetGridHeight() { return nz_; }

int Properties::GetGridWidth() { return nx_; }

double Properties::GetEpsilon1() { return eps1_; }

double Properties::GetEpsilon2() { return eps2_; }

int Properties::GetMaxIterations() { return max_iter_; }

int Properties::GetTimeLayers() { return time_layers_; }

std::ostream &operator<<(std::ostream &out, const Properties &properties) {
  out << "--- Method ---" << std::endl;
  out << "Delta t = " << properties.delta_t_ << std::endl;
  out << "Eps1 = " << properties.eps1_ << std::endl;
  out << "Eps2 = " << properties.eps2_ << std::endl;
  out << "Max iter = " << properties.max_iter_ << std::endl;
  out << "Time layers = " << properties.time_layers_ << std::endl;
  out << "Grid height = " << properties.nz_ << std::endl;
  out << "Grid width = " << properties.nx_ << std::endl;
  out << std::endl;
  out << "--- Plate ---" << std::endl;
  out << "Height = " << properties.plate_height_ << std::endl;
  out << "Length = " << properties.plate_lenght_ << std::endl;
  out << "Init temp = " << properties.plate_init_temp_ << std::endl;
  out << "Density = "
            << properties.materials_[properties.plate_material_i_].GetDensity().
              ApproximateAt(properties.plate_init_temp_)
            << std::endl;
  out << "Heat capacity = "
            << properties.materials_[properties.plate_material_i_].GetHeatCapacity().
              ApproximateAt(properties.plate_init_temp_)
            << std::endl;
  out
      << "Thermal cond = "
      << properties.materials_[properties.plate_material_i_].GetThermalConductivity().
        ApproximateAt(properties.plate_init_temp_)
      << std::endl;
  out << std::endl;
  out << "--- Backing ---" << std::endl;
  out << "Height = " << properties.backing_height_ << std::endl;
  out << "Length = " << properties.backing_length_ << std::endl;
  out << "Init temp = " << properties.backing_init_temp_ << std::endl;
  out << "Density = "
            << properties.materials_[properties.backing_material_i_].GetDensity().
              ApproximateAt(properties.backing_init_temp_)
            << std::endl;
  out << "Heat capacity = "
            << properties.materials_[properties.backing_material_i_].GetHeatCapacity().
              ApproximateAt(properties.backing_init_temp_)
            << std::endl;
  out
      << "Thermal cond = "
      << properties.materials_[properties.backing_material_i_].GetThermalConductivity().
        ApproximateAt(properties.backing_init_temp_)
      << std::endl;
  out << std::endl;
  out << "--- Tool ---" << std::endl;
  out << "Height = " << properties.tool_height_ << std::endl;
  out << "Length = " << properties.tool_radius_ << std::endl;
  out << "Penetration = " << properties.tool_penetration_depth_ << std::endl;
  out << "Angular Velo = " << properties.tool_angular_velo_ << std::endl;
  out << "Friction = " << properties.friction_coef_ << std::endl;
  out << "Fx = " << properties.f_x_ << std::endl;
  out << "Fz = " << properties.f_z_ << std::endl;
  out << "Init temp = " << properties.tool_init_temp_ << std::endl;
  out << "Density = "
            << properties.materials_[properties.tool_material_i_].GetDensity().
              ApproximateAt(properties.tool_init_temp_)
            << std::endl;
  out << "Heat capacity = "
            << properties.materials_[properties.tool_material_i_].GetHeatCapacity().
              ApproximateAt(properties.tool_init_temp_)
            << std::endl;
  out
      << "Thermal cond = "
      << properties.materials_[properties.tool_material_i_].GetThermalConductivity().
        ApproximateAt(properties.tool_init_temp_)
      << std::endl;
  out << std::endl;
  out << "--- Heat Exchange ---" << std::endl;
  out << "alpha1 = " << properties.alpha1_ << std::endl;
  out << "alpha2 = " << properties.alpha2_ << std::endl;
  out << "alpha3 = " << properties.alpha3_ << std::endl;
  out << "alpha4 = " << properties.alpha4_ << std::endl;
  out << "alpha6 = " << properties.alpha6_ << std::endl;
  out << "alpha7 = " << properties.alpha7_ << std::endl;
  out << "Out temp = " << properties.out_temp_ << std::endl;

  return out;
}
