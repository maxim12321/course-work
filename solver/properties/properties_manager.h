#pragma once

#include "material.h"
#include "../utils/matrix.h"
#include "../utils/computation_grid.h"

#include <functional>
#include <map>
#include <utility>
#include <vector>

class Properties;

class PropertiesLoader {
  enum class PropertiesType {
    Tool,
    Plate,
    Backing,
    Method,
    HeatExchange,
  };
  static const std::map<std::string, PropertiesType> name_to_type_;

public:
  PropertiesLoader() = delete;
  static Properties Load(const std::string &properties_filepath);

private:
  static std::map<std::string, double *>
  CreateFieldsMap(Properties *properties);
  static void LoadProperties(std::ifstream &file, Properties *properties);
  static void LoadNamedProperties(
      std::ifstream &file,
      const std::map<std::string, double *> &property_name_to_value,
      size_t count);

  // Some properties came in non-ISU units of measurement
  static void NormalizePropertiesMeasurmentUnits(Properties *properties);

private:
  static const size_t tool_properties_count_ = 9;
  static const size_t plate_properties_count_ = 4;
  static const size_t backing_properties_count_ = 4;
  static const size_t method_properties_count_ = 7;
  static const size_t heat_exchange_properties_count_ = 7;
};

class Properties {
  friend class PropertiesLoader;
  friend class PropertiesWrapper;
  friend std::ostream &operator<<(std::ostream &out, const Properties &properties);
  
  Properties() = default;

public:
  // Matrix InitializeGrids(int nx, int nz);
  ComputationGrid InitComputationGrid(bool log);

  double GetDeltaX(int i, int k);
  double GetDeltaZ(int i, int k);

  double GetDensity(int x, int z, double temp);
  double GetHeatCapacity(int x, int z, double temp);
  double GetThermalConductivity(int x, int z, double temp);

  double GetToolHeight();
  double GetToolPenetration();
  double GetToolWaveHeight();
  double GetToolInitTemperature();

  // double GetHeatOutputX();
  // double GetHeatOutputZ(int x);

  // double GetHeatX(int x, int z);
  // double GetHeatZ(int x, int z);

  int GetGridHeight();
  int GetGridWidth();
  double GetEpsilon1();
  double GetEpsilon2();
  int GetMaxIterations();
  int GetTimeLayers();

private:
  std::vector<int> DistributeNodes(const std::vector<double> &lens, int total_nodes_count);
  std::vector<double> LogSegmentSplit(double segment_length, int parts);
  std::vector<double> UniformSegmentSplit(double segment_length, int parts);
  void ComputeDeltas(bool log);
  
private:
  std::vector<Material> materials_;

  // Material ID for each cell
  std::vector<std::vector<int>> materials_grid_;

  // Plate properties
  double plate_lenght_;
  double plate_height_;
  double plate_init_temp_;
  size_t plate_material_i_;

  // Backing properties
  double backing_length_;
  double backing_height_;
  double backing_init_temp_;
  size_t backing_material_i_;

  double total_height_;
  double height_without_penetration_;
  double tool_start_;
  double tool_finish_;

  // Tool properties
  double tool_radius_;
  double tool_height_;
  double tool_penetration_depth_;
  double tool_wave_height_;
  double tool_init_temp_;
  double tool_angular_velo_;
  size_t tool_material_i_;
  double friction_coef_;
  double f_z_;
  double f_x_;

  // nodes in [z_plate_start, z_plate_finish) belongs to plate zone
  int z_plate_finish_;
  int z_plate_start_;

  int z_tool_start_;
  
  // nodes in [x_tool_start_, x_tool_finish_) belongs to instrument zone
  int x_tool_start_;
  int x_tool_finish_;

  // Heat exchange properties
  double alpha1_;
  double alpha2_;
  double alpha3_;
  double alpha4_;
  double alpha6_;
  double alpha7_;
  double out_temp_;

  // Method properties
  double nx_;
  double nz_;
  double delta_t_;
  double eps1_;
  double eps2_;
  double max_iter_;
  double time_layers_;

  double magic_split_;

  // store deltas in vector<vector<>> because z_deltas in (z_plate_finish - 1) row differ depending on x  
  std::vector<std::vector<double>> new_delta_x_;
  std::vector<std::vector<double>> new_delta_z_;

  std::vector<double> delta_x_;
  std::vector<double> delta_z_;

  // For each cell stores x position of its center
  std::vector<double> x_position_;

  const double kLogFactor = 1.25;
};

std::ostream &operator<<(std::ostream &out, const Properties &properties);
