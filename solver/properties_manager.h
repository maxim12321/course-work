#pragma once

#include "material.h"
#include "utils/matrix.h"
#include "utils/vector.h"

#include <functional>
#include <map>
#include <utility>
#include <vector>

class PropertiesManager;

class PropertiesManagerFabric {
  enum class PropertiesType {
    Tool,
    Plate,
    Backing,
    Method,
    HeatExchange,
  };
  static const std::map<std::string, PropertiesType> name_to_type_;

public:
  PropertiesManagerFabric() = delete;
  static PropertiesManager Create(const std::string &properties_filepath);

private:
  static std::map<std::string, double *>
  CreateFieldsMap(PropertiesManager *manager);
  static void LoadProperties(std::ifstream &file, PropertiesManager *manager);
  static void LoadNamedProperties(
      std::ifstream &file,
      const std::map<std::string, double *> &property_name_to_value,
      size_t count);

  // Some properties came in non-ISU units of measurement
  static void NormalizePropertiesMeasurmentUnits(PropertiesManager *manager);

private:
  static const size_t tool_properties_count_ = 9;
  static const size_t plate_properties_count_ = 4;
  static const size_t backing_properties_count_ = 4;
  static const size_t method_properties_count_ = 8;
  static const size_t heat_exchange_properties_count_ = 6;
};

class PropertiesManager {
  friend class PropertiesManagerFabric;
  PropertiesManager() = default;

public:
  void PrintAllProperties();

  double GetAlpha1();
  double GetAlpha2();
  double GetAlpha3();
  double GetAlpha4();
  double GetAlpha4Tool();
  double GetOutTemperature();

  Matrix InitializeGrids(int nx, int nz);

  double GetDeltaT();

  double GetDeltaX(int i);
  double GetDeltaBackX(int i);
  double GetDeltaZ(int i);
  double GetDeltaBackZ(int i);

  double GetDensity(int x, int z, double temp);
  double GetHeatCapacity(int x, int z, double temp);
  double GetThermalConductivity(int x, int z, double temp);

  double GetToolHeight();
  double GetToolPenetration();
  double GetToolWaveHeight();
  double GetToolInitTemperature();
  int GetToolStartI();
  int GetToolFinishI();
  int GetBackingStartI();

  double GetHeatOutputX();
  double GetHeatOutputZ(int x);

  double GetHeatX(int x, int z);
  double GetHeatZ(int x, int z);

  int GetGridHeight();
  int GetGridWidth();
  double GetEpsilon1();
  double GetEpsilon2();
  int GetMaxIterations();
  int GetTimeLayers();

private:
  std::vector<int> ComputeDeltas(Vector &delta, Vector borders);

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
  double total_length_;
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

  // (Nx * Is) in doc
  int i_tool_start_;
  // (Nx * If) in doc
  int i_tool_finish_;
  // Useful for detecting nodes that adjacent to tool
  int i_backing_start_;
  // It's z coordinate for first node under the tool
  int i_tool_bottom_start_;
  // It's z coordinate for first node under the plate
  int i_plate_start_;

  // Heat exchange properties
  double alpha1_;
  double alpha2_;
  double alpha3_;
  double alpha4_;
  double alpha4_tool_;
  double out_temp_;

  // Method properties
  double nx_;
  double nz_;
  double delta_t_;
  double eps1_;
  double eps2_;
  double max_iter_;
  double time_layers_;
  double tool_words_;

  Vector delta_x_;
  Vector delta_z_;

  // For each cell stores x position of its center
  Vector x_position_;
};
