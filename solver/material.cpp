#include "material.h"
#include "utils/strings.h"

#include <cassert>
#include <fstream>
#include <iostream>

Material::Material(std::ifstream &file, int id) : id_(id) {
  ReadTableValues(file, &density_);
  ReadTableValues(file, &heat_capacity_);
  ReadTableValues(file, &thermal_conductivity_);
}

int Material::GetId() const { return id_; }

const TableValue &Material::GetDensity() const { return density_; }

const TableValue &Material::GetHeatCapacity() const { return heat_capacity_; }

const TableValue &Material::GetThermalConductivity() const {
  return thermal_conductivity_;
}

void Material::InitializeMaterial() {}

void Material::ReadTableValues(std::ifstream &file, TableValue *table_value) {
  std::string keys_line = ReadLine(file);
  std::vector<std::string> keys = SplitString(keys_line, " ");

  std::string values_line = ReadLine(file);
  std::vector<std::string> values = SplitString(values_line, " ");

  assert(keys.size() == values.size());
  assert(!keys.empty());

  for (size_t i = 0; i < keys.size(); ++i) {
    double key = std::stod(keys[i]);
    double value = std::stod(values[i]);
    table_value->AddTableValue(key, value);
  }
}
