#pragma once

#include <fstream>
#include <string>

#include "table_value.h"

class Material {
public:
    Material() = default;
    // Material properties the properties of the metric come in the input file
    explicit Material(std::ifstream& file, int id);

    int GetId() const;

    const TableValue& GetDensity() const;
    const TableValue& GetHeatCapacity() const;
    const TableValue& GetThermalConductivity() const;

private:
    void InitializeMaterial();

    void ReadTableValues(std::ifstream& file, TableValue* table_value);

private:
    int id_;

    TableValue density_;
    TableValue heat_capacity_;
    TableValue thermal_conductivity_;
};
