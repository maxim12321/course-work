#pragma once

#include <QFile>
#include <QString>

#include "table_value.h"

class Material {
public:
    // Material properties stored at resources/properties/materials/*id*.dat
    explicit Material(int id);

    int GetId() const;

    const TableValue& GetDensity() const;
    const TableValue& GetHeatCapacity() const;
    const TableValue& GetThermalConductivity() const;

private:
    const QString kMaterialsPath = "://resources/properties/materials/";

private:
    void InitializeMaterial();

    void ReadTableValues(QFile* file, TableValue* table_value);

private:
    int id_;

    TableValue density_;
    TableValue heat_capacity_;
    TableValue thermal_conductivity_;
};
