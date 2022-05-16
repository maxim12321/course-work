#include "material.h"

#include <QDebug>
#include <QStringList>
#include <assert.h>

Material::Material(int id) : id_(id) {
    InitializeMaterial();
}

int Material::GetId() const {
    return id_;
}

const TableValue& Material::GetDensity() const {
    return density_;
}

const TableValue& Material::GetHeatCapacity() const {
    return heat_capacity_;
}

const TableValue& Material::GetThermalConductivity() const {
    return thermal_conductivity_;
}

void Material::InitializeMaterial() {
    QFile file(kMaterialsPath + QString::number(id_) + ".dat");
    if(!file.open(QIODevice::ReadOnly)) {
        qDebug() << "error opening file for material " << id_;
        return;
    }

    ReadTableValues(&file, &density_);
    ReadTableValues(&file, &heat_capacity_);
    ReadTableValues(&file, &thermal_conductivity_);

    file.close();
}

void Material::ReadTableValues(QFile* file, TableValue* table_value) {
    assert(!file->atEnd());
    QString keys_line = file->readLine().trimmed();
    QStringList keys = keys_line.split(" ");

    assert(!file->atEnd());
    QString values_line = file->readLine().trimmed();
    QStringList values = values_line.split(" ");

    assert(keys.size() == values.size());
    assert(!keys.empty());

    for (size_t i = 0; i < keys.size(); ++i) {
        table_value->AddTableValue(keys[i].toDouble(), values[i].toDouble());
    }
}
