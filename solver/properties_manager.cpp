#include "properties_manager.h"

#include <QDebug>
#include <QFile>
#include <QStringList>
#include <QString>

PropertiesManager::PropertiesManager() {
    LoadMaterials();
}

void PropertiesManager::LoadMaterials() {
    QFile file(kMaterialsConfigFile);
    if(!file.open(QIODevice::ReadOnly)) {
        qDebug() << "error opening file: " << file.error();
        return;
    }

    file.readLine();
    while (!file.atEnd()) {
        QString line = file.readLine().trimmed();
        QStringList properties = line.split(" ");

        if (properties.empty() || (properties.size() == 1 && properties[0].trimmed().isEmpty())) {
            continue;
        }
        if (properties.size() != 4) {
            qDebug() << "Wrong line in materials file\n";
            continue;
        }

        int material_id = properties[0].toInt();
        materials_[material_id].density = properties[1].toDouble();
        materials_[material_id].heat_capacity = properties[2].toDouble();
        materials_[material_id].thermal_conductivity = properties[3].toDouble();
    }

    file.close();
}

Matrix PropertiesManager::InitializeGrids(int nx, int nz) {
    double length = std::max(plate_lenght_, backing_length_);
    // + 1 according to enumeration in doc
    delta_x_ = Vector(nx + 1, length / nx);
    delta_x_[0] = 0;

    double height = plate_height_ + backing_height_;
    delta_z_ = Vector(nz + 1, height / nz);
    delta_z_[0] = 0;

    // + 2 because there are nodes on borders;
    init_temp_ = Matrix(nx + 2, nz + 2);
    for (size_t x = 0; x < nx + 2; ++x) {
        for (size_t z = 0; z < nz + 2; ++z) {
        // init temp in node using its position in space and init temps of obejcts;
        }
    }

    return init_temp_;
}

void PropertiesManager::SetPlateProperties(double length, double height, double init_temp, int material) {
    if (materials_.count(material) == 0) {
        qDebug() << "Unknown plate material!";
        return;
    }
    plate_material_ = materials_[material];
    plate_lenght_ = length;
    plate_height_ = height;
    plate_init_temp_ = init_temp;
}

void PropertiesManager::SetBackingProperties(double length, double height, double init_temp, int material) {
    if (materials_.count(material) == 0) {
        qDebug() << "Unknown backing material!";
        return;
    }
    backing_material_ = materials_[material];
    backing_length_ = length;
    backing_height_ = height;
    backing_init_temp_ = init_temp;
}

void PropertiesManager::SetToolProperties(double radius, double height, double penetration_depth, double init_temp,
                                         double angular_velo, double friction_coef, double f_z, double f_x, int material) {
    if (materials_.count(material) == 0) {
        qDebug() << "Unknown tool material!";
        return;
    }
    tool_material_ = materials_[material];
    tool_radius_ = radius;
    tool_height_ = height;
    tool_penetration_depth_ = penetration_depth;
    tool_init_temp_ = init_temp;
    tool_angular_velo_ = angular_velo;
    friction_coef_ = friction_coef;
    f_z_ = f_z;
    f_x_ = f_x;
}

void PropertiesManager::SetMethodProperties(double delta_t, double eps1, double eps2, int max_iter_count, int time_layers_count, int tool_words) {
    delta_t_ = delta_t;
    eps1_ = eps1;
    eps2_ = eps2;
    max_iter_ = max_iter_count;
    tool_words_ = tool_words;
}

void PropertiesManager::SetHeatExchangePropeties(double alpha_1, double alpha_2, double alpha_3, double alpha_4, double alpha_4_tool, double out_temp) {
    alpha1_ = alpha_1;
    alpha2_ = alpha_2;
    alpha3_ = alpha_3;
    alpha4_ = alpha_4;
    alpha4_tool_ = alpha_4_tool;
    out_temp_ = out_temp;
}

double PropertiesManager::GetAlpha1() {
    return alpha1_;
}

double PropertiesManager::GetAlpha2() {
    return alpha2_;
}

double PropertiesManager::GetAlpha3() {
    return alpha3_;
}

double PropertiesManager::GetAlpha4() {
    return alpha4_;
}

double PropertiesManager::GetAlpha4Tool() {
    return alpha4_tool_;
}

double PropertiesManager::GetOutTemperature() {
    return out_temp_;
}


double PropertiesManager::GetInitTemperature(double x, double z) {
    return init_temp_[0][0];  // Not implemented
}

double PropertiesManager::GetDensity(double x, double z) {
    return 0;  // Not implemented
}

double PropertiesManager::GetThermalConductivity(double x, double z) {
    return 0;  // Not implemented
}

double PropertiesManager::GetHeatCapacity(double x, double z) {
    return tool_material_.heat_capacity;  // Not implemented
}
