#include "properties_manager.h"

#include <QDebug>
#include <QFile>
#include <QStringList>
#include <QString>
#include <cmath>

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

std::pair<double, double> PropertiesManager::GridToObjectPosition(int x, int z) {
    double left_x = 0;
    double right_x = 0;
    double down_z = 0;
    double up_z;

    for (int i = 0; i <= x; ++i) {
        left_x = right_x;
        right_x += i < delta_x_.GetSize() ? delta_x_[i] : 0;
    }
    for (int i = 0; i <= z; ++i) {
        down_z = up_z;
        up_z += i < delta_z_.GetSize() ? delta_z_[i] : 0;
    }

    return std::make_pair((left_x + right_x) / 2, (down_z + up_z) / 2);
}

PropertiesManager::Object PropertiesManager::GetObjectInPositon(double x, double z) {
    int floor_x = std::floor(x);
    int ceil_x = std::ceil(x);
    int floor_z = std::floor(z);
    int ceil_z = std::ceil(z);

    double dx = x - static_cast<double>(floor_x);
    double dz = z - static_cast<double>(floor_z);

    auto floor_pos = GridToObjectPosition(floor_x, floor_z);
    std::pair<double, double> ceil_pos;
    if (floor_x == ceil_x && floor_z == ceil_z) {
        ceil_pos = floor_pos;
    } else {
        ceil_pos = GridToObjectPosition(ceil_x, ceil_z);
    }


    double pos_x = floor_pos.first * (1 - dx) + ceil_pos.first * dx;
    double pos_z = floor_pos.second * (1 - dz) + ceil_pos.second * dz;

    if (pos_z < backing_height_ || std::fabs(backing_height_ - pos_z) < 1e-9) {
        return Object::kBacking;
    } else if (pos_z < height_without_penetration_ || std::fabs(height_without_penetration_ - pos_z) < 1e-9){
        return Object::kPlate;
    } else if (pos_x < tool_start_ || pos_x > tool_finish_) {
        return Object::kPlate;
    } else {
        return Object::kTool;
    }
}

Matrix PropertiesManager::InitializeGrids(int nx, int nz) {
    total_length_ = std::max(plate_lenght_, backing_length_);
    // + 1 according to enumeration in doc
    delta_x_ = Vector(nx + 1, total_length_ / nx);
    delta_x_[0] = 0;

    total_height_ = plate_height_ + backing_height_;
    delta_z_ = Vector(nz + 1, total_height_ / nz);
    delta_z_[0] = 0;

    height_without_penetration_ = total_height_ - tool_penetration_depth_;
    tool_start_ = total_length_ / 2 - tool_radius_;
    tool_finish_ = total_length_ / 2 + tool_radius_;

    // + 2 because there are nodes on borders;
    init_temp_ = Matrix(nx + 2, nz + 2);

    for (size_t x = 0; x < nx + 2; ++x) {
        for (size_t z = 0; z < nz + 2; ++z) {
            auto object = GetObjectInPositon(x, z);

            switch (object) {
            case Object::kPlate:
                init_temp_[x][z] = plate_init_temp_;
                break;
            case Object::kBacking:
                init_temp_[x][z] = backing_init_temp_;
                break;
            case Object::kTool:
                init_temp_[x][z] = tool_init_temp_;
                break;
            }

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


double PropertiesManager::GetInitTemperature(int x, int z) {
    return init_temp_[x][z];
}

double PropertiesManager::GetDensity(double x, double z) {
    auto object = GetObjectInPositon(x, z);

    switch (object) {
    case Object::kPlate:
        return plate_material_.density;
    case Object::kBacking:
        return backing_material_.density;
    case Object::kTool:
        return tool_material_.density;
    }
}

double PropertiesManager::GetThermalConductivity(double x, double z) {
    auto object = GetObjectInPositon(x, z);

    switch (object) {
    case Object::kPlate:
        return plate_material_.thermal_conductivity;
    case Object::kBacking:
        return backing_material_.thermal_conductivity;
    case Object::kTool:
        return tool_material_.thermal_conductivity;
    }
}

double PropertiesManager::GetHeatCapacity(double x, double z) {
    auto object = GetObjectInPositon(x, z);

    switch (object) {
    case Object::kPlate:
        return plate_material_.heat_capacity;
    case Object::kBacking:
        return backing_material_.heat_capacity;
    case Object::kTool:
        return tool_material_.heat_capacity;
    }
}

double PropertiesManager::GetDeltaT() {
    return delta_t_;
}

double PropertiesManager::GetDeltaX(int i) {
    return delta_x_[i];
}
double PropertiesManager::GetDeltaZ(int i) {
    return delta_z_[i];
}
