#include "properties_manager.h"

#include <QDebug>
#include <QFile>
#include <QStringList>
#include <QString>
#include <cmath>
#include <utility>

PropertiesManager::PropertiesManager() {
    LoadMaterials();
}

void PropertiesManager::LoadMaterials() {
//    QFile file(kMaterialsConfigFile);
//    if(!file.open(QIODevice::ReadOnly)) {
//        qDebug() << "error opening file: " << file.error();
//        return;
//    }

//    file.readLine();
//    while (!file.atEnd()) {
//        QString line = file.readLine().trimmed();
//        QStringList properties = line.split(" ");

//        if (properties.empty() || (properties.size() == 1 && properties[0].trimmed().isEmpty())) {
//            continue;
//        }
//        if (properties.size() != 4) {
//            qDebug() << "Wrong line in materials file\n";
//            continue;
//        }

//        int material_id = properties[0].toInt();
//        materials_[material_id].density = properties[1].toDouble();
//        materials_[material_id].heat_capacity = properties[2].toDouble();
//        materials_[material_id].thermal_conductivity = properties[3].toDouble();
//    }

//    file.close();

    materials_[backing_material_.GetId()] = &backing_material_;
    materials_[plate_material_.GetId()] = &plate_material_;
    materials_[tool_material_.GetId()] = &tool_material_;
}

std::pair<int, int> PropertiesManager::ComputeDeltas(int n, long double dx, Vector &delta, Vector borders) {
    long double remainder = 0;
    long double curr_pos = 0;
    
    int j = 0;
    for (int i = 1; i <= n; ++i) {
        delta[i] = std::min(dx + remainder, borders[j] - curr_pos);
        if (delta[i] < 1e-9) {
            borders[j] = i - 1;
            ++j;
            if (j == borders.GetSize()) {
                assert(i == n);
            }
            --i;
            continue;
        }
        remainder = dx - delta[i];
        curr_pos += delta[i];
    }
    return std::make_pair(static_cast<int>(borders[0]), static_cast<int>(borders[1]));
}

Matrix PropertiesManager::InitializeGrids(int nx, int nz) {
    total_length_ = std::max(plate_lenght_, backing_length_);
    total_height_ = plate_height_ + backing_height_;    
    height_without_penetration_ = total_height_ - tool_penetration_depth_;
    
    tool_start_ = total_length_ / 2 - tool_radius_;
    tool_finish_ = total_length_ / 2 + tool_radius_;
    
    // + 2 according to enumeration in doc
    delta_x_ = Vector(nx + 2);
    auto hor_indxs = ComputeDeltas(nx, total_length_ / nx, delta_x_, {tool_start_, tool_finish_, total_length_});
    i_tool_start_ = hor_indxs.first;
    i_tool_finish_ = hor_indxs.second;
    
    delta_z_ = Vector(nz + 2);
    auto vert_indxs = ComputeDeltas(nz, total_height_ / nz, delta_z_, {backing_height_, height_without_penetration_, total_height_});
    i_plate_start_ = vert_indxs.first;
    i_tool_bottom_start_ = vert_indxs.second;
    
    tool_wave_height_ = tool_height_ - tool_penetration_depth_ + 0.25 * delta_z_[nz];

    // + 2 because there are nodes on borders;
    Matrix init_temp(nx + 2, nz + 2);
    materials_grid_.fill(QVector<int>(nz + 2), nx + 2);

    for (size_t x = 0; x < nx + 2; ++x) {
        for (size_t z = 0; z < nz + 2; ++z) {
            if (z <= i_plate_start_) {
                init_temp[x][z] = backing_init_temp_;
                materials_grid_[x][z] = backing_material_.GetId();
            } else if (z <= i_tool_bottom_start_) {
                init_temp[x][z] = plate_init_temp_;
                materials_grid_[x][z] = plate_material_.GetId();
            } else if (x <= i_tool_start_ || x > i_tool_finish_) {
                init_temp[x][z] = plate_init_temp_;
                materials_grid_[x][z] = plate_material_.GetId();
            } else {
                init_temp[x][z] = tool_init_temp_;
                materials_grid_[x][z] = tool_material_.GetId();
            }
        }
    }

    return init_temp;
}

// Properties setters
void PropertiesManager::SetPlateProperties(double length, double height, double init_temp, int /*material*/) {
//    if (materials_.count(material) == 0) {
//        qDebug() << "Unknown plate material!";
//        return;
//    }
//    plate_material_ = materials_[material];
    plate_lenght_ = length;
    plate_height_ = height;
    plate_init_temp_ = init_temp;
}

void PropertiesManager::SetBackingProperties(double length, double height, double init_temp, int /*material*/) {
//    if (materials_.count(material) == 0) {
//        qDebug() << "Unknown backing material!";
//        return;
//    }
//    backing_material_ = materials_[material];
    backing_length_ = length;
    backing_height_ = height;
    backing_init_temp_ = init_temp;
}

void PropertiesManager::SetToolProperties(double radius, double height, double penetration_depth, double init_temp,
                                         double angular_velo, double friction_coef, double f_z, double f_x, int /*material*/) {
//    if (materials_.count(material) == 0) {
//        qDebug() << "Unknown tool material!";
//        return;
//    }
//    tool_material_ = materials_[material];
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
    time_layers_ = time_layers_count;
}

void PropertiesManager::SetHeatExchangePropeties(double alpha_1, double alpha_2, double alpha_3, double alpha_4, double alpha_4_tool, double out_temp) {
    alpha1_ = alpha_1;
    alpha2_ = alpha_2;
    alpha3_ = alpha_3;
    alpha4_ = alpha_4;
    alpha4_tool_ = alpha_4_tool;
    out_temp_ = out_temp;
}

void PropertiesManager::PrintAllProperties() {
    qDebug() << "--- Method ---";
    qDebug() << "Delta t = " << delta_t_;
    qDebug() << "Eps1 = " << eps1_;
    qDebug() << "Eps2 = " << eps2_;
    qDebug() << "Max iter = " << max_iter_;
    qDebug() << "Time layers = " << time_layers_;
    qDebug() << "Tool layers = " << tool_words_;
    qDebug() << "";
    qDebug() << "--- Plate ---";
    qDebug() << "Height = " << plate_height_;
    qDebug() << "Length = " << plate_lenght_;
    qDebug() << "Init temp = " << plate_init_temp_;
    qDebug() << "Density = " << plate_material_.GetDensity().ApproximateAt(plate_init_temp_);
    qDebug() << "Heat capacity = " << plate_material_.GetHeatCapacity().ApproximateAt(plate_init_temp_);
    qDebug() << "Thermal cond = " << plate_material_.GetThermalConductivity().ApproximateAt(plate_init_temp_);
    qDebug() << "";
    qDebug() << "--- Backing ---";
    qDebug() << "Height = " << backing_height_;
    qDebug() << "Length = " << backing_length_;
    qDebug() << "Init temp = " << backing_init_temp_;
    qDebug() << "Density = " << backing_material_.GetDensity().ApproximateAt(backing_init_temp_);
    qDebug() << "Heat capacity = " << backing_material_.GetHeatCapacity().ApproximateAt(backing_init_temp_);
    qDebug() << "Thermal cond = " << backing_material_.GetThermalConductivity().ApproximateAt(backing_init_temp_);
    qDebug() << "";
    qDebug() << "--- Tool ---";
    qDebug() << "Height = " << tool_height_;
    qDebug() << "Length = " << tool_radius_;
    qDebug() << "Penetration = " << tool_penetration_depth_;
    qDebug() << "Angular Velo = " << tool_angular_velo_;
    qDebug() << "Friction = " << friction_coef_;
    qDebug() << "Fx = " << f_x_;
    qDebug() << "Fz = " << f_z_;
    qDebug() << "Init temp = " << tool_init_temp_;
    qDebug() << "Density = " << tool_material_.GetDensity().ApproximateAt(tool_init_temp_);
    qDebug() << "Heat capacity = " << tool_material_.GetHeatCapacity().ApproximateAt(tool_init_temp_);
    qDebug() << "Thermal cond = " << tool_material_.GetThermalConductivity().ApproximateAt(tool_init_temp_);
    qDebug() << "";
    qDebug() << "--- Heat Exchange ---";
    qDebug() << "alpha1 = " << alpha1_;
    qDebug() << "alpha2 = " << alpha2_;
    qDebug() << "alpha3 = " << alpha3_;
    qDebug() << "alpha4 = " << alpha4_;
    qDebug() << "alpha4_tool = " << alpha4_tool_;
    qDebug() << "Out temp = " << out_temp_;
}

// Alpha getters
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

// Physical properties getters
double PropertiesManager::GetDensity(int x, int z, const Matrix& temperatures) {
    Material* material = materials_[materials_grid_[x][z]];
    double temperature = temperatures[x][z];
    return material->GetDensity().ApproximateAt(temperature);
}

double PropertiesManager::GetHeatCapacity(int x, int z, const Matrix& temperatures) {
    Material* material = materials_[materials_grid_[x][z]];
    double temperature = temperatures[x][z];
    return material->GetHeatCapacity().ApproximateAt(temperature);
}

double PropertiesManager::GetThermalConductivity(double x, double z, const Matrix& temperatures) {
    int x1 = std::floor(x);
    int x2 = std::ceil(x);

    int z1 = std::floor(z);
    int z2 = std::ceil(z);

    double int_part;
    if (std::modf(x, &int_part) == 0.0) {
        x1 = static_cast<int>(int_part);
        x2 = x1;
    } else if (std::modf(z, &int_part) == 0.0) {
        z1 = static_cast<int>(int_part);
        z2 = z1;
    } else {
        qDebug() << "Oh, something wrong in indexes for thermal condactivity: " << x << z;
    }

    Material* prev_material = materials_[materials_grid_[x1][z1]];
    double prev_temp = temperatures[x1][z1];
    double prev_value = prev_material->GetThermalConductivity().ApproximateAt(prev_temp);

    Material* next_material = materials_[materials_grid_[x2][z2]];
    double next_temp = temperatures[x2][z2];
    double next_value = next_material->GetThermalConductivity().ApproximateAt(next_temp);

    return 2 * prev_value * next_value / (prev_value + next_value);
}

// Delta getters
double PropertiesManager::GetDeltaT() {
    return delta_t_;
}

double PropertiesManager::GetDeltaX(int i) {
    return delta_x_[i];
}

double PropertiesManager::GetDeltaBackX(int i) {
    double current = delta_x_[i];
    double next = i + 1 < delta_x_.GetSize() ? delta_x_[i + 1] : current;
    return (current + next) / 2;
}

double PropertiesManager::GetDeltaZ(int i) {
    return delta_z_[i];
}

double PropertiesManager::GetDeltaBackZ(int i) {
    double current = delta_z_[i];
    double next = i + 1 < delta_z_.GetSize() ? delta_z_[i + 1] : current;
    return (current + next) / 2;
}

// Tool getters
double PropertiesManager::GetToolHeight() {
    return tool_height_;
}

double PropertiesManager::GetToolPenetration() {
    return tool_penetration_depth_;
}

double PropertiesManager::GetToolWaveHeight() {
    return tool_wave_height_;
}

double PropertiesManager::GetToolInitTemperature() {
    return tool_init_temp_;
}

int PropertiesManager::GetToolStartI() {
    return i_tool_start_;
}

int PropertiesManager::GetToolFinishI() {
    return i_tool_finish_;
}

// Heat getters
double PropertiesManager::GetHeatOutput1() {
    return 0;
}

double PropertiesManager::GetHeatOutput2() {
    return 0;
}

double PropertiesManager::GetHeatOutput3() {
    return 0;
}

double PropertiesManager::GetHeatX(int x, int z) {
    // Ignoring node in tool
    if (z < i_tool_bottom_start_) {
        return 0;
    }

    if (x == i_tool_start_) {
        return GetHeatOutput1() / delta_x_[x];
    }
    if (x == i_tool_finish_ + 1) {
        return GetHeatOutput2() / delta_x_[x];
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

    return GetHeatOutput3() / delta_z_[z];
}

// Method getters
double PropertiesManager::GetEpsilon1() {
    return eps1_;
}

double PropertiesManager::GetEpsilon2() {
    return eps2_;
}

int PropertiesManager::GetMaxIterations() {
    return max_iter_;
}
