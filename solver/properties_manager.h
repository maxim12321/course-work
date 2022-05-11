#pragma once

#include "utils/matrix.h"
#include "utils/vector.h"

#include <QObject>
#include <QMap>
#include <utility>

class PropertiesManager {
    struct Material {
        double density;
        double heat_capacity;
        double thermal_conductivity;
    };

    enum class Object {
        kBacking,
        kPlate,
        kTool,
    };

public:
    PropertiesManager();

    void SetPlateProperties(double length, double height, double init_temp, int material);
    void SetBackingProperties(double length, double height, double init_temp, int material);
    void SetToolProperties(double radius, double height, double penetration_depth, double init_temp,
                           double angular_velo, double friction_coef, double f_z, double f_x, int material);
    void SetMethodProperties(double delta_t, double eps1, double eps2, int max_iter_count, int time_layers_count, int tool_words);

    void SetHeatExchangePropeties(double alpha_1, double alpha_2, double alpha_3, double alpha_4, double alpha_4_tool, double out_temp);
    double GetAlpha1();
    double GetAlpha2();
    double GetAlpha3();
    double GetAlpha4();
    double GetAlpha4Tool();
    double GetOutTemperature();

    Matrix InitializeGrids(int nx, int nz);

    double GetDeltaT();

    double GetDeltaX(int i);
    double GetDeltaZ(int i);

    double GetInitTemperature(int x, int z);

    double GetDensity(double x, double z);
    double GetHeatCapacity(double x, double z);
    double GetThermalConductivity(double x, double z);

private:
    std::pair<double, double> GridToObjectPosition(int x, int z);
    Object GetObjectInPositon(double x, double z);

    void LoadMaterials();

private:
    const QString kMaterialsConfigFile = "://resources/properties/materials";
    QMap<int, Material> materials_;

    // Plate properties
    double plate_lenght_;
    double plate_height_;
    double plate_init_temp_;
    Material plate_material_;

    // Backing properties
    double backing_length_;
    double backing_height_;
    double backing_init_temp_;
    Material backing_material_;

    double total_height_;
    double height_without_penetration_;
    double total_length_;
    double tool_start_;
    double tool_finish_;

    // Tool properties
    double tool_radius_;
    double tool_height_;
    double tool_penetration_depth_;
    double tool_init_temp_;
    double tool_angular_velo_;
    Material tool_material_;
    double friction_coef_;
    double f_z_;
    double f_x_;

    // Heat exchange properties
    double alpha1_;
    double alpha2_;
    double alpha3_;
    double alpha4_;
    double alpha4_tool_;
    double out_temp_;

    // Method properties
    double delta_t_;
    double eps1_;
    double eps2_;
    int max_iter_;
    int time_layers_;
    int tool_words_;

    // Grids for method
    Matrix density_grid_;
    Matrix heat_capacity_grid_;
    Matrix thermal_conductivity_grid_;
    Matrix init_temp_;
    Vector delta_x_;
    Vector delta_z_;
};
