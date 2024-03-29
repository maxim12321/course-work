#pragma once

#include <QTreeView>
#include "properties_item.h"

class PropertiesManager;

class PropertiesWidget : public QTreeView {
    Q_OBJECT

public:
    explicit PropertiesWidget(QWidget* parent = nullptr);

    ~PropertiesWidget() = default;

    void SaveToFile(const QString& file_name);

    int GetTimeLayersCount();
    double GetTimeStep();

private:
    void SaveMaterialPropertiesToFile(QTextStream& file_stream, int material_id);

private:
    const QString kMaterialsConfigPath = "://resources/properties/materials/";

    const QString kPlateConfigFileName = "plate_material_properties.dat";
    PropertiesItem* plate_properties_;

    const QString kBackingConfigFileName = "backing_material_properties.dat";
    PropertiesItem* backing_properties_;

    const QString kToolConfigFileName = "tool_properties.dat";
    PropertiesItem* tool_properties_;

    const QString kHeatExchConfigFileName = "heat_exchange_properties.dat";
    PropertiesItem* heat_exch_properties_;

    const QString kMethodConfigFileName = "numerical_method_properties.dat";
    PropertiesItem* method_properties_;
};
