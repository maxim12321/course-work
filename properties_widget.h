#pragma once

#include <QTreeView>
#include "properties_item.h"

class PropertiesWidget : public QTreeView {
    Q_OBJECT

public:
    explicit PropertiesWidget(QWidget* parent = nullptr);

    ~PropertiesWidget() = default;

    void SaveToFile(const QString& file_name);

public slots:
    void CreateInputForSolver();

private:
    const QList<QStringList> kPropertyConfigFileNames = {
        {"tool", "tool_properties.dat", "inp_Tool.dat"},
        {"plate", "plate_material_properties.dat", "inp_Plast.dat"},
        {"backing", "backing_material_properties.dat", "inp_Sub.dat"},
        {"heat_exchange", "heat_exchange_properties.dat", "inp_CoeffToBound.dat"},
        {"method", "numerical_method_properties.dat", "inp_method.dat"},
    };

    QMap<QString, PropertiesItem*> property_to_item_;
};
