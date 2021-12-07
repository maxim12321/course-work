#pragma once

#include <QTreeView>
#include "properties_item.h"

class PropertiesWidget : public QTreeView {
    Q_OBJECT

public:
    explicit PropertiesWidget(QWidget* parent = nullptr);

    ~PropertiesWidget() = default;

    void SaveToFile(const QString& file_name);

private:
    const QList<QPair<QString, QString>> kPropertyConfigFileNames = {
        {"tool", "tool_properties.dat"},
        {"plate", "plate_material_properties.dat"},
        {"backing", "backing_material_properties.dat"},
        {"heat_exchange", "heat_exchange_properties.dat"},
        {"method", "numerical_method_properties.dat"},
    };

    QMap<QString, PropertiesItem*> property_to_item_;
};
