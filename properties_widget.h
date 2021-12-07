#pragma once

#include <QTreeView>
#include "properties_item.h"

class PropertiesWidget : public QTreeView {
    Q_OBJECT

public:
    explicit PropertiesWidget(QWidget* parent = nullptr);

    ~PropertiesWidget() = default;

    void SaveToFile(const QString& file_name);
    void LoadFromFile(const QString& file_name);

private:
    const QList<QPair<QString, QString>> kPropertyConfigFileNames = {
        {"tool", "tool_properties.dat"}
    };

    QMap<QString, PropertiesItem*> property_to_item_;
};
