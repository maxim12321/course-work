#pragma once

#include <QTreeView>

class PropertiesWidget : public QTreeView {
    Q_OBJECT

public:
    explicit PropertiesWidget(QWidget* parent = nullptr);

    ~PropertiesWidget() = default;

private:
    const QStringList kPropertyFileNames = {"tool_properties.dat"};
};
