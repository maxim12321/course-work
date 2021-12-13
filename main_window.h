#pragma once

#include <QMainWindow>
#include "menu.h"
#include "properties_widget.h"
#include "grid_data_processor.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

    ~MainWindow() = default;

private slots:
    void SaveProperties();
    void ComputeTask();

private:
    Menu* menu_;
    PropertiesWidget* properties_;
    GridDataProcessor* processor_;
};
