#pragma once

#include <QMainWindow>
#include "menu.h"
#include "properties_widget.h"
#include "grid_data_processor.h"
#include "heatmap.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

    ~MainWindow() = default;

private slots:
    void SaveProperties();

private:
    constexpr static int kWidth = 100;
    constexpr static int kHeight = 100;

private:
    QWidget* main_widget_;
    QPushButton* compute_button_;
    Heatmap* heatmap_;
    GridDataProcessor* processor_;
    Menu* menu_;
    PropertiesWidget* properties_;
};
