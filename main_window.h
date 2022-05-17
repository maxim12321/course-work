#pragma once

#include <QMainWindow>
#include "menu.h"
#include "properties_widget.h"
#include "heatmap.h"
#include "solver/properties_manager.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

    ~MainWindow() = default;

private slots:
    void SaveProperties();

private:
    constexpr static int kWidth = 10;
    constexpr static int kHeight = 10;

private:
    QWidget* main_widget_;
    Heatmap* heatmap_;
    Menu* menu_;
    PropertiesWidget* properties_;
    PropertiesManager properties_manager_;
};
