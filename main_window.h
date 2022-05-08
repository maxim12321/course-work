#pragma once

#include <QMainWindow>
#include "menu.h"
#include "properties_widget.h"
#include "heatmap.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

    ~MainWindow() = default;

private slots:
    void SaveProperties();

private:
    QWidget* main_widget_;
    Heatmap* heatmap_;
    Menu* menu_;
    PropertiesWidget* properties_;
};
