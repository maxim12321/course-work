#pragma once

#include <QMainWindow>
#include "menu.h"
#include "properties_widget.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

    ~MainWindow() = default;

private slots:
    void SaveProperties();

private:
    Menu* menu_;
    PropertiesWidget* properties_;
};
