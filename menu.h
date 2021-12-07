#pragma once

#include <QMenuBar>
#include <QMainWindow>
#include <QApplication>

class MainWindow;

class Menu : public QMainWindow {
public:
    explicit Menu(MainWindow* parent = nullptr);

    ~Menu() = default;

private:
};
