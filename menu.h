#pragma once

#include <QMenuBar>
#include <QMainWindow>
#include <QApplication>

class Menu : public QMainWindow {
public:
    explicit Menu(QMainWindow* parent = nullptr);

    ~Menu() = default;

private:
};
