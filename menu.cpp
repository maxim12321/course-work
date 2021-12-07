#include "menu.h"

#include <QMenu>
#include <QMenuBar>

#include "main_window.h"

Menu::Menu(MainWindow *parent) : QMainWindow(parent) {
    QAction *save_properties = new QAction("&Сохранить параметры", parent);

    QMenu *file;
    file = parent->menuBar()->addMenu("&Файл");
    file->addAction(save_properties);

    connect(save_properties, SIGNAL(triggered()), parent, SLOT(SaveProperties()));
}
