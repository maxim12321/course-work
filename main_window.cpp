#include "main_window.h"

#include <QBoxLayout>
#include <QPushButton>

#include "properties_widget.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    layout->addWidget(new PropertiesWidget(this));

    QWidget* main_widget = new QWidget(this);
    main_widget->setLayout(layout);
    setCentralWidget(main_widget);

    setGeometry(400, 200, 640, 480);
}
