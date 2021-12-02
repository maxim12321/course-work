#include "main_window.h"

#include <QBoxLayout>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPushButton>

#include "grid_value_rect_item.h"
#include "properties_widget.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    QGraphicsScene* scene = new QGraphicsScene(0, 0, 400, 400, this);
    scene->addItem(new GridValueRectItem(0, 0, 400, 400));
    layout->addWidget(new QGraphicsView(scene, this));

    layout->addWidget(new PropertiesWidget(this));

    QWidget* main_widget = new QWidget(this);
    main_widget->setLayout(layout);
    setCentralWidget(main_widget);

    setGeometry(400, 200, 640, 480);
}
