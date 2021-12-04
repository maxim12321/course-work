#include "main_window.h"

#include <QBoxLayout>
#include <QGraphicsScene>
#include <QGraphicsView>

#include "grid_value_rect_item.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QGraphicsScene* scene = new QGraphicsScene(0, 0, 400, 400, this);
    scene->addItem(new GridValueRectItem(0, 0, 400, 400));

    QGraphicsView* view = new QGraphicsView(scene, this);

    QHBoxLayout* layout = new QHBoxLayout();
    layout->addWidget(view);

    QWidget* main_widget = new QWidget(this);
    main_widget->setLayout(layout);
    setCentralWidget(main_widget);

    setGeometry(400, 200, 640, 480);
}
