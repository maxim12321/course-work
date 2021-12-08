#include "main_window.h"

#include <QBoxLayout>
#include <QFileDialog>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPushButton>
#include <QPushButton>

#include "grid_value_rect_item.h"
#include "menu.h"
#include "properties_widget.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* main_layout = new QHBoxLayout();

    QVBoxLayout* graphics_layout = new QVBoxLayout();

    QGraphicsScene* scene = new QGraphicsScene(0, 0, 400, 400, this);
    scene->addItem(new GridValueRectItem(0, 0, 400, 400));
    graphics_layout->addWidget(new QGraphicsView(scene, this));

    player_ = new PlayerWidget(this);
    graphics_layout->addWidget(player_);

    main_layout->addItem(graphics_layout);

    properties_ = new PropertiesWidget(this);
    main_layout->addWidget(properties_);

    menu_ = new Menu(this);

    QWidget* main_widget = new QWidget(this);

    main_widget->setLayout(main_layout);
    setCentralWidget(main_widget);

    setGeometry(400, 200, 840, 480);
}

void MainWindow::SaveProperties() {
    QString filename = QFileDialog::getSaveFileName(this, "Properties saving");
    if (filename.size() == 0) {
        return;
    }
    properties_->SaveToFile(filename);
}
