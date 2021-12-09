#include "main_window.h"

#include <QBoxLayout>
#include <QFileDialog>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPushButton>

#include "grid_data_processor.h"
#include "grid_value_rect_item.h"
#include "menu.h"
#include "properties_widget.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    QGraphicsScene* scene = new QGraphicsScene(0, 0, 100, 100, this);
    GridValueRectItem* grid_item = new GridValueRectItem(0, 0, 100, 100);
    scene->addItem(grid_item);

    QGraphicsView* view = new QGraphicsView(scene, this);
    view->scale(4, 4);
    layout->addWidget(view);

    QVBoxLayout* properties_layout = new QVBoxLayout();
    properties_ = new PropertiesWidget(this);
    properties_layout->addWidget(properties_);

    QPushButton* compute_button = new QPushButton("Вычислить", this);
    connect(compute_button, SIGNAL (released()), properties_, SLOT (CreateInputForSolver()));
    properties_layout->addWidget(compute_button);

    layout->addLayout(properties_layout);

    menu_ = new Menu(this);

    QWidget* main_widget = new QWidget(this);

    main_widget->setLayout(layout);
    setCentralWidget(main_widget);

    setGeometry(400, 200, 840, 480);

    GridDataProcessor* processor = new GridDataProcessor(grid_item, 99, 100);
    processor->Start();
}

void MainWindow::SaveProperties() {
    QString filename = QFileDialog::getSaveFileName(this, "Properties saving");
    if (filename.size() == 0) {
        return;
    }
    properties_->SaveToFile(filename);
}
