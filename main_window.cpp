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
#include "solution_runner.h"
#include "plots_widget.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    QVBoxLayout* visual_layout = new QVBoxLayout();

    QGraphicsScene* scene = new QGraphicsScene(0, 0, 512, 512, this);
    GridValueRectItem* grid_item = new GridValueRectItem(0, 0, 512, 512);
    scene->addItem(grid_item);
    QGraphicsView* view = new QGraphicsView(scene, this);

    view->scale(2, 2);
    visual_layout->addWidget(view);

    PlotsWidget* plots = new PlotsWidget();
//    visual_layout->addWidget(plots);

    layout->addLayout(visual_layout);

    QVBoxLayout* properties_layout = new QVBoxLayout();
    properties_ = new PropertiesWidget(this);
    properties_layout->addWidget(properties_);

    QPushButton* compute_button = new QPushButton("Вычислить", this);
    connect(compute_button, SIGNAL (released()), this, SLOT (ComputeTask()));
    properties_layout->addWidget(compute_button);

    layout->addLayout(properties_layout);

    menu_ = new Menu(this);

    QWidget* main_widget = new QWidget(this);

    main_widget->setLayout(layout);
    setCentralWidget(main_widget);

    setGeometry(400, 200, 840, 480);

    processor_ = new GridDataProcessor(grid_item, plots, 402, 202, 99, 100);
}

void MainWindow::SaveProperties() {
    QString filename = QFileDialog::getSaveFileName(this, "Properties saving");
    if (filename.size() == 0) {
        return;
    }
    properties_->SaveToFile(filename);
}

void MainWindow::ComputeTask() {
    properties_->CreateInputForSolver();
    SolutionRunner::Run();
    processor_->Start();
}
