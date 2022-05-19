#include "main_window.h"

#include <QBoxLayout>
#include <QFileDialog>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPushButton>

#include "solver/solver.h"
#include "grid_data_processor.h"
#include "menu.h"
#include "properties_widget.h"
#include "solution_runner.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), properties_manager_() {
    QHBoxLayout* layout = new QHBoxLayout();

    heatmap_ = new Heatmap(kWidth + 2, kHeight + 2, this);
    heatmap_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    layout->addWidget(heatmap_);

    QVBoxLayout* properties_layout = new QVBoxLayout();
    properties_ = new PropertiesWidget(properties_manager_, this);
    properties_layout->addWidget(properties_);

    compute_button_ = new QPushButton("Вычислить", this);
    connect(compute_button_, SIGNAL(released()), this, SLOT(SaveProperties()));
    properties_layout->addWidget(compute_button_);

    layout->addLayout(properties_layout);

    menu_ = new Menu(this);

    main_widget_ = new QWidget(this);

    main_widget_->setLayout(layout);
    setCentralWidget(main_widget_);

    setGeometry(400, 200, 840, 480);

    GridDataProcessor* processor = new GridDataProcessor(heatmap_, 402, 202, 99, 100);
    processor->Start();
}

void MainWindow::SaveProperties() {
    properties_->ConfigManager();

    Solver solver(kWidth, kHeight, &properties_manager_, [&](const Matrix& matrix) {
        heatmap_->SetValues(matrix);
    });

    compute_button_->setEnabled(false);
    solver.Start();
    compute_button_->setEnabled(true);
}
