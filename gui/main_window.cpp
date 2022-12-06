#include "main_window.h"

#include <QBoxLayout>
#include <QFileDialog>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPushButton>

#include "grid_data_processor.h"
#include "menu.h"
#include "properties_widget.h"
#include "solution_runner.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    heatmap_ = new Heatmap(kWidth + 2, kHeight + 2, this);
    heatmap_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    layout->addWidget(heatmap_);

    QVBoxLayout* properties_layout = new QVBoxLayout();
    properties_ = new PropertiesWidget(this);
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

    processor_ = new GridDataProcessor(heatmap_);
}

void MainWindow::SaveProperties() {
    compute_button_->setEnabled(false);
    SolutionRunner::Run(properties_);

    int layers = properties_->GetTimeLayersCount();
    double step_s = properties_->GetTimeStep();
    int step_ms = static_cast<int>(step_s * 1000);
    processor_->ProcessSolverOutput(layers, step_ms);
    processor_->Start();

    compute_button_->setEnabled(true);
}
