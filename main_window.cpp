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

    heatmap_ = new Heatmap(402, 202, this);
    heatmap_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    layout->addWidget(heatmap_);

    QVBoxLayout* properties_layout = new QVBoxLayout();
    properties_ = new PropertiesWidget(this);
    properties_layout->addWidget(properties_);

    QPushButton* compute_button = new QPushButton("Вычислить", this);
    connect(compute_button, SIGNAL (released()), properties_, SLOT (CreateInputForSolver()));
    properties_layout->addWidget(compute_button);

    layout->addLayout(properties_layout);

    menu_ = new Menu(this);

    main_widget_ = new QWidget(this);

    main_widget_->setLayout(layout);
    setCentralWidget(main_widget_);

    setGeometry(400, 200, 840, 480);

    SolutionRunner::Run();

    GridDataProcessor* processor = new GridDataProcessor(heatmap_, 402, 202, 99, 100);
    processor->Start();
}

void MainWindow::SaveProperties() {
    QString filename = QFileDialog::getSaveFileName(this, "Properties saving");
    if (filename.size() == 0) {
        return;
    }
    properties_->SaveToFile(filename);
}
