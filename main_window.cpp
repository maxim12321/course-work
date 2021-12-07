#include "main_window.h"

#include <QBoxLayout>
#include <QFileDialog>
#include <QPushButton>
#include "menu.h"
#include "properties_widget.h"

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    properties_ = new PropertiesWidget(this);
    layout->addWidget(properties_);

    menu_ = new Menu(this);

    QWidget* main_widget = new QWidget(this);

    main_widget->setLayout(layout);
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
