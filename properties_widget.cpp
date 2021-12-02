#include "properties_widget.h"

#include "properties_item.h"

PropertiesWidget::PropertiesWidget(QWidget* parent) : QTreeView(parent) {
    QStandardItemModel* model = new QStandardItemModel(0, 1);

    for (const QString& property_file_name: kPropertyFileNames) {
        model->appendRow(new PropertiesItem(property_file_name));
    }

    model->setHorizontalHeaderItem(0, new QStandardItem("Параметр"));
    model->setHorizontalHeaderItem(1, new QStandardItem("Значение"));

    setModel(model);
    resizeColumnToContents(0);
}
