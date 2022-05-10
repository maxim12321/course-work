#include "properties_widget.h"
#include "solver/properties_manager.h"

#include <QDebug>
#include <QFile>

PropertiesWidget::PropertiesWidget(PropertiesManager& manager, QWidget* parent) : QTreeView(parent), manager_(manager) {
    QStandardItemModel* model = new QStandardItemModel(0, 1);

    plate_properties_ = new PropertiesItem(kPlateConfigFileName);
    backing_properties_ = new PropertiesItem(kBackingConfigFileName);
    tool_properties_ = new PropertiesItem(kToolConfigFileName);
    heat_exch_properties_ = new PropertiesItem(kHeatExchConfigFileName);
    method_properties_ = new PropertiesItem(kMethodConfigFileName);

    model->appendRow(plate_properties_);
    model->appendRow(backing_properties_);
    model->appendRow(tool_properties_);
    model->appendRow(heat_exch_properties_);
    model->appendRow(method_properties_);

    model->setHorizontalHeaderItem(0, new QStandardItem("Параметр"));
    model->setHorizontalHeaderItem(1, new QStandardItem("Значение"));

    setModel(model);
    resizeColumnToContents(0);
}

void PropertiesWidget::SaveToFile(const QString& file_name) {
    QFile file(file_name);
    if(!file.open(QIODevice::WriteOnly)) {
        qDebug() << "error opening file: " << file.error();
        return;
    }
    QTextStream file_stream(&file);

    file_stream << "plate\n";
    plate_properties_->SaveToFile(file_stream);

    file_stream << "backing\n";
    backing_properties_->SaveToFile(file_stream);

    file_stream << "tool\n";
    tool_properties_->SaveToFile(file_stream);

    file_stream << "method\n";
    method_properties_->SaveToFile(file_stream);

    file_stream << "heat exchange\n";
    heat_exch_properties_->SaveToFile(file_stream);

    file.close();
}

void PropertiesWidget::ConfigManager() {
    auto plate_props = plate_properties_->GetValues();
    manager_.SetPlateProperties(plate_props[0], plate_props[1], plate_props[3], plate_props[2]);

    auto backing_props = backing_properties_->GetValues();
    manager_.SetBackingProperties(backing_props[0], backing_props[1], backing_props[3], backing_props[2]);

    auto tool_props = tool_properties_->GetValues();
    manager_.SetToolProperties(tool_props[0], tool_props[1], tool_props[2], tool_props[4], tool_props[5], tool_props[6], tool_props[7], tool_props[8], tool_props[3]);

    auto method_props = method_properties_->GetValues();
    manager_.SetMethodProperties(method_props[0], method_props[3], method_props[4], method_props[5], method_props[1], method_props[2]);

    auto heat_exch_props = heat_exch_properties_->GetValues();
    manager_.SetHeatExchangePropeties(heat_exch_props[0], heat_exch_props[1], heat_exch_props[2], heat_exch_props[3], heat_exch_props[4], heat_exch_props[5]);
}
