#include "properties_widget.h"

#include <QDebug>
#include <QFile>

PropertiesWidget::PropertiesWidget(QWidget* parent) : QTreeView(parent) {
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

int PropertiesWidget::GetTimeLayersCount() {
    return method_properties_->GetValue("ksloi_fin");
}

double PropertiesWidget::GetTimeStep() {
    return method_properties_->GetValue("tau");
}

double PropertiesWidget::GetOutTemp() {
    return heat_exch_properties_->GetValue("TeOut");
}

void PropertiesWidget::SaveToFile(const QString& file_name) {
    QFile file(file_name);
    if(!file.open(QIODevice::WriteOnly)) {
        qDebug() << "[SaveToFile] error opening file(" << file_name << "):" << file.error();
        return;
    }
    QTextStream file_stream(&file);

    file_stream << "plate\n";
    int plate_material = static_cast<int>(plate_properties_->GetValue("mat_plast"));
    SaveMaterialPropertiesToFile(file_stream, plate_material);
    plate_properties_->SaveToFile(file_stream);

    file_stream << "backing\n";
    int backing_material = static_cast<int>(backing_properties_->GetValue("mat_sub"));
    SaveMaterialPropertiesToFile(file_stream, backing_material);
    backing_properties_->SaveToFile(file_stream);

    file_stream << "tool\n";
    int tool_material = static_cast<int>(tool_properties_->GetValue("mat_tool"));
    SaveMaterialPropertiesToFile(file_stream, tool_material);
    tool_properties_->SaveToFile(file_stream);

    file_stream << "method\n";
    method_properties_->SaveToFile(file_stream);

    file_stream << "heat_exchange\n";
    heat_exch_properties_->SaveToFile(file_stream);

    file.close();
}

void PropertiesWidget::SaveMaterialPropertiesToFile(QTextStream& file_stream, int id) {
    QFile file(kMaterialsConfigPath + QString::number(id) + ".dat");
    if(!file.open(QIODevice::ReadOnly)) {
        qDebug() << "error opening file for material " << id;
        return;
    }

    auto material_properties = file.readAll();
    file_stream << material_properties;

    file.close();
}
