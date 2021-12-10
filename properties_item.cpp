#include "properties_item.h"

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QStandardItem>
#include <QString>

PropertiesItem::PropertiesItem(const QString& config_file_name, const QString& solver_input_file_name)
    : QStandardItem(), config_file_name_(config_file_name), solver_input_file_name_(solver_input_file_name) {
    setEditable(false);

    CreateFromFile();
}

void PropertiesItem::CreateFromFile() {
    QFile file(kConfigFilePath + config_file_name_);
    if(!file.open(QIODevice::ReadOnly)) {
        qDebug() << "error opening file: " << file.error();
        return;
    }

    setText(file.readLine().trimmed());

    int max_name_width = 0;
    while (!file.atEnd()) {
        QString params_line = file.readLine().trimmed();
        QStringList params = params_line.split('|');

        if (params.empty() || (params.size() == 1 && params[0].trimmed().isEmpty())) {
            continue;
        }
        if (params.size() != 3) {
            qDebug() << "Wrong number of parameters in file\n";
            continue;
        }

        QList<QStandardItem*> property_items;

        QStandardItem* property_name = new QStandardItem(params[1].trimmed());
        property_name->setEditable(false);
        property_items.append(property_name);

        max_name_width = std::max(max_name_width, params[1].trimmed().size());

        QStandardItem* property_value = new QStandardItem(params[2].trimmed());
        property_label_to_item_[params[0].trimmed()] = property_value;
        property_value->setEditable(true);
        property_items.append(property_value);

        appendRow(property_items);
    }

    setSizeHint(QSize(max_name_width * 10, 10));
}

void PropertiesItem::CreateInputForSolver() {
    if (!QDir(kInputFilePath).exists()) {
        QDir().mkdir(kInputFilePath);
    }
    QFile file(kInputFilePath + solver_input_file_name_);
    if(!file.open(QIODevice::WriteOnly)) {
        qDebug() << "error opening file: " << file.error();
        return;
    }
    QTextStream file_stream(&file);

    for (int i = 0; i < rowCount(); ++i) {
        file_stream << child(i, 1)->text() << "\n";
    }

    file.close();
}

void PropertiesItem::SaveToFile(QTextStream& file_stream) {
    for (const auto& label : property_label_to_item_.keys()) {
        file_stream << label << " | " << property_label_to_item_[label]->text() << "\n";
    }
}
