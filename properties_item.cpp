#include "properties_item.h"

#include <QDebug>
#include <QFile>
#include <QStandardItem>
#include <QString>

PropertiesItem::PropertiesItem(const QString& file_name)
    : QStandardItem(), file_name_(file_name) {
    setEditable(false);

    CreateFromFile();
}

void PropertiesItem::CreateFromFile() {
    QFile file(kFilePath + file_name_);
    if(!file.open(QIODevice::ReadOnly)) {
        qDebug() << "error opening file: " << file.error();
        return;
    }

    setText(file.readLine().trimmed());

    int max_width = 0;
    while (!file.atEnd()) {
        QString params_line = file.readLine().trimmed();
        QStringList params = params_line.split('|');

        if (params.empty() || (params.size() == 1 && params[0].trimmed().isEmpty())) {
            continue;
        }
        if (params.size() > 2) {
            qDebug() << "More than two parameters in file\n";
            continue;
        }

        QList<QStandardItem*> property_items;

        QStandardItem* property_name = new QStandardItem(params[0].trimmed());
        property_name->setEditable(false);
        property_items.append(property_name);

        max_width = std::max(max_width, params[0].size());

        if (params.size() > 1) {
            QStandardItem* property_value = new QStandardItem(params[1].trimmed());
            property_value->setEditable(true);
            property_items.append(property_value);
        }

        appendRow(property_items);
    }

    setSizeHint(QSize(max_width * 10, 10));
}
