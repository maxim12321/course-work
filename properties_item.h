#pragma once

#include <QMap>
#include <QStandardItem>

class PropertiesItem : public QStandardItem {
public:
    explicit PropertiesItem(const QString& file_name);

    void SaveToFile(QTextStream& file_stream);

private:
    void CreateFromFile();

private:
    const QString kFilePath = "://resources/properties/";

private:
    QString file_name_;
    QMap<QString, QStandardItem*> property_label_to_item_;
};
