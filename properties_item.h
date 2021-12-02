#pragma once

#include <QStandardItem>

class PropertiesItem : public QStandardItem {
public:
    explicit PropertiesItem(const QString& file_name);



private:
    void CreateFromFile();

private:
    const QString kFilePath = "://resources/properties/";

private:
    QString file_name_;
};
