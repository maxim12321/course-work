#pragma once

#include <QMap>
#include <QStandardItem>
#include <QTextStream>

class PropertiesItem : public QStandardItem {
public:
    explicit PropertiesItem(const QString& config_file_name);

    void SaveToFile(QTextStream& file_stream);

    double GetValue(const QString& name);
    QList<double> GetValues();

private:
    void CreateFromFile();

private:
    const QString kConfigFilePath = "://resources/properties/";
private:
    QString config_file_name_;
    QString solver_input_file_name_;
    QMap<QString, QStandardItem*> property_label_to_item_;
};
