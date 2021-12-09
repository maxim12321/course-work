#pragma once

#include <QMap>
#include <QStandardItem>
#include <QTextStream>

class PropertiesItem : public QStandardItem {
public:
    explicit PropertiesItem(const QString& config_file_name, const QString& solver_input_file_name);

    void SaveToFile(QTextStream& file_stream);

    void CreateInputForSolver();

private:
    void CreateFromFile();

private:
    const QString kConfigFilePath = "://resources/properties/";
    const QString kInputFilePath = "solver_input/";
private:
    QString config_file_name_;
    QString solver_input_file_name_;
    QMap<QString, QStandardItem*> property_label_to_item_;
};
