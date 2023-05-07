#pragma once

#include <QString>
#include <QVector>

class GridDataLoader {
public:
    GridDataLoader() = default;

    void LoadData(int size, double out_temp);

    int GetGridWidth();
    int GetGridHeight();
    const QVector<QVector<double>>& GetGridData(int layer);

private:
    const QString kOutputFileName = "./solver_input/output.txt";

private:
    QVector<QVector<QVector<double>>> data_;
    int width_;
    int height_;
};
