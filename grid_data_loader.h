#pragma once

#include <QString>
#include <QVector>

class GridDataLoader {
public:
    GridDataLoader(int width, int height);

    QVector<QVector<qreal>> LoadData() const;

private:
    const QString kOutputFileName = "./solver_input/Te_C.dat";

private:
    int width_;
    int height_;
};
