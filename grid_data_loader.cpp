#include "grid_data_loader.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>

GridDataLoader::GridDataLoader(int width, int height)
    : width_(width), height_(height) {}

QVector<QVector<qreal>> GridDataLoader::LoadGridData() const {
    QFile file(kGridFileName);
    if(!file.open(QIODevice::ReadOnly)) {
        throw std::runtime_error("error opening file: " + file.error());
    }

    QTextStream file_stream(&file);
    QVector<QVector<qreal>> data(height_, QVector<qreal>(width_));

    for (int y = 0; y < height_; y++) {
        for (int x = 0; x < width_; x++) {
            file_stream >> data[y][x];
        }
    }
    return data;
}

QVector<QPointF> GridDataLoader::LoadPlotData() const {
    QFile file(kPlotFileName);
    if(!file.open(QIODevice::ReadOnly)) {
        throw std::runtime_error("error opening file: " + file.error());
    }

    QList<QPointF> data;

    while (!file.atEnd()) {
        QString params_line = file.readLine().trimmed();
        QStringList params = params_line.split("      ");

        double x = params[0].trimmed().toDouble();
        double temperature = params[1].trimmed().toDouble();

        data.append(QPointF(x, temperature));
    }

    file.close();
    return data;
}
