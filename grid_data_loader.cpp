#include "grid_data_loader.h"

#include <QFile>
#include <QTextStream>

GridDataLoader::GridDataLoader(int width, int height)
    : width_(width), height_(height) {}

QVector<QVector<qreal>> GridDataLoader::LoadData() const {
    QFile file(kOutputFileName);
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
