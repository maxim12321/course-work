#include "grid_data_loader.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <cassert>

int GridDataLoader::LoadData() {
    QFile file(kOutputFileName);
    if(!file.open(QIODevice::ReadOnly)) {
        throw std::runtime_error("error opening file: " + file.error());
    }

    QTextStream file_stream(&file);

    int size = file_stream.readLine().toInt();

    data_.clear();

    int height = 0;
    int width = 0;
    QVector<QVector<double>> matrix;
    QVector<double> vals;
    for (int i = 0; i < size; i++) {
        QString line_sizes = file_stream.readLine();
        if (line_sizes.isEmpty()) {
            i--;
            continue;
        }

        QStringList sizes_list = line_sizes.trimmed().split(' ');
        assert(sizes_list.size() == 2);
        int h = sizes_list[0].toInt();
        int w = sizes_list[1].toInt();
        assert((height == 0 || h == height) && (width == 0 || w == width));
        height = h;
        width = w;
        matrix.resize(height);
        vals.resize(width);

        for (int row = 0; row < height; row++) {
            QString s = file_stream.readLine().trimmed();
            assert(s.front() == '[' && s.back() == ']');

            s.remove(0, 1);
            if (row == 0) {
                assert(s.front() == '[');
                s.remove(0, 1);
            }

            s.chop(1);
            if (row + 1 == height) {
                assert(s.back() == ']');
                s.chop(1);
            }

            QStringList vals_str = s.split(' ');
            assert(vals_str.size() == width);
            for (int column = 0; column < width; column++) {
                vals[column] = vals_str[column].toDouble();
            }
            matrix[row] = vals;
        }

        data_.push_back(matrix);
    }
    height_ = height;
    width_ = width;

    file.close();
    return size;
}

int GridDataLoader::GetGridWidth() {
    return width_;
}

int GridDataLoader::GetGridHeight() {
    return height_;
}

const QVector<QVector<double>>& GridDataLoader::GetGridData(int layer) {
    assert(layer >= 0 && layer < data_.size());
    return data_[layer];
}
