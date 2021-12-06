#include "grid_value_rect_item.h"

#include <QtMath>
#include <QVector2D>

#include <QDebug>

GridValueRectItem::GridValueRectItem(int x, int y, int width, int height,
                                     QGraphicsItem* parent)
    : QGraphicsRectItem(x, y, width, height, parent),
      grid_values_(2, QVector<qreal>(2, 0)),
      grid_width_(2),
      grid_height_(2) {
    // random colors to test
    grid_values_[0][1] = kMaxValue;
    grid_values_[1][0] = kMaxValue / 2;
}

void GridValueRectItem::paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) {
    painter->save();

    for (int x = 0; x < rect().width(); x++) {
        for (int y = 0; y < rect().height(); y++) {
            QPoint grid_nearest = FindNearestOnGrid(QPoint(x, y));
            QVector<QColor> colors;

            for (int add_x = 0; add_x <= 1; add_x++) {
                for (int add_y = 0; add_y <= 1; add_y++) {
                    QPoint grid_point = grid_nearest + QPoint(add_x, add_y);
                    if (grid_point.x() >= grid_width_ || grid_point.y() >= grid_height_) {
                        continue;
                    }

                    QPointF distance = DistanceFromGrid(grid_point, QPoint(x, y));

                    // TODO: maybe adjust this formula somehow
                    QColor color = GetColorByGridValue(grid_point);
                    color.setAlphaF((1 - distance.x()) * (1 - distance.y()));
                    colors.push_back(color);
                }
            }

            painter->setPen(BlendColors(colors));
            painter->drawPoint(QPoint(x, y));
        }
    }

    painter->restore();
}

void GridValueRectItem::SetGridValues(const QVector<QVector<qreal>>& grid_values) {
    if (grid_values.size() < 2 || grid_values[0].size() < 2) {
        throw std::runtime_error("Grid size must be >= 2");
    }
    grid_values_ = grid_values;
    grid_height_ = grid_values.size();
    grid_width_ = grid_values[0].size();
}

QColor GridValueRectItem::BlendColors(const QVector<QColor>& colors) {
    QVector2D total(0, 0);
    for (const QColor& color: colors) {
        total += QVector2D(1, 0) * qCos(color.hueF() * 2 * M_PI) * color.alphaF();
        total += QVector2D(0, 1) * qSin(color.hueF() * 2 * M_PI) * color.alphaF();
    }
    total.normalize();

    qreal hue = qAtan2(total.y(), total.x());
    hue /= M_PI * 2;
    if (hue < 0) {
        hue += 1;
    }
    return QColor::fromHsvF(hue, 1, 1);
}

QPoint GridValueRectItem::FindNearestOnGrid(const QPoint& to) {
    if (grid_values_.size() < 2 || grid_values_[0].size() < 2) {
        throw std::runtime_error("Grid size is less than 2");
    }
    qreal cell_height = rect().height() / (grid_values_.size() - 1);
    qreal cell_width = rect().width() / (grid_values_[0].size() - 1);

    return QPoint(qFloor(to.x() / cell_width), qFloor(to.y() / cell_height));
}

QPointF GridValueRectItem::DistanceFromGrid(const QPoint& grid_point, const QPoint& to) {
    if (grid_values_.size() < 2 || grid_values_[0].size() < 2) {
        throw std::runtime_error("Grid size is less than 2");
    }
    qreal cell_height = rect().height() / (grid_values_.size() - 1);
    qreal cell_width = rect().width() / (grid_values_[0].size() - 1);

    qreal grid_x = grid_point.x() * cell_width;
    qreal grid_y = grid_point.y() * cell_height;

    return QPointF(qAbs(grid_x - to.x()) / cell_width,
                   qAbs(grid_y - to.y()) / cell_height);
}

QColor GridValueRectItem::GetColorByGridValue(const QPoint& grid_point) {
    qreal value = grid_values_[grid_point.y()][grid_point.x()];
    qreal color_coefficient = (value - kMinValue) / (kMaxValue - kMinValue);

    QColor min_color = kMinValueColor;
    min_color.setAlphaF(1 - color_coefficient);

    QColor max_color = kMaxValueColor;
    max_color.setAlphaF(color_coefficient);

    return BlendColors({min_color, max_color});
}
