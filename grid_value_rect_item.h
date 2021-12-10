#pragma once

#include <QColor>
#include <QGraphicsRectItem>
#include <QPainter>
#include <QVector>

class GridValueRectItem : public QGraphicsRectItem {
public:
    GridValueRectItem(int x, int y, int width, int height,
                      QGraphicsItem* parent = nullptr);

    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;

    void SetGridValues(const QVector<QVector<qreal>>& grid_values);

private:
    QColor BlendColors(const QVector<QColor>& colors);

    QPoint FindNearestOnGrid(const QPoint& to);

    QPointF DistanceFromGrid(const QPoint& grid_point, const QPoint& to);

    QColor GetColorByGridValue(const QPoint& grid_point);

private:
    const QColor kMinValueColor = QColor::fromHsv(230, 255, 255);
    const QColor kMaxValueColor = QColor::fromHsv(0, 255, 255);

private:
    QVector<QVector<qreal>> grid_values_;
    int grid_width_;
    int grid_height_;

    qreal min_value_ = 0;
    qreal max_value_ = 1;
};
