#pragma once

#include <QObject>
#include <QTimer>

#include "grid_value_rect_item.h"

class GridDataProcessor : public QObject {
    Q_OBJECT

public:
    GridDataProcessor(GridValueRectItem* grid_item,
                      int total_steps, int timer_interval = 1000 / 30);

signals:
    // 0 <= time_passed <= 1 -- current progress
    void DataUpdated(qreal time_passed);

public slots:
    void Start();

    void Resume();
    void Pause();

    void SetTimeScale(qreal time_scale);

    // 0 <= time_passed <= 1
    void SetCurrentProgress(qreal time_passed);

private slots:
    void UpdateData();

private:
    QVector<QVector<qreal>> GetGridData(int step);
    void DisplayData();
    void RerunTimer();

private:
    GridValueRectItem* grid_item_;
    int total_steps_;
    int default_timer_interval_;

    QTimer* timer_;
    int current_step_ = 0;
    qreal time_scale_ = 1;
};
