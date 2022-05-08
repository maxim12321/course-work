#pragma once

#include <QObject>
#include <QTimer>

#include "grid_data_loader.h"
#include "heatmap.h"

class GridDataProcessor : public QObject {
    Q_OBJECT

public:
    GridDataProcessor(Heatmap* heatmap,
                      int width, int height,
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
    Heatmap* heatmap_;
    int total_steps_;
    int default_timer_interval_;

    QTimer* timer_;
    int current_step_ = 0;
    qreal time_scale_ = 1;

    GridDataLoader data_loader_;
    QVector<QVector<qreal>> grid_data_;
};
