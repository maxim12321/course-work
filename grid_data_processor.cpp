#include "grid_data_processor.h"

#include <QDebug>

GridDataProcessor::GridDataProcessor(GridValueRectItem* grid_item, PlotsWidget* plots,
                                     int width, int height,
                                     int total_steps, int timer_interval)
    : grid_item_(grid_item),
      plots_(plots),
      total_steps_(total_steps),
      default_timer_interval_(timer_interval),
      data_loader_(width, height) {
    timer_ = new QTimer(this);
    connect(timer_, SIGNAL(timeout()), this, SLOT(UpdateData()));
}

void GridDataProcessor::Start() {
    if (timer_->isActive()) {
        timer_->stop();
    }
    current_step_ = 0;
    plots_->MakePlot(data_loader_.LoadPlotData());
    grid_data_ = data_loader_.LoadGridData();
    DisplayData();

    timer_->start(default_timer_interval_ / time_scale_);
}

void GridDataProcessor::Resume() {
    if (timer_->isActive()) {
        throw std::runtime_error("Resuming processor, that's already running");
    }
    timer_->start(default_timer_interval_ / time_scale_);
}

void GridDataProcessor::Pause() {
    timer_->stop();
}

void GridDataProcessor::SetTimeScale(qreal time_scale) {
    if (time_scale < 0.001) {
        throw std::runtime_error("Time scale should be > 0");
    }

    time_scale_ = time_scale;
    RerunTimer();
}

void GridDataProcessor::SetCurrentProgress(qreal time_passed) {
    if (time_passed < 0 || time_passed > 1) {
        throw std::runtime_error("Time passed should be in range [0; 1]");
    }
    current_step_ = qRound(time_passed * total_steps_);
    DisplayData();
    RerunTimer();
}

void GridDataProcessor::UpdateData() {
    current_step_++;
    if (current_step_ >= total_steps_) {
        current_step_ = total_steps_;
        Pause();
    }

    DisplayData();

    emit DataUpdated(1. * current_step_ / total_steps_);
}

QVector<QVector<qreal>> GridDataProcessor::GetGridData(int) {
    return grid_data_;
}

void GridDataProcessor::DisplayData() {
    grid_item_->SetGridValues(grid_data_);
}

void GridDataProcessor::RerunTimer() {
    if (timer_->isActive()) {
        Pause();
        Resume();
    }
}
