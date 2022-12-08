#include "grid_data_processor.h"

#include <QDebug>

GridDataProcessor::GridDataProcessor(Heatmap *heatmap) : heatmap_(heatmap) {
    timer_ = new QTimer(this);
    connect(timer_, SIGNAL(timeout()), this, SLOT(UpdateData()));
}

void GridDataProcessor::ProcessSolverOutput(int total_steps, int time_interval) {
    qDebug() << "ProcessSolverOutput";
    total_steps_ = total_steps;
    default_timer_interval_ = std::max(1, time_interval);
    data_loader_.LoadData(total_steps);
    heatmap_->Resize(data_loader_.GetGridWidth(), data_loader_.GetGridHeight());
}

void GridDataProcessor::Start() {
    if (timer_->isActive()) {
        timer_->stop();
    }
    current_step_ = 0;
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
        current_step_ = total_steps_ - 1;
        Pause();
    }

    DisplayData();

    emit DataUpdated(1. * current_step_ / total_steps_);
}

QVector<QVector<qreal>> GridDataProcessor::GetGridData(int) {
    return grid_data_;
}

void GridDataProcessor::DisplayData() {
    heatmap_->SetValues(data_loader_.GetGridData(current_step_));
}

void GridDataProcessor::RerunTimer() {
    if (timer_->isActive()) {
        Pause();
        Resume();
    }
}
