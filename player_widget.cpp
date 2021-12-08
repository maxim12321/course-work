#include "player_widget.h"

#include <QBoxLayout>

PlayerWidget::PlayerWidget(QWidget* parent) : QWidget(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    play_button_ = new QPushButton("Play", this);
    layout->addWidget(play_button_);

    stop_button_ = new QPushButton("Stop", this);
    layout->addWidget(stop_button_);

    time_line_ = new QSlider(Qt::Horizontal, this);
    layout->addWidget(time_line_);

    playback_rate_ = new QComboBox(this);
    playback_rate_->addItems(kPlaybackRates);
    layout->addWidget(playback_rate_);

    setLayout(layout);
}

void PlayerWidget::ChangeTimeLineState() {}
void PlayerWidget::ChangePlaybackRate() {}
