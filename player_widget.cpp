#include "player_widget.h"

#include <QBoxLayout>

#include "grid_data_processor.h"

PlayerWidget::PlayerWidget(GridDataProcessor* processor, QWidget* parent) : QWidget(parent) {
    QHBoxLayout* layout = new QHBoxLayout();

    play_button_ = new QPushButton("Play", this);
    layout->addWidget(play_button_);
    connect(play_button_, SIGNAL(released()), processor, SLOT(Start()));

    stop_button_ = new QPushButton("Stop", this);
    layout->addWidget(stop_button_);
    connect(stop_button_, SIGNAL(released()), processor, SLOT(Pause()));

    time_line_ = new QSlider(Qt::Horizontal, this);
    layout->addWidget(time_line_);

    playback_rate_ = new QComboBox(this);
    playback_rate_->addItems(kPlaybackRates);
    layout->addWidget(playback_rate_);

    setLayout(layout);
}

void PlayerWidget::ChangeTimeLineState() {}
void PlayerWidget::ChangePlaybackRate() {}
