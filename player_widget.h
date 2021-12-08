#pragma once

#include <QWidget>
#include <QPushButton>
#include <QSlider>
#include <QComboBox>

class GridDataProcessor;

class PlayerWidget : public QWidget {
    Q_OBJECT

public:
    PlayerWidget(GridDataProcessor* processor, QWidget* parent = nullptr);

    ~PlayerWidget() = default;

signals:
    void PlaybackRateChanged(double rate);

private slots:
    void ChangeTimeLineState();
    void ChangePlaybackRate();

private:
    QPushButton* play_button_;
    QPushButton* stop_button_;
    QSlider* time_line_;
    QComboBox* playback_rate_;

    const QStringList kPlaybackRates = {
        "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2"
    };
};
