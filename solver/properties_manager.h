#pragma once

#include <QObject>

class PropertiesManager {
public:
    PropertiesManager() = default;

    // Returns either lower, uppper or tool heat capacity
    qreal GetHeatCapacity(qreal x, qreal z);

private:
    qreal lower_heat_capacity_;
    qreal upper_heat_capacity_;
    qreal tool_heat_capacity_;

    // Other properties here...
};
