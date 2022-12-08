#pragma once

#include <QString>

#include "properties_widget.h"

class SolutionRunner {
public:
    static void Run(PropertiesWidget* properties);

private:
    static void PrepareSolutionInput(PropertiesWidget* properties);
private:
    static const QString kSolverFilePath;
    static const QString kSolverInputFilePath;

    static const QString kSolverFileName;
    static const QString kRunnerFileName;
};
