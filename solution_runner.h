#pragma once

#include <QString>

class SolutionRunner {
public:
    static void Run();

private:
    static const QString kSolverFilePath;
    static const QString kSolverInputFilePath;

    static const QString kSolverFileName;
    static const QString kRunnerFileName;
};
