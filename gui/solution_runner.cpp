#include "solution_runner.h"

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QProcess>

const QString SolutionRunner::kSolverFilePath = "://resources/solver/";
const QString SolutionRunner::kSolverInputFilePath = "./solver_input/";
const QString SolutionRunner::kSolverFileName = "solver";
const QString SolutionRunner::kRunnerFileName = "runner.sh";

void SolutionRunner::Run(PropertiesWidget* properties) {
    PrepareSolutionInput(properties);
    for (const QString& file_name: {kSolverFileName, kRunnerFileName}) {
        if (QFile::exists(kSolverInputFilePath + file_name)) {
            QFile::remove(kSolverInputFilePath + file_name);
        }
        QFile::copy(kSolverFilePath + file_name, kSolverInputFilePath + file_name);
        QFile file(kSolverInputFilePath + file_name);
        file.setPermissions(QFile::ReadOwner | QFile::WriteOwner | QFile::ExeOwner);
    }

    QProcess* process = new QProcess();
    process->setWorkingDirectory(kSolverInputFilePath);
    process->start("./" + kRunnerFileName);

    process->waitForFinished(-1);

    qDebug() << "STDOUT:\n" << process->readAllStandardOutput() << "\n";
    if (process->exitCode() != 0) {
        qDebug() << "Error during running solver...\n";
        qDebug() << "STDERR:\n" << process->readAllStandardError() << "\n";
    }
    qDebug() << "Solver is finished";

    delete process;
}

void SolutionRunner::PrepareSolutionInput(PropertiesWidget* properties) {
    QDir input_dir(kSolverInputFilePath);
    if (!input_dir.exists()) {
        input_dir.mkpath(".");
    }

    properties->SaveToFile(kSolverInputFilePath + "input.txt");
}
