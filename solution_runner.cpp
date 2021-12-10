#include "solution_runner.h"

#include <QDebug>
#include <QFile>
#include <QProcess>

const QString SolutionRunner::kSolverFilePath = "://resources/solver/";
const QString SolutionRunner::kSolverInputFilePath = "./solver_input/";
const QString SolutionRunner::kSolverFileName = "solver.exe";
const QString SolutionRunner::kRunnerFileName = "runner.sh";

void SolutionRunner::Run() {
    for (const QString& file_name: {kSolverFileName, kRunnerFileName}) {
        QFile::copy(kSolverFilePath + file_name, kSolverInputFilePath + file_name);
        QFile file(kSolverInputFilePath + file_name);
        file.setPermissions(QFile::ReadOwner | QFile::WriteOwner | QFile::ExeOwner);
    }

    QProcess* process = new QProcess();
    process->setWorkingDirectory(kSolverInputFilePath);
    process->start("./" + kRunnerFileName);

    process->waitForFinished();

    if (process->exitCode() != 0) {
        qDebug() << "Error during running solver...\n";
        qDebug() << "STDOUT:\n" << process->readAllStandardOutput() << "\n";
        qDebug() << "STDERR:\n" << process->readAllStandardError() << "\n";
    }
    delete process;
}
