#pragma once

#include "matrix.h"

#include <string>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <fstream>
#include <atomic>
#include <thread>

class ResultSaver {
public:
    ResultSaver(const std::string &file_name);
    ~ResultSaver();

    void Save(const Matrix &matrix);

private:
    void WorkerRoutine();

private:
    std::thread worker_;

    std::atomic<bool> is_finished_;
    std::mutex mutex_;
    std::condition_variable empty_;
    std::queue<Matrix> results_;

    std::ofstream output_;
};