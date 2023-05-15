#include "result_saver.h"

#include <iostream>
#include <stdexcept>

ResultSaver::ResultSaver(const std::string &file_name) 
    : is_finished_(false), mutex_(), output_(file_name, std::ios::out) {
  if (output_.fail()) {
    std::cerr << "An error occurred while opening the output file" << std::endl;
    throw std::runtime_error("failed to open output file");
  }

  worker_ = std::thread([&]() {
    WorkerRoutine();
  });
}

ResultSaver::~ResultSaver() {
    is_finished_.store(true);
    empty_.notify_all();
    worker_.join();
}

void ResultSaver::Save(const Matrix &matrix) {
    std::lock_guard guard{mutex_};
    results_.push(matrix);
    empty_.notify_all();
}

void ResultSaver::WorkerRoutine() {
    Matrix result;
    while (true) {
        {
            std::unique_lock lock{mutex_};
            while (!is_finished_.load() && results_.empty()) {
                empty_.wait(lock);
            }

            if (results_.empty()) {
                return;
            }
            
            result = std::move(results_.front());
            results_.pop();
        }
        output_ << result << "\n";
    }
}