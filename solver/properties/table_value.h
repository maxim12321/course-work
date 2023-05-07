#pragma once

#include <map>
#include <vector>

class TableValue {
public:
  void AddTableValue(int key, double value);

  // Returns weighted sum of nearest table values,
  // or just nearest value (if *key* is outside the table)
  double ApproximateAt(double key) const;

private:
  std::map<int, double> table_;
};
