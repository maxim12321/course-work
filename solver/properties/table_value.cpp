#include "table_value.h"

void TableValue::AddTableValue(int key, double value) { table_[key] = value; }

double TableValue::ApproximateAt(double key) const {
  auto next = table_.upper_bound(key);

  // key is less than first table value
  if (next == table_.begin()) {
    return next->second;
  }

  // key is bigger than last table value
  if (next == table_.end()) {
    return table_.rbegin()->second;
  }

  auto prev = next;
  prev--;
  double next_coef = (key - prev->first) / (next->first - prev->first);

  return next->second * next_coef + prev->second * (1 - next_coef);
}
