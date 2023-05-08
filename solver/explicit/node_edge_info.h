#pragma once

enum class EdgeType {
  kAir,
  kMaterial,
  kMixed,
  kNone
};

struct NodeEdgeInfo {
  double width = 0;
  double height = 0;

  EdgeType left_edge = EdgeType::kNone;
  EdgeType right_edge = EdgeType::kNone;
  EdgeType top_edge = EdgeType::kNone;
  EdgeType bottom_edge = EdgeType::kNone;

  int i = -1;
  int k = -1;
};
