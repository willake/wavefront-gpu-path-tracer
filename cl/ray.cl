// Define ray
// __attribute__((aligned(64)))
typedef struct __attribute__((aligned(128))) {
  float3 O, D, rD;    // 36 bytes
  float t;            // 4 bytes
  float2 barycentric; // 8 bytes
  int objIdx;         // 4 bytes
  int triIdx;         // 4 bytes
  int traversed;
  int tested;
  bool inside; // 1 bytes
} Ray;         // total 64 bytes