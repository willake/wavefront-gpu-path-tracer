typedef struct {
  float3 aabbMin, aabbMax;  // 24 bytes
  uint leftFirst, triCount; // 8 bytes; total: 32 bytes
} BVHNode;