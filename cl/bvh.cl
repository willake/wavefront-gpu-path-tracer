typedef struct {
  float minx, miny, minz;
  float maxx, maxy, maxz;
  uint leftFirst, triCount; // total: 32 bytes
} BVHNode;

typedef struct {
  float minx, miny, minz;
  float maxx, maxy, maxz;
  uint leftRight; // 2x16 bits
  uint BLAS;
} TLASNode;

typedef struct {
  float v0x, v0y, v0z;
  float v1x, v1y, v1z;
  float v2x, v2y, v2z;
  float cx, cy, cz;
} Tri;

typedef struct {
  float2 uv0, uv1, uv2;
  float N0x, N0y, N0z;
  float N1x, N1y, N1z;
  float N2x, N2y, N2z;
  float dummy;
} TriEx;

typedef struct {
  uint dummy1, dummy2;
  uint idx;
  float16 transform;
  float16 invTransform; // inverse transform
  uint dummy[6];
} BVHInstance;