struct Tri {
  float3 vertex0, vertex1, vertex2; // 48 bytes
  float3 normal0, normal1, normal2; // 48 bytes
  float2 uv0, uv1, uv2;             // 24 bytes
  float3 centroid;                  // 16 bytes
  int objIdx;                       // 4 bytes
  // 140 in total
};