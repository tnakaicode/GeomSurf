SetFactory("OpenCASCADE");

a() = ShapeFromFile("wing001.step");

Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 1.0;
