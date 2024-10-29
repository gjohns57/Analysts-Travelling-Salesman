cd igraph
rm -r build
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="../.."
cmake --build .
cmake --build . --target check
cmake --install .