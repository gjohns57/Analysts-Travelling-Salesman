WORKSPACE=$(pwd)
cd build
cmake .. -DCMAKE_PREFIX_PATH="$WORKSPACE"
cmake --build .


for file in *; do
    if [ -f "$file" ] && [ "$file" != "cmake_install.cmake" ] && [ "$file" != "CMakeCache.txt" ] && [ "$file" != "Makefile" ]; then
        cp "$file" ../bin/
    fi
done