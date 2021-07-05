# Go into the directory where this bash script is contained.
cd `dirname $0`

# remove the build directory
rm -rf build

# make a new build directory and navigate into it
mkdir build && cd build

# compile the cpp code
cmake .. && make

# run
./particle_filter