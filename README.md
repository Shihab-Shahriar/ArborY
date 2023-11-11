Minimum seperation distance functions between many 3d geometrical shapes, 0.0 distance implies overlap.

## Compilation
+ There are two executables: oen called `app` that demonstrates the use of these functions. The other is called `cgal`, mainly for testing purposes. It calls the many of the same functions using CGAL library to compare output against `app`.
+ To compile: from a `build` directory, run `bash ../compile.sh` and then `make`. Comment out the last few lines if CGAL is not installed.

## Usage
+ Most of the functions are in `include/distances.h`. 
+ `vector.h` is a quick implementation of a 3d vector class.
+ I expect many of these functions to be rewritten based on the geometric representation mundy/stk uses for these primitive shapes.


## Notes
+ Testing is done against CGAL for all except capsule pairs. But I only chose few ad-hoc pairs of test cases. In short, testing isn't comprehensive, but it's robust for the cases I've tested.
+ 