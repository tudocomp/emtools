# Sample STXXL Project

This is a sample project using the STXXL library via CMake project files.

This method should be used for "finished" projects, where the source code is
shipped without including the STXXL library. Thus when building, a binary
installation must be detected. Using CMake project files, the STXXL library can
be detected via the standard `find_package()` mechanism.

Note that for **prototype** projects and **early development**, STXXL can be
used and included via simpler methods. See the STXXL documentation for details.

See http://stxxl.sourceforge.net

Written 2013-10-26 Timo Bingmann
