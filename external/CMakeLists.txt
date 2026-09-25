include(delaunator.cmake)
include(dkm.cmake)
include(FastNoiseLite.cmake)
include(hmm.cmake)
include(libnpy.cmake)
include(nn-c.cmake)
include(mixbox.cmake)

include(noise.cmake)
include(terrain-descriptors.cmake)

add_subdirectory(CLWrapper)
add_subdirectory(PointSampler)

include(nanoflann.cmake)

include(lightusd.cmake)
