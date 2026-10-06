#!/bin/bash

# directories to be formatted (recursive search)
DIRS="HighMap/include HighMap/src examples tests benchmarks"
FORMAT_CMD="clang-format -style=file:scripts/clang_style -i"

CPP_ONLY=0

usage() {
    echo "Usage: $0 [--cpp-only] [-h|--help]"
    echo "  --cpp-only   format only *.cpp files"
}

for ARG in "$@"; do
    case ${ARG} in
        --cpp-only) CPP_ONLY=1 ;;
        -h|--help)  usage; exit 0 ;;
        *)          echo "Unknown option: ${ARG}"; usage; exit 1 ;;
    esac
done

echo "- clang-format"

if [ ${CPP_ONLY} -eq 1 ]; then
    # format only *.cpp files
    for D in ${DIRS}; do
        for F in `find ${D}/. -type f -iname \*.cpp`; do
            echo ${F}
            ${FORMAT_CMD} ${F}
        done
    done
    exit 0
fi

# format opencl kernels
for D in ${DIRS}; do
    for F in `find ${D}/. -type f -iname \*.cl`; do
        echo ${F}
        sed '1d;$d' ${F} > ${F}_tmp
        ${FORMAT_CMD} ${F}_tmp
        sed -i '1s/^/R""(\n/' ${F}_tmp
        echo ')""' >> ${F}_tmp
        mv ${F}_tmp ${F}
    done
done

# format C++
for D in ${DIRS}; do
    for F in `find ${D}/. -type f \( -iname \*.hpp -o -iname \*.inl -o -iname \*.cpp \)`; do
        echo ${F}
        ${FORMAT_CMD} ${F}
    done
done

# fix headers
for D in ${DIRS}; do
    for F in `find ${D}/. -type f \( -iname \*.hpp -o -iname \*.inl \)`; do
        echo ${F}
        scripts/fix_format_doxygen.py ${F}
    done
done

# format cmake files
echo "- cmake-format"

cmake-format -i CMakeLists.txt HighMap/CMakeLists.txt external/CMakeLists.txt external/*.cmake
