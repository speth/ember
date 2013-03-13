#include "mathUtils.h"

inline VecMap map_vector(double* data, index_t size, index_t stride=1)
{
    return VecMap(data, size, Stride1X(1, stride));
}

inline MatrixMap map_matrix(double* data, index_t rows, index_t cols,
                            index_t inner, index_t outer)
{
    return MatrixMap(data, rows, cols, StrideXX(inner, outer));
}
