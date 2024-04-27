/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file gpu.h
 * @brief Contain gpu constants
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#define MAX_GPU_STRANDS    1000
#define LAST_STRAND_INDEX  MAX_GPU_STRANDS - 2 // because on 0 is length 
#define MAX_GPU_STRAND_LEN 1000
#define GPU_DEPTH          5
#define GPU_BACK_DEPTH     2
#define GPU_FRONT_DEPTH    4

#define GPU_MEM_SIZE       (MAX_GPU_STRANDS * MAX_GPU_STRAND_LEN * GPU_DEPTH + 1)
#define GPU_BACK_MEM_SIZE  (MAX_GPU_STRANDS * MAX_GPU_STRAND_LEN * GPU_BACK_DEPTH + 1)
#define GPU_FRONT_MEM_SIZE  (MAX_GPU_STRANDS * MAX_GPU_STRAND_LEN * GPU_FRONT_DEPTH + 1)

// definition of depth
#define GPU_PARTNER_STRAND 0
#define GPU_PARTNER_ATOM   1 
#define GPU_DOMAIN         2
#define GPU_NUCLEOTIDES    3
#define GPU_NOTE           4

// on first position is alwais length
#define GPU_ELEMENT(x, y, z)    (1 + y + (x + 1) * MAX_GPU_STRAND_LEN + z * MAX_GPU_STRAND_LEN * MAX_GPU_STRANDS)
#define GET_GPU_COUNT_COL(y, z) GPU_ELEMENT(-1, y, z)
#define GET_GPU_COUNT_ROW()     0
