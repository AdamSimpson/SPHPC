#ifndef SPH_SRC_ALLOC_H_
#define SPH_SRC_ALLOC_H_

#include "stdio.h"
#include "stdlib.h"

#define SAFE_ALLOC(num_elements, element_size) \
        SafeAlloc(num_elements, element_size, __FILE__, __LINE__)

static inline void* SafeAlloc(const size_t num_elements,
                              const size_t element_size,
                              const char* file_name,
                              const size_t line_number) {
  void *ptr = calloc(num_elements, element_size);
  if (!ptr) {
    printf("Failed to allocate memory in %s at %zu\n", file_name, line_number);
    exit(-1);
  }
  else
    return ptr;
}

#endif
