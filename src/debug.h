#ifndef SPH_SRC_DEBUG_H_
#define SPH_SRC_DEBUG_H_

#define DEBUG

// Macro to print file name, line number, and function information along
// with standard printf formated string
#ifdef DEBUG
#define DEBUG_PRINT(str, args...) printf("DEBUG: %s:%d:%s(): " str, \
 __FILE__, __LINE__, __func__, ##args)
#else
#define DEBUG_PRINT(str, args...) do {} while (0)
#endif

#endif
