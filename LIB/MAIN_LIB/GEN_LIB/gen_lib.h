#ifndef gen_lib_H
#define gen_lib_H

#include "mainlib.h"

//--------------- DATA DIRECTORY + FILE CREATION -----------------// 

// Creation of directory with DATA from simulations //
// dirName = Directory name
// fileName = File name
// num_file = file number fp[num_file] to allow multiple writing with fopen
// mode = 0: destruction + creation + write => mode of fopen "w"
// mode = 1: open + write => mode of fopen "a"

extern FILE * fp[20];
void DIR_FILE(char dirName[], const char *fileName, int num_file, int mode);

void createNDimensionalArray(void ***array_name, size_t n, const size_t *size, size_t type_size);
void freeNDimensionalArray(void *array_name, size_t n, const size_t *size);


#endif /* gen_lib_H */