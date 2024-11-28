#include "gen_lib.h"

FILE *fp[20];

//--------------- DATA DIRECTORY + FILE CREATION -----------------// 

// Creation of directory with DATA from simulations //
// dirName = Directory name
// fileName = File name
// num_file = file number fp[num_file] to allow multiple writing with fopen
// mode = 0: destruction + creation + write => mode of fopen "w"
// mode = 1: open + write => mode of fopen "a"

void DIR_FILE(char dirName[], const char *fileName, int num_file, int mode)
{
	char file[100];
	
    if (mkdir(dirName, 0777) == -1) 
	{
        if (errno != EEXIST) {
            perror("Error creating directory");
            exit(EXIT_FAILURE);
        }
    }
	
	sprintf(file, "%s/%s",dirName,fileName);
	
	if (mode == 0)
	{
		fp[num_file] = fopen(file, "w" );
	}
	else if (mode == 1)
	{
		fp[num_file] = fopen(file, "a" );
	}
	
    if (fp[num_file] == NULL)
    {
        perror("Erreur lors de l'ouverture du fichier");
        exit(EXIT_FAILURE);
    }
	
}

// Create an N-dimensional array of size size[i] for each dimension i, 
// the type of elements is indicated by type_size and the name of the array is given by array_name

void createNDimensionalArray(void ***array_name, size_t n, const size_t *size, size_t type_size)
{
	assert(n > 0);
	assert(type_size > 0);
	assert(array_name != NULL);
	assert(size != NULL);
	
	if(n == 1)
	{
		*array_name = calloc(size[0], type_size);
		return;
	}
	
	*array_name = calloc(size[0], sizeof(void*));
	assert(*array_name != NULL);

	for (size_t i = 0; i < size[0]; i++)
	{
		createNDimensionalArray((void***)(*array_name + i), n-1, size+1, type_size);
	}
}

void freeNDimensionalArray(void *array_name, size_t n, const size_t *size)
{
	assert(array_name != NULL);
	assert(n > 0);
	assert(size != NULL);
	
	if (n == 1)
	{
		free(array_name);
		return;
	}
	
	for (size_t i = 0; i < size[0]; i++)
	{
		freeNDimensionalArray(*((void**)array_name + i), n-1, size+1);
	}
	
	free(array_name);
}
