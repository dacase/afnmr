#include <stdio.h>
#include <stdlib.h>

// Amber's portable getline() for systems that lack it

ssize_t getline(char **lineptr, size_t *n, FILE *stream);
