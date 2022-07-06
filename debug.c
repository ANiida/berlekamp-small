/* Obtain a backtrace and print it to stdout. */
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>

void print_trace(void)   //// 呼ばれない
{
    void *array[10];
    int size = backtrace(array, 10);
    char ** strings = backtrace_symbols(array, size);

    if (strings != NULL) {
        printf("Obtained %d stack frames.\n", size);
        for (int i = 0; i < size; i++)
            printf("%s\n", strings[i]);
    }
    free(strings);
}

