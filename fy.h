#ifndef FY_H
#define FY_H
/// function prototype in fy.c
void merge_rand(unsigned short a[], int n);

#define SWAP(type, a, b)   { type temp = a; a = b; b = temp; }

#define MIN_TO_WORK     {         \
        d_type temp;              \
        if (*array_i <= *array_j) \
            temp = *(array_i++);  \
        else                      \
            temp = *(array_j++);  \
        *(work_k++) = temp;       \
    }

#define MAX_TO_WORK     {         \
        d_type temp;              \
        if (*array_x > *array_y)  \
            temp = *(array_x--);  \
        else                      \
            temp = *(array_y--);  \
        *(work_z--) = temp;       \
    }

#endif /// FY_H
