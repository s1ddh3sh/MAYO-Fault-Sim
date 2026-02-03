// #include "mayo_example.h"
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include "mayo_sim3.c"

int main(void)
{
    // const mayo_params_t *params = &MAYO_1;
    int n = 100;
    // unsigned char msg[32] = "Hello";
    // size_t msglen = 32;
    // memset(msg + strlen((char *)msg), 0, msglen - strlen((char *)msg));
    // return example_mayo(params, msg, msglen);
    while (n-- && !example_fault_sim(NULL))
        ;
}