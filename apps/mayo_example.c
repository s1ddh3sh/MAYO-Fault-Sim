#include <mayo.h>
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mayo_example.h"

int example_mayo(const mayo_params_t *p, const unsigned char *msg, size_t msglen)
{
    int res = MAYO_OK;
    // size_t msglen = 32;
    size_t smlen = PARAM_sig_bytes(p) + msglen;

    unsigned char *pk = calloc(PARAM_cpk_bytes(p), 1);
    unsigned char *sk = calloc(PARAM_csk_bytes(p), 1);

    uint64_t *epk = calloc(1, sizeof(pk_t));
    sk_t *esk = calloc(1, sizeof(sk_t));

    unsigned char *sig = calloc(PARAM_sig_bytes(p) + msglen, 1);

    // unsigned char msg[32] = {0};
    unsigned char msg2[32] = {0};

    // printf("Enter msg : ");
    // fflush(stdout);
    // if (fgets((char *)msg, sizeof(msg), stdin) == NULL)
    // {
    //     printf("Error\n");
    //     goto err;
    // }
    // printf("\n");
    // size_t inp_len = strlen((char *)msg);
    // if (inp_len < msglen)
    // {
    //     memset(msg + inp_len, 0, msglen - inp_len);
    // }
    printf("Input message: \n");
    for (size_t i = 0; i < msglen; ++i)
    {
        printf("%02x", msg[i]);
    }
    printf("\n");
    printf("Example with %s\n", PARAM_name(p));

    printf("mayo_keypair -> ");
    res = mayo_keypair(p, pk, sk);
    if (res != MAYO_OK)
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        printf("OK\n");
    }

    // Print public key (hex)
    printf("Public key: \n");
    for (int i = 0; i < PARAM_cpk_bytes(p); ++i)
    {
        printf("%02x", pk[i]);
    }
    printf("\n");

    // Print secret key (hex, first 32 bytes for brevity)
    printf("Secret key (first 32 bytes): ");
    for (int i = 0; i < (PARAM_csk_bytes(p) < 32 ? PARAM_csk_bytes(p) : 32); ++i)
    {
        printf("%02x", sk[i]);
    }
    printf("%s\n", PARAM_csk_bytes(p) > 32 ? "..." : "");

    printf("mayo_expand_sk -> ");
    res = mayo_expand_sk(p, sk, esk);
    if (res != MAYO_OK)
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        printf("OK\n");
    }

    printf("mayo_expand_pk -> ");
    res = mayo_expand_pk(p, pk, epk);
    if (res != MAYO_OK)
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        printf("OK\n");
    }

    printf("mayo_sign -> ");
    res = mayo_sign(p, sig, &smlen, msg, msglen, sk);
    if (res != MAYO_OK)
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        printf("OK\n");
    }

    // Print signature (hex)
    printf("Signature: \n");
    for (int i = 0; i < PARAM_sig_bytes(p); ++i)
    {
        printf("%02x", sig[i]);
    }
    printf("\n");
    printf("mayo_open (with correct signature) -> ");
    res = mayo_open(p, msg2, &msglen, sig, smlen, pk);
    if (res != MAYO_OK || memcmp(msg, msg2, msglen))
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        res = MAYO_OK;
        printf("OK\n");
    }

    printf("mayo_verify (with correct signature) -> ");
    res = mayo_verify(p, msg, msglen, sig, pk);
    if (res != MAYO_OK)
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        res = MAYO_OK;
        printf("OK\n");
    }

    printf("mayo_open (with altered signature) -> ");
    sig[0] = ~sig[0];
    memset(msg2, 0, msglen);
    res = mayo_open(p, msg2, &msglen, sig, smlen, pk);
    if (res != MAYO_ERR || !memcmp(msg, msg2, msglen))
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        res = MAYO_OK;
        printf("OK\n");
    }

    printf("mayo_verify (with altered signature) -> ");
    res = mayo_verify(p, msg, msglen, sig, pk);
    if (res == MAYO_OK)
    {
        printf("FAIL\n");
        res = -1;
        goto err;
    }
    else
    {
        res = MAYO_OK;
        printf("OK\n");
    }

err:
    free(pk);
    free(epk);
    mayo_secure_free(sk, PARAM_csk_bytes(p));
    mayo_secure_free(esk, sizeof(sk_t));
    free(sig);
    return res;
}

// int main(void)
// {
//     const mayo_params_t *params = &MAYO_1;
//     unsigned char msg[32] = "Hello";
//     size_t msglen = 32;
//     memset(msg + strlen((char *)msg), 0, msglen - strlen((char *)msg));
//     return example_mayo(params, msg, msglen);
// }