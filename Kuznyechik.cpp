#include <iostream>
#include <chrono>
#include <emmintrin.h>
#include <cstring>
#include "Kuznyechik.h"

uint8_t mul_table[256][256];
uint8_t lut_table[Kuznyechik::BLOCK_SIZE][256][Kuznyechik::BLOCK_SIZE];

Kuznyechik::Kuznyechik(const unsigned char *key) {
    memset(round_keys, 0, sizeof(round_keys));

    expand_key(key, key + BLOCK_SIZE);

    gen_mul_table();
    gen_lut_table();
}

void Kuznyechik::X2(const unsigned char *a, const unsigned char *b, unsigned char *c) {
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        c[i] = a[i] ^ b[i];
    }
}

void Kuznyechik::F(const uint8_t *in_key_1, const uint8_t *in_key_2,
                   uint8_t *out_key_1, uint8_t *out_key_2,
                   uint8_t *iter_const) {
    uint8_t tmp[BLOCK_SIZE];
    memcpy(out_key_2, in_key_1, KEY_SIZE);

    X2(in_key_1, iter_const, tmp);
    S(tmp, tmp);
    L(tmp, tmp);
    X2(tmp, in_key_2, out_key_1);
}

void Kuznyechik::expand_key(const uint8_t *key_1, const uint8_t *key_2) {
    uint8_t iter_C[32][BLOCK_SIZE];
    uint8_t iter_num[32][BLOCK_SIZE];

    for (int i = 0; i < 32; ++i) {
        memset(iter_num[i], 0, BLOCK_SIZE);
        iter_num[i][15] = i + 1;
    }

    for (int i = 0; i < 32; ++i) {
        L(iter_num[i], iter_C[i]);
    }

    uint8_t iter_1[64];
    uint8_t iter_2[64];
    uint8_t iter_3[64];
    uint8_t iter_4[64];

    memcpy(round_keys[0], key_1, 64);
    memcpy(round_keys[1], key_2, 64);
    memcpy(iter_1, key_1, 64);
    memcpy(iter_2, key_2, 64);

    for (int i = 0; i < 4; ++i) {
        F(iter_1, iter_2, iter_3, iter_4, iter_C[0 + 8 * i]);
        F(iter_3, iter_4, iter_1, iter_2, iter_C[1 + 8 * i]);
        F(iter_1, iter_2, iter_3, iter_4, iter_C[2 + 8 * i]);
        F(iter_3, iter_4, iter_1, iter_2, iter_C[3 + 8 * i]);
        F(iter_1, iter_2, iter_3, iter_4, iter_C[4 + 8 * i]);
        F(iter_3, iter_4, iter_1, iter_2, iter_C[5 + 8 * i]);
        F(iter_1, iter_2, iter_3, iter_4, iter_C[6 + 8 * i]);
        F(iter_3, iter_4, iter_1, iter_2, iter_C[7 + 8 * i]);
        memcpy(round_keys[2 * i + 2], iter_1, 64);
        memcpy(round_keys[2 * i + 3], iter_2, 64);
    }
}

uint8_t Kuznyechik::GF_mul(uint8_t a, uint8_t b) {
    uint8_t c = 0;
    uint8_t hi_bit;

    for (int i = 0; i < 8; ++i) {
        if (b & 1) {
            c ^= a;
        }

        hi_bit = a & 0x80;
        a <<= 1;

        if (hi_bit) {
            a ^= 0xc3;
        }

        b >>= 1;
    }

    return c;
}

Kuznyechik::Matrix Kuznyechik::sqr_matrix(const Matrix &mat) {
    Matrix res{};

    for (size_t i = 0; i < BLOCK_SIZE; ++i)
        for (size_t j = 0; j < BLOCK_SIZE; ++j)
            for (size_t k = 0; k < BLOCK_SIZE; ++k)
                res[i][j] ^= mul_table[mat[i][k]][mat[k][j]];

    return res;
}

void Kuznyechik::gen_mul_table() {
    for (unsigned i = 0; i < 256; ++i)
        for (unsigned j = 0; j < 256; ++j)
            mul_table[i][j] = GF_mul(i, j);
}

void Kuznyechik::gen_lut_table() {
    Matrix l_matrix;
    for (size_t i = 0; i < BLOCK_SIZE; ++i)
        for (size_t j = 0; j < BLOCK_SIZE; ++j)
            if (i == 0)
                l_matrix[i][j] = rcon[j];
            else if (i == j + 1)
                l_matrix[i][j] = 1;
            else
                l_matrix[i][j] = 0;

    for (unsigned i = 0; i < 4; ++i)
        l_matrix = sqr_matrix(l_matrix);

    for (size_t i = 0; i < BLOCK_SIZE; ++i)
        for (size_t j = 0; j < 256; ++j)
            for (size_t k = 0; k < BLOCK_SIZE; ++k)
                lut_table[i][j][k] = mul_table[sbox[j]][l_matrix[k][i]];
}

void Kuznyechik::X(unsigned char *a, const unsigned char *b) {
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        a[i] = a[i] ^ b[i];
    }
}

void Kuznyechik::S(const unsigned char *input, unsigned char *output) {
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        output[i] = sbox[input[i]];
    }
}

void Kuznyechik::reverse_S(const unsigned char *input, unsigned char *output) {
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        output[i] = inv_sbox[input[i]];
    }
}

void Kuznyechik::R(uint8_t *state) {
    uint8_t a_0 = 0;
    uint8_t tmp[BLOCK_SIZE];

    for (int i = BLOCK_SIZE - 1; i >= 0; --i) {
        if (i > 0) {
            tmp[i] = state[i - 1];
        }
        a_0 ^= GF_mul(state[i], rcon[i]);
    }

    tmp[0] = a_0;
    memcpy(state, tmp, BLOCK_SIZE);
}

void Kuznyechik::reverse_R(uint8_t *state) {
    uint8_t a_15 = state[0];
    uint8_t tmp[BLOCK_SIZE];

    for (int i = 0; i < BLOCK_SIZE - 1; ++i) {
        tmp[i] = state[i + 1];
        a_15 ^= GF_mul(tmp[i], rcon[i]);
    }

    tmp[15] = a_15;
    memcpy(state, tmp, BLOCK_SIZE);
}

void Kuznyechik::L(const uint8_t *input, uint8_t *output) {
    uint8_t tmp[BLOCK_SIZE];
    memcpy(tmp, input, BLOCK_SIZE);

    for (int i = 0; i < BLOCK_SIZE; ++i) {
        R(tmp);
    }

    memcpy(output, tmp, BLOCK_SIZE);
}

void Kuznyechik::reverse_L(const uint8_t *input, uint8_t *output) {
    uint8_t tmp[BLOCK_SIZE];
    memcpy(tmp, input, BLOCK_SIZE);

    for (int i = 0; i < 16; ++i) {
        reverse_R(tmp);
    }

    memcpy(output, tmp, BLOCK_SIZE);
}

void Kuznyechik::encrypt(const unsigned char *input, unsigned char *output) {
    memcpy(output, input, BLOCK_SIZE);

    for (int i = 0; i < 9; ++i) {
        X(output, round_keys[i]);
        S(output, output);
        L(output, output);
    }

    X(output, round_keys[9]);
}

void Kuznyechik::LSX(uint8_t *input, const uint8_t *key) {
    X(input, key);

    unsigned char tmp[BLOCK_SIZE]{};

    for (size_t i = 0; i < BLOCK_SIZE; ++i) {
        X(tmp, lut_table[i][input[i]]);
    }

    memcpy(input, tmp, BLOCK_SIZE);
}


void Kuznyechik::encrypt2(unsigned char *input, unsigned char *output) const {
    memcpy(output, input, BLOCK_SIZE);

    for (int i = 0; i < 9; ++i) {
        LSX(output, round_keys[i]);
    }

    X(output, round_keys[9]);
}

void Kuznyechik::LSX_intr(uint8_t *input, const uint8_t *key) {
    auto *x = reinterpret_cast<__m128i *>(input);
    const auto *k = reinterpret_cast<const __m128i *>(key);

    *x = _mm_xor_si128(*x, *k);

    uint8_t tmp[BLOCK_SIZE]{};
    uint8_t *tmp_p{tmp};

    for (size_t i = 0; i < BLOCK_SIZE; ++i) {
        auto *t = reinterpret_cast<__m128i *>(tmp_p);

        uint8_t *lut_p{lut_table[i][input[i]]};
        const auto *lut = reinterpret_cast<const __m128i *>(lut_p);

        *t = _mm_xor_si128(*t, *lut);
    }

    memcpy(input, tmp, BLOCK_SIZE);
}

void Kuznyechik::encrypt3(unsigned char *input, unsigned char *output) const {
    memcpy(output, input, BLOCK_SIZE);

    for (int i = 0; i < 9; ++i) {
        LSX_intr(output, round_keys[i]);
    }

    auto *x = reinterpret_cast<__m128i *>(output);
    const auto *k = reinterpret_cast<const __m128i *>(round_keys[9]);

    *x = _mm_xor_si128(*x, *k);
}

void Kuznyechik::decrypt(const unsigned char *input, unsigned char *output) {
    memcpy(output, input, BLOCK_SIZE);
    X(output, round_keys[9]);

    for (int i = 8; i >= 0; --i) {
        reverse_L(output, output);
        reverse_S(output, output);
        X(output, round_keys[i]);
    }
}

int main() {
    unsigned char key[Kuznyechik::KEY_SIZE] = {
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10,
        0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef
    };

    Kuznyechik kuz{key};

    /*unsigned char s_input[Kuznyechik::BLOCK_SIZE] = {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
                                                     0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x00};

    unsigned char s_output[Kuznyechik::BLOCK_SIZE];
    Kuznyechik::S(s_input, s_output);

    std::cout << "S:" << '\n';

    for (unsigned char i: s_output) {
        std::cout << std::hex << (int) i << ' ';
    }

    std::cout << '\n';

    unsigned char r_input[Kuznyechik::BLOCK_SIZE] = {0x94, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01};
    Kuznyechik::R(r_input);

    std::cout << "R:" << '\n';

    for (unsigned char i: r_input) {
        std::cout << std::hex << (int) i << ' ';
    }

    std::cout << '\n';

    unsigned char l_input[Kuznyechik::BLOCK_SIZE] = {0x64, 0xa5, 0x94, 0x00, 0x00, 0x00, 0x00, 0x00,
                                                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

    unsigned char l_output[Kuznyechik::BLOCK_SIZE];
    Kuznyechik::L(l_input, l_output);

    std::cout << "L:" << '\n';

    for (unsigned char i: l_output) {
        std::cout << std::hex << (int) i << ' ';
    }

    std::cout << '\n';

    std::cout << "ROUND_KEYS:" << '\n';

    for (auto & round_key : kuz.round_keys) {
        for (int j = 0; j < Kuznyechik::BLOCK_SIZE; ++j) {
            uint8_t el = round_key[j];
            std::cout << std::hex << (int) el << ' ';
        }
        std::cout << '\n';
    }*/

    /*unsigned char a[Kuznyechik::BLOCK_SIZE] = { 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x00,
                                                0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88 };

    unsigned char a_output[Kuznyechik::BLOCK_SIZE]{};

    unsigned char a_decrypt[Kuznyechik::BLOCK_SIZE]{};

    kuz.encrypt3(a, a_output);
    kuz.decrypt(a_output, a_decrypt);

    std::cout << "ENCRYPT:" << '\n';

    for (unsigned char el : a_output) {
        std::cout << std::hex << (int)el << ' ';
    }

    std::cout << '\n';

    std::cout << "DECRYPT:" << '\n';

    for (unsigned char el : a_decrypt) {
        std::cout << std::hex << (int)el << ' ';
    }

    std::cout << '\n';*/


    unsigned char b[Kuznyechik::BLOCK_SIZE] = {
        0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x00,
        0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88
    };

    //unsigned char c[Kuznyechik::BLOCK_SIZE]{};

    double num_mb = 100;
    double SIZE = static_cast<double>(1024 * 1024) / static_cast<double>(Kuznyechik::BLOCK_SIZE); // 1 MB
    double num_for = SIZE * num_mb;

    std::cout << "Num of for: " << static_cast<int>(num_for) << '\n';

    int num = static_cast<int>(num_for);

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num; ++i) {
        kuz.encrypt3(b, b);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Speed: " << std::dec << num_mb * 1000 / static_cast<double>(duration.count()) << " MB/s" <<
            std::endl;

    std::cout << "FOR_ENCRYPT:" << '\n';

    for (unsigned char el: b) {
        std::cout << std::hex << static_cast<int>(el) << ' ';
    }

    std::cout << '\n';

    return 0;
}
