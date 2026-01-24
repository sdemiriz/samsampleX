/*
 * xxHash - Extremely Fast Hash algorithm
 * Header File
 * Copyright (C) 2012-2021 Yann Collet
 *
 * BSD 2-Clause License (https://www.opensource.org/licenses/bsd-license.php)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following disclaimer
 *      in the documentation and/or other materials provided with the
 *      distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This is a minimal single-header implementation of xxHash32 for samsampleX.
 * For the full implementation, see: https://github.com/Cyan4973/xxHash
 */

#ifndef XXHASH_H
#define XXHASH_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* xxHash32 constants */
#define XXH_PRIME32_1  0x9E3779B1U
#define XXH_PRIME32_2  0x85EBCA77U
#define XXH_PRIME32_3  0xC2B2AE3DU
#define XXH_PRIME32_4  0x27D4EB2FU
#define XXH_PRIME32_5  0x165667B1U

/* Rotate left */
static inline uint32_t XXH_rotl32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

/* Read 32-bit little-endian */
static inline uint32_t XXH_read32(const void *ptr) {
    uint32_t val;
    memcpy(&val, ptr, sizeof(val));
    return val;
}

/* Round function */
static inline uint32_t XXH32_round(uint32_t acc, uint32_t input) {
    acc += input * XXH_PRIME32_2;
    acc = XXH_rotl32(acc, 13);
    acc *= XXH_PRIME32_1;
    return acc;
}

/* Avalanche function */
static inline uint32_t XXH32_avalanche(uint32_t h32) {
    h32 ^= h32 >> 15;
    h32 *= XXH_PRIME32_2;
    h32 ^= h32 >> 13;
    h32 *= XXH_PRIME32_3;
    h32 ^= h32 >> 16;
    return h32;
}

/*
 * XXH32 - Calculate 32-bit hash
 *
 * @input: pointer to input data
 * @len: length of input in bytes
 * @seed: seed value for hash
 * @return: 32-bit hash value
 */
static inline uint32_t XXH32(const void *input, size_t len, uint32_t seed) {
    const uint8_t *p = (const uint8_t *)input;
    const uint8_t *const bEnd = p + len;
    uint32_t h32;

    if (len >= 16) {
        const uint8_t *const limit = bEnd - 15;
        uint32_t v1 = seed + XXH_PRIME32_1 + XXH_PRIME32_2;
        uint32_t v2 = seed + XXH_PRIME32_2;
        uint32_t v3 = seed + 0;
        uint32_t v4 = seed - XXH_PRIME32_1;

        do {
            v1 = XXH32_round(v1, XXH_read32(p)); p += 4;
            v2 = XXH32_round(v2, XXH_read32(p)); p += 4;
            v3 = XXH32_round(v3, XXH_read32(p)); p += 4;
            v4 = XXH32_round(v4, XXH_read32(p)); p += 4;
        } while (p < limit);

        h32 = XXH_rotl32(v1, 1) + XXH_rotl32(v2, 7) +
              XXH_rotl32(v3, 12) + XXH_rotl32(v4, 18);
    } else {
        h32 = seed + XXH_PRIME32_5;
    }

    h32 += (uint32_t)len;

    /* Process remaining bytes */
    while (p + 4 <= bEnd) {
        h32 += XXH_read32(p) * XXH_PRIME32_3;
        h32 = XXH_rotl32(h32, 17) * XXH_PRIME32_4;
        p += 4;
    }

    while (p < bEnd) {
        h32 += (*p++) * XXH_PRIME32_5;
        h32 = XXH_rotl32(h32, 11) * XXH_PRIME32_1;
    }

    return XXH32_avalanche(h32);
}

/*
 * Hash a string and return a fraction in [0, 1)
 * This is useful for probabilistic sampling.
 */
static inline double XXH32_fraction(const char *str, uint32_t seed) {
    uint32_t hash = XXH32(str, strlen(str), seed);
    return (double)hash / (double)UINT32_MAX;
}

#endif /* XXHASH_H */

