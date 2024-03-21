// PUBLIC IMPORT HEADER
#ifndef _UINT256_H_
#define _UINT256_H_
#include "uint256_t_config.include"
#define UINT256_T_EXTERN _UINT256_T_IMPORT
typedef __uint128_t uint128_t;
const uint128_t uint128_0(0);
const uint128_t uint128_1(1);
#ifndef __LITTLE_ENDIAN__
#ifndef __BIG_ENDIAN__
#define __LITTLE_ENDIAN__ 1
#endif
#endif
#include "uint256_t.include"
#endif
