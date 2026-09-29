/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Vinicius P. M. Miguel
 */
#ifndef MFCC_H
#define MFCC_H

#include <stdint.h>

#include "dct.h"
#include "fft_fp.h"
#include "mel.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MFCC_FRAME_SIZE_SEC 0.025
#define MFCC_FRAME_STEP_SEC 0.010

#define MFCC_NFFT         NFFT
#define MFCC_NUM_FILTERS  NUM_FILTERS
#define MFCC_NUM_CEPS     NUM_CEPS

typedef struct {
  int F_PRE;
  int F_HAMMING;
  int F_FFT;
  int F_MEL;
  int F_DCT;

  int TRUNCATE_PRE;
  int TRUNCATE_HAMMING;
  int TRUNCATE_FFT;
  int TRUNCATE_MEL;
  int TRUNCATE_DCT;
} mfcc_config_t;

typedef struct {
  int64_t *pre_emphasis;
  int num_samples;

  int32_t *window;
  int frame_size;

  complex_t *twiddles;
  int32_t **filterbank;
  int16_t max_width_mel;

  int64_t **frames;
  int64_t **hamming_frames;
  int64_t **power_spectrum;
  int32_t **energies;
  int32_t **coefficients;

  int num_frames;
  int num_ceps;
} mfcc_result_t;

int mfcc_compute(const int16_t *samples, int num_samples, int sample_rate,
                 const mfcc_config_t *config, mfcc_result_t *result);

void mfcc_free_result(mfcc_result_t *result);

#ifdef __cplusplus
}
#endif

#endif