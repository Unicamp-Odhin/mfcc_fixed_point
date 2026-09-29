#include <stdint.h>
#include <stdlib.h>

#include "mfcc/mfcc.h"

#include "dct.h"
#include "fft.h"
#include "fft_fp.h"
#include "mel.h"
#include "process.h"

int mfcc_compute(const int16_t *samples, int num_samples, int sample_rate,
                 const mfcc_config_t *config, mfcc_result_t *output) {
  if (!samples || !config || !output)
    return -1;

  int frame_size = (int)(sample_rate * MFCC_FRAME_SIZE_SEC);

  int frame_step = (int)(sample_rate * MFCC_FRAME_STEP_SEC);

  int num_frames = (int)((double)(num_samples - frame_size) / frame_step) + 1;

  if (num_frames <= 0)
    return -1;

  // PRE-EMPHASIS
  int64_t *samples_64bit = malloc(sizeof(int64_t) * num_samples);

  if (!samples_64bit)
    return -1;

  pre_emphasis(samples, num_samples, samples_64bit, config->F_PRE);

  if (config->TRUNCATE_PRE) {
    int64_t mask = ~((1LL << config->F_PRE) - 1);

    for (int i = 0; i < num_samples; i++)
      samples_64bit[i] &= mask;
  }

  output->pre_emphasis = samples_64bit;
  output->num_samples = num_samples;

  // TABLES
  int32_t *window = malloc(frame_size * sizeof(int32_t));

  if (!window) {
    free(samples_64bit);
    return -1;
  }

  generate_hamming_window(window, frame_size, config->F_HAMMING);

  output->window = window;
  output->frame_size = frame_size;

  complex_t *twiddles = malloc((MFCC_NFFT / 2) * sizeof(complex_t));

  if (!twiddles) {
    free(samples_64bit);
    free(window);
    return -1;
  }

  generate_twiddles(twiddles, MFCC_NFFT, config->F_FFT);

  output->twiddles = twiddles;

  int32_t **filterbank = malloc(MFCC_NUM_FILTERS * sizeof(int32_t *));

  if (!filterbank) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    return -1;
  }

  int16_t max_width_mel =
      create_op_filterbank(filterbank, sample_rate, config->F_MEL);

  output->filterbank = filterbank;
  output->max_width_mel = max_width_mel;

  init_cos_lut(config->F_DCT);

  // FRAMES
  int64_t **frames = frame_signal_int(samples_64bit, num_samples, frame_size,
                                      frame_step, &num_frames);

  output->frames = frames;

  int64_t **hamming_frames = malloc(num_frames * sizeof(int64_t *));

  if (!hamming_frames) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    free(filterbank);
    free(frames);
    return -1;
  }

  int64_t **power_spectrum = malloc(num_frames * sizeof(int64_t *));

  if (!power_spectrum) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    free(filterbank);
    free(frames);
    free(hamming_frames);
    return -1;
  }

  int32_t **energies = malloc(num_frames * sizeof(int32_t *));

  if (!energies) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    free(filterbank);
    free(frames);
    free(hamming_frames);
    free(power_spectrum);
    return -1;
  }

  // OUTPUT
  int32_t **ceps = malloc(num_frames * sizeof(int32_t *));

  if (!ceps) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    free(filterbank);
    free(frames);
    free(hamming_frames);
    free(power_spectrum);
    free(energies);
    return -1;
  }

  for (int i = 0; i < num_frames; i++) {

    ceps[i] = malloc(MFCC_NUM_CEPS * sizeof(int32_t));

    if (!ceps[i]) {
      return -1;
    }

    hamming_frames[i] = malloc(frame_size * sizeof(int64_t));

    if (!hamming_frames[i]) {
      return -1;
    }

    power_spectrum[i] = malloc(MFCC_NFFT * sizeof(int64_t));

    if (!power_spectrum[i]) {
      return -1;
    }

    energies[i] = malloc(MFCC_NUM_FILTERS * sizeof(int32_t));

    if (!energies[i]) {
      return -1;
    }

    for (int j = 0; j < frame_size; j++)
      hamming_frames[i][j] = frames[i][j];

    // HAMMING
    hamming_window_fixed(hamming_frames[i], window, frame_size,
                         config->F_HAMMING, config->F_PRE);

    if (config->TRUNCATE_HAMMING) {

      int64_t mask = ~((1LL << config->F_HAMMING) - 1);

      for (int j = 0; j < frame_size; j++)
        hamming_frames[i][j] &= mask;
    }

    // FFT
    power_spectrum[i][0] = 0;

    fft_real_power(hamming_frames[i], frame_size, power_spectrum[i], twiddles,
                   config->F_FFT, config->F_HAMMING);

    if (config->TRUNCATE_FFT) {

      int64_t mask = ~((1LL << config->F_FFT) - 1);

      for (int j = 0; j < MFCC_NFFT; j++)
        power_spectrum[i][j] &= mask;
    }

    // MEL
    apply_op_filterbank(power_spectrum[i], energies[i], sample_rate, filterbank,
                        config->F_MEL, 16, config->F_FFT);

    if (config->TRUNCATE_MEL) {

      int64_t mask = ~((1LL << 16) - 1);

      for (int j = 0; j < MFCC_NUM_FILTERS; j++)
        energies[i][j] &= mask;
    }

    // DCT
    dct_fixed(energies[i], MFCC_NUM_FILTERS, ceps[i], 16, config->F_DCT);

    if (config->TRUNCATE_DCT) {

      int64_t mask = ~((1LL << config->F_DCT) - 1);

      for (int j = 0; j < MFCC_NUM_CEPS; j++)
        ceps[i][j] &= mask;
    }
  }

  output->hamming_frames = hamming_frames;
  output->power_spectrum = power_spectrum;
  output->energies = energies;
  output->coefficients = ceps;
  output->num_frames = num_frames;
  output->num_ceps = MFCC_NUM_CEPS;

  return 0;
}

void mfcc_free_result(mfcc_result_t *output) {
  if (!output)
    return;

  if (output->coefficients) {
    for (int i = 0; i < output->num_frames; i++)
      free(output->coefficients[i]);
    free(output->coefficients);
  }

  if (output->energies) {
    for (int i = 0; i < output->num_frames; i++)
      free(output->energies[i]);
    free(output->energies);
  }

  if (output->power_spectrum) {
    for (int i = 0; i < output->num_frames; i++)
      free(output->power_spectrum[i]);
    free(output->power_spectrum);
  }

  if (output->hamming_frames) {
    for (int i = 0; i < output->num_frames; i++)
      free(output->hamming_frames[i]);
    free(output->hamming_frames);
  }

  if (output->frames) {
    for (int i = 0; i < output->num_frames; i++)
      free(output->frames[i]);
    free(output->frames);
  }

  if (output->filterbank) {
    for (int i = 0; i < MFCC_NUM_FILTERS; i++)
      free(output->filterbank[i]);
    free(output->filterbank);
  }

  free(output->twiddles);
  free(output->window);
  free(output->pre_emphasis);

  output->coefficients = NULL;
  output->energies = NULL;
  output->power_spectrum = NULL;
  output->hamming_frames = NULL;
  output->frames = NULL;
  output->filterbank = NULL;
  output->twiddles = NULL;
  output->window = NULL;
  output->pre_emphasis = NULL;
  output->num_frames = 0;
  output->num_ceps = 0;
  output->num_samples = 0;
  output->frame_size = 0;
  output->max_width_mel = 0;
}