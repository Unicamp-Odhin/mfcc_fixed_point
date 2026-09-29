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

  //PRE-EMPHASIS
  int64_t *samples_64bit = malloc(sizeof(int64_t) * num_samples);

  if (!samples_64bit)
    return -1;

  pre_emphasis(samples, num_samples, samples_64bit, config->F_PRE);

  if (config->TRUNCATE_PRE) {
    int64_t mask = ~((1LL << config->F_PRE) - 1);

    for (int i = 0; i < num_samples; i++)
      samples_64bit[i] &= mask;
  }

  //TABLES
  int32_t *window = malloc(frame_size * sizeof(int32_t));

  if (!window) {
    free(samples_64bit);
    return -1;
  }

  generate_hamming_window(window, frame_size, config->F_HAMMING);

  complex_t *twiddles = malloc((MFCC_NFFT / 2) * sizeof(complex_t));

  if (!twiddles) {
    free(samples_64bit);
    free(window);
    return -1;
  }

  generate_twiddles(twiddles, MFCC_NFFT, config->F_FFT);

  int32_t **filterbank = malloc(MFCC_NUM_FILTERS * sizeof(int32_t *));

  if (!filterbank) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    return -1;
  }

  create_op_filterbank(filterbank, sample_rate, config->F_MEL);

  init_cos_lut(config->F_DCT);

  //FRAMES
  int64_t **frames = frame_signal_int(samples_64bit, num_samples, frame_size,
                                      frame_step, &num_frames);

  //OUTPUT
  int32_t **ceps = malloc(num_frames * sizeof(int32_t *));

  if (!ceps) {
    free(samples_64bit);
    free(window);
    free(twiddles);
    free(filterbank);
    free(frames);
    return -1;
  }

  for (int i = 0; i < num_frames; i++) {

    ceps[i] = malloc(MFCC_NUM_CEPS * sizeof(int32_t));

    if (!ceps[i]) {
      return -1;
    }

    //HAMMING
    hamming_window_fixed(frames[i], window, frame_size, config->F_HAMMING,
                         config->F_PRE);

    if (config->TRUNCATE_HAMMING) {

      int64_t mask = ~((1LL << config->F_HAMMING) - 1);

      for (int j = 0; j < frame_size; j++)
        frames[i][j] &= mask;
    }

    //FFT
    int64_t *power_spectrum = malloc(MFCC_NFFT * sizeof(int64_t));
    int32_t *energies = malloc(MFCC_NUM_FILTERS * sizeof(int32_t));

    if (!power_spectrum || !energies)
      return -1;

    power_spectrum[0] = 0;

    fft_real_power(frames[i], frame_size, power_spectrum, twiddles,
                   config->F_FFT, config->F_HAMMING);

    if (config->TRUNCATE_FFT) {

      int64_t mask = ~((1LL << config->F_FFT) - 1);

      for (int j = 0; j < MFCC_NFFT; j++)
        power_spectrum[j] &= mask;
    }

    //MEL
    apply_op_filterbank(power_spectrum, energies, sample_rate, filterbank,
                        config->F_MEL, 16, config->F_FFT);

    if (config->TRUNCATE_MEL) {

      int64_t mask = ~((1LL << 16) - 1);

      for (int j = 0; j < MFCC_NUM_FILTERS; j++)
        energies[j] &= mask;
    }

    //DCT
    dct_fixed(energies, MFCC_NUM_FILTERS, ceps[i], 16, config->F_DCT);

    if (config->TRUNCATE_DCT) {

      int64_t mask = ~((1LL << config->F_DCT) - 1);

      for (int j = 0; j < MFCC_NUM_CEPS; j++)
        ceps[i][j] &= mask;
    }

    free(power_spectrum);
    free(energies);
  }

  output->coefficients = ceps;
  output->num_frames = num_frames;
  output->num_ceps = MFCC_NUM_CEPS;

  free(samples_64bit);
  free(window);
  free(twiddles);
  free(filterbank);
  free(frames);

  return 0;
}

void mfcc_free_result(mfcc_result_t *output) {
  if (!output || !output->coefficients)
    return;

  for (int i = 0; i < output->num_frames; i++)
    free(output->coefficients[i]);

  free(output->coefficients);

  output->coefficients = NULL;
  output->num_frames = 0;
  output->num_ceps = 0;
}
