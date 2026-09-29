# MFCC Fixed-Point (C)

C implementation of the Mel-Frequency Cepstral Coefficients (MFCC) algorithm using configurable fixed-point arithmetic. The fractional bit width (F) can be independently adjusted at each processing stage, enabling flexible precision and dynamic control of numerical representation.

## Layout
```bash
.
├── include/mfcc/mfcc.h # public API
├── src/ # implementation
│ ├── mfcc.c
│ ├── fft.c
│ ├── fft_fp.c
│ ├── mel.c
│ ├── dct.c
│ └── process.c
├── build/ # compiled objects + libmfcc.a
├── Makefile
└── shell.nix
```

## Build

```bash
make
```
Produces `build/libmfcc.a`.

## API

int mfcc_compute(const int16_t *samples, int num_samples, int sample_rate,
const mfcc_config_t *config, mfcc_result_t *result);
void mfcc_free_result(mfcc_result_t *result);

See include/mfcc/mfcc.h for `mfcc_config_t` and `mfcc_result_t`.
Configuration

Each stage has its own fractional width and optional truncation flag:
| Field | Stage |
|---|---|
| `F_PRE` / `TRUNCATE_PRE` | Pre-emphasis |
| `F_HAMMING` / `TRUNCATE_HAMMING` | Hamming window |
| `F_FFT` / `TRUNCATE_FFT` | FFT / power spectrum |
| `F_MEL` / `TRUNCATE_MEL` | Mel filterbank |
| `F_DCT` / `TRUNCATE_DCT` | DCT |