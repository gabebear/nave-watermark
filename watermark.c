#include "watermark.h"

#define WINDOW_SIZE              256
#define WINDOW_INCREMENT         4
#define NEIGHBORHOOD_START_INDEX 45
#define NEIGHBORHOOD_BIN_SIZE    25

// frequency indices we care about
#define I1    48
#define IMID  55
#define I0    62

#define REF_I 53 // reference frequency index 5 Khz

#define SHIFT 5

typedef double complex fft_t;

void print_fft(fft_t *r);

int pick_bit(fft_t *, int);
void calculate_fft(short *, size_t, fftw_complex *);
void sliding_fft(fftw_complex *, short, short);
int test_triple_tone(fft_t *);
void convolute_window_function(fftw_complex *, fft_t *);
double bin_power(fft_t);
void max_min_bin(fft_t *, int, int *, double *, int *, double *);

// status information
typedef struct sis {
  int pcs;
  int ji;
  int gc;
  int da[15];
  int op[10]; // 5 bits each
} sis_t;

int jump[] = { 2, 5, 1, 4, 3, 2, 5 };

// status information array
sis_t SIS[64];

int main(int argc, const char *argv[])
{
  int sample = 0;
  int i, ret, tt, p = 0;
  SNDFILE *snd_file;
  SF_INFO sf_info;
  fifo_t *fifo;
  short samples[WINDOW_SIZE], oldest_samples[WINDOW_INCREMENT];
  fftw_complex *fft  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * WINDOW_SIZE);

  memset(&SIS, 0, sizeof(sis_t) * 64);

  // initialize fifo queue
  fifo = fifo_init(WINDOW_SIZE, sizeof(short));

  snd_file = sf_open("./samples/homeland2.aiff", SFM_READ, &sf_info);

  ret = sf_read_short(snd_file, (short *)&samples, WINDOW_SIZE);

  if(ret != WINDOW_SIZE) {
    puts("not enough samples in file");
    return -1;
  }

  // add all samples to the fifo
  for(i = 0; i < WINDOW_SIZE; i++)
    fifo_push(fifo, &samples[i]);

  // perform first FFT calculation
  calculate_fft((short *)&samples, WINDOW_SIZE, fft);

  sample += 256;

  // start incremental mode..
  do {
    ret = sf_read_short(snd_file, (short *)&samples, WINDOW_INCREMENT);

    if(ret == 0)
      break;

    // remove oldest samples, add newest samples
    for(i = 0; i < WINDOW_INCREMENT; i++) {
      fifo_pop(fifo, &oldest_samples[i]);
      fifo_push(fifo, &samples[i]);
      sliding_fft(fft, oldest_samples[i], samples[i]);
    }
    sample += 4;

    // TODO: clean up, get rid of fft_complex data structure, just use complex double typedef
    fft_t r[WINDOW_SIZE];
    convolute_window_function(fft, (fft_t *)&r);

    if(SIS[p].pcs > 0) {
      int ji = jump[SIS[p].ji % 7];
      tt = pick_bit((fft_t *)&r, ji);

      if(tt >= 0) {
        SIS[p].da[SIS[p].ji] = tt;

        printf("p = %d, value = ", p);
        for(i = 0; i < SIS[p].ji; i++) {
          printf("%d", SIS[p].da[i]);
        }
        puts("");
        SIS[p].ji++;
      } else {
        SIS[p].pcs = 0;
        SIS[p].ji = 0;
      }
    } else {
      tt = test_triple_tone((fft_t *)&r);
      if(tt >= 0) {
        SIS[p].pcs = 1;
        SIS[p].ji = 1;
        SIS[p].da[0] = tt;
      }
    }

    if(++p == 64)
      p = 0;

  } while(ret > 0);

  return 0;
}

void print_fft(fft_t *r)
{
  int i;
  printf("[\n");
  for(i = 0; i < WINDOW_SIZE; i++) {
    printf("%f + %fi,\n", creal(r[i]), cimag(r[i]));
  }
  printf("]\n\n");
}

int pick_bit(fft_t *r, int ji)
{
  int i0, i1;
  int index, imax, imin;
  double fmax, fmin;

  i0 = i1 = -1;
  index = REF_I + ji - SHIFT; // i1
  max_min_bin(r, index, &imax, &fmax, &imin, &fmin);

  if(index == imax || abs(index - imax) == 1) {
    i1 = 1;
  } else if(index == imin || abs(index - imin) == 1) {
    i1 = 0;
  }

  index = REF_I + ji + SHIFT; // i0
  max_min_bin(r, index, &imax, &fmax, &imin, &fmin);

  if(index == imax || abs(index - imax) == 1) {
    i0 = 1;
  } else if(index == imin || abs(index - imin) == 1) {
    i0 = 0;
  }

  if(i1 == 1 && i0 == 0)
    return 1;

  if(i0 == 1 && i1 == 0)
    return 0;

  return -1;
}

int test_triple_tone(fft_t *fft)
{
  int i0, imid, i1, i0min, imidmin, i1min;
  double f0, fmid, f1, f0min, f1min, fmidmin;

  max_min_bin(fft, I1, &i1, &f1, &i1min, &f1min);
  max_min_bin(fft, IMID, &imid, &fmid, &imidmin, &fmidmin);
  max_min_bin(fft, I0, &i0, &f0, &i0min, &f0min);

  if(i1 == I1 && imid == IMID && i0 == I0) {
    return f1 > f0 ? 1 : 0;
  }

/*
  if(i1min == I1 && imidmin == IMID && i0min == I0) {
    return f1min > f0min ? 1 : 0;
  }
*/
  return -1;
}

void max_min_bin(fft_t *fft, int base_index, int *max_index, double *max_freq, int *min_index, double *min_freq)
{
  int i, mini, maxi;
  double value, minf, maxf;

  maxf = 0;
  minf = INT16_MAX;

  for(i = base_index - 2; i <= base_index + 2; i++) {
    value = bin_power(fft[i]);

    if(value >= maxf) {
      maxf = value;
      maxi = i;
    }

    if(value <= minf) {
      minf = value;
      mini = i;
    }
  }

  *max_freq = maxf;
  *min_freq = minf;
  *max_index = maxi;
  *min_index = mini;
}

void sliding_fft(fftw_complex *fft, short old, short new)
{
  int i;
  fft_t p, e, result;

  // goertzel sliding DFT
  for(i = NEIGHBORHOOD_START_INDEX; i <= NEIGHBORHOOD_START_INDEX + NEIGHBORHOOD_BIN_SIZE; i++) {
    p = fft[i][0] + fft[i][1] * I;

    e = cos(2 * M_PI * i / WINDOW_SIZE) + sin(2 * M_PI * i / WINDOW_SIZE) * I;
    result = (p - old + new) * e;

    fft[i][0] = creal(result);
    fft[i][1] = cimag(result);
  }
}

void
convolute_window_function(fftw_complex *fft, fft_t *r)
{
  int i;

  for(i = 0; i < WINDOW_SIZE; i++) {
    r[i] = fft[i][0] + fft[i][1] * I;
  }

  // apply hanning window in frequency domain around bins that contain our signal
  for(i = NEIGHBORHOOD_START_INDEX; i <= NEIGHBORHOOD_START_INDEX + NEIGHBORHOOD_BIN_SIZE; i++) {
    fft_t a, b, e;

    a = fft[i-1][0] + fft[i-1][1] * I;
    b = fft[i+1][0] + fft[i+1][1] * I;

    e = cos(-2 * M_PI * i / WINDOW_SIZE) + sin(-2 * M_PI * i / WINDOW_SIZE) * I;
    a = -0.5 * a;// * e;

    e = cos(2 * M_PI * i / WINDOW_SIZE) + sin(2 * M_PI * i / WINDOW_SIZE) * I;
    b = -0.5 * b;// * e;

    r[i] = 0.5 * (r[i] + a + b);
  }
}

void calculate_fft(short *samples, size_t total, fftw_complex *result)
{
  fftw_complex    *data;
  fftw_plan       plan_forward;
  int             i;

  data        = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * total );
  plan_forward  = fftw_plan_dft_1d( total, data, result, FFTW_FORWARD, FFTW_ESTIMATE );

  /* populate input data */
  for( i = 0 ; i < total ; i++ ) {
    data[i][0] = samples[i];
    data[i][1] = 0.0;
  }

  fftw_execute( plan_forward );
}

double
bin_power(fft_t bin)
{
  return sqrt((creal(bin) * creal(bin)) + (cimag(bin) * cimag(bin)));
}
