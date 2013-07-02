#include "clink.h"
#include <math.h>
#include <assert.h>

socklen_t salen = sizeof (Sockaddr);

char *host;                  /* name of the destination machine */
char *dump_file = NULL;      /* name of the dump file or NULL */
FILE *dump_file_fp = NULL;   /* file pointer for dump file */
int dump_minima = 0;         /* flag: should we write mins files? */
char *input_file = NULL;     /* name of input file or NULL */
FILE *input_file_fp = NULL;  /* file pointer for dump file */
int min_ttl = 1;             /* starting ttl */
int max_ttl = 30;            /* ending ttl, unless we reach dest first */
int kernel_timestamps = 0;   /* flag: should we use kernel-level timestamps? */
int num_probes = 8;          /* how many probes at each size */
int low = 28;                /* smallest packet size (in Bytes) */
int high = 1500;             /* largest packet size */
int step = 16;               /* interval between packet sizes */
int tos = 0;                 /* type of service in outgoing probes */
int wait_time = 2;           /* how many seconds to wait before timeout */
int dns_resolve = 1;         /* flag: should we use DNS to resolve addrs? */
int verbose = 0;             /* flag: should we print verbose output */
/* how long should we wait between probes */
double inter_sample_ratio = 10.0;   /* a multiple of the previous rtt? */
double inter_sample_time = 0.0;     /* or a fixed interval in microsecs */

int probe_count;     /* how many probes have we sent to this link so far? */

int num_fiddles;     /* how many times we fiddle with the data */

/* paths: the big, hairy data structure that contains all the data
          we've collected.  paths is an array of Path structures;
	  each one contains alternate paths, but at the end we only
	  characterize the sequence of links contained in paths[0] */

#define MAX_ALT_PATH 10
Path *paths[MAX_ALT_PATH];

/* make_sample */

Sample *make_sample (int ttl, int packet_size, int space)
{
  Sample *new = (Sample *) malloc (sizeof (Sample));

  new->ttl = ttl;
  new->size = packet_size;
  if (space != 0) {
    new->x = (double *) malloc (space * sizeof (double));
  }
  new->n = 0;
  new->space = space;
  
  return new;
}

/* make_link */

Link *make_link (int ttl)
{
  int i;
  Link *new = (Link *) malloc (sizeof (Link));

  new->ttl = ttl;
  for (i=0; i<MAX_SIZE; i++) {
    new->sample[i] = NULL;
  }
  return new;
}

/* make_fit */

Fit *make_fit ()
{
  Fit *new = (Fit *) malloc (sizeof (Fit));

  new->b0 = 0.0;
  new->b1 = 0.0;
  return new;
}

/* make_path */

Path *make_path ()
{
  int i;
  Path *new = (Path *) malloc (sizeof (Path));

  for (i=0; i<MAX_LINKS; i++) {
    new->both[i] = NULL;
    new->evens[i] = NULL;
    new->odds[i] = NULL;
    bzero (&new->addr[i], salen);
  }
  return new;
}

/* make_sizes */

Sizes *make_sizes (int low, int high, int step)
{
  int i;
  Sizes *sizes = (Sizes *) malloc (sizeof (Sizes));
  int n = (high-low) / step + 1;
  int size = low;

  sizes->size = (int *) Malloc (n * sizeof (int));
  sizes->permute = (int *) Malloc (n * sizeof (int));
  sizes->n = n;

  for (i=0; i<n; i++) {
    sizes->size[i] = size;
    size += step;
  }
  return sizes;
}

/* random_int: return a random integer between low and
   high, including low and not including high */

int random_int (int low, int high)
{
  int x = random () % (high-low);
  return low + x;
}

/* swap_elts: swap two elements of an integer array (like sizes) */

void swap_elts (int *a, int i, int j)
{
  int x = a[i];
  a[i] = a[j];
  a[j] = x;
}

/* permute_arrays: put the numbers 1 to n-1 in an array and shuffle
   them.  This array is used to generate the elements of the size
   array in random order */

void permute_sizes (Sizes *s)
{
  int i, j;

  for (i=0; i<s->n; i++) {
    s->permute[i] = i;
  }

  for (i=0; i<s->n; i++) {
    j = random_int (i, s->n);
    swap_elts (s->permute, i, j);
  }
}

/* calc_moments: find the first three moments of the sample */

void calc_moments (Sample *s)
{
  int i;
  double d;
  double sx = 0.0;
  double sx2 = 0.0;
  double sx3 = 0.0;

  for (i=0; i < s->n; i++) {
    sx += s->x[i];
    sx2 += s->x[i] * s->x[i];
  }
  s->m1p = sx/s->n;
  s->m2p = sx2/s->n;

  sx2 = 0.0;
  s->m1 = s->m1p;
  for (i=0; i < s->n; i++) {
    d = s->x[i] - s->m1;
    sx2 += d*d;
    sx3 += d*d*d;
  }
  s->m2 = sx2/s->n;
  s->m3 = sx3/s->n;

  /* printf ("moms: m1 = %lf\tm2 = %lf\n", s->m1, s->m2); */
}

/* calc_log_moments: apply a log transform to the sample and find
   the moments of the transformed data */

void calc_log_moments (Sample *s)
{
  int i;
  int n = s->n;
  double *x = s->x;
  double d;
  double sx = 0.0;
  double sx2 = 0.0;
  double m1, m2;

  for (i=0; i < n; i++) {
    assert (x[i] >= 0.0);
    sx += log (x[i]);
  }
  s->zeta = sx/n;
 
  for (i=0; i < n; i++) {
    d = log(x[i]) - s->zeta;
    sx2 += d*d;
  }
  s->sig2 = sx2/n;

  m1 = exp (s->zeta + s->sig2/2);
  m2 = exp (2*s->zeta + 2*s->sig2) - m1*m1;
  /* printf ("logn: m1 = %lf\tm2 = %lf\n", m1, m2); */
}

/* compar: takes two pointers to integers and tells which is
   bigger.  Used by sort_array. */

int compar (const void *px, const void *py)
{
  double x = *(double *) px;
  double y = *(double *) py;
  if (x < y) return -1;
  if (x > y) return 1;
  return 0;
}

/* sort_array: this is a veneer on the library implementation
   of quicksort, used to generate order statistics */

void sort_array (double x[], int n)
{
  void *base = (void *) x;
  size_t nmemb = n;
  size_t size = sizeof (double);
  
  qsort (base, nmemb, size, compar);
}

/* calc_order_stats: finds the minimum and whatever other order
   stats are of interest */

void calc_order_stats (Sample *s)
{
  int n = s->n;
  double *x = s->x;
  int i, j, k;
  double pi, pj, pk;

  sort_array (x, n);
  s->min = x[0];

  /* calculate the linearly extrapolated 1st percentile */
  i = 0;
  pi = (double) (i+1) / (n+1);

  k = n/4-1;
  pk = (double) (k+1) / (n+1);

  s->lep = x[i] - pi * (x[k] - x[i]) / (pk - pi);
}


/* fit: find a linear least-squares fit for the arrays and put
   the result in the given Fit structure

   this started out as the routine from Numerical Recipes in C,
   Press et al., but I have modified it a bit, and added the R^2 stuff */

void fit (double x[], double y[], int n, Fit *f)
{
  int i;
  double t, meanx, meany;
  double sx = 0.0;
  double sx2 = 0.0;
  double sy = 0.0;
  double sy2 = 0.0;
 
  /* find the means of x and y */
  for (i=0; i < n; i++) {
    sx += x[i];
    sy += y[i];
  }
  meanx = sx / n;
  meany = sy / n;

  /* find the variance of x and y, and the parameters b0 and b1 */
  f->b1 = 0.0;
  for (i=0; i < n; i++) {
    t = y[i] - meany;
    sy2 += t*t;
    t = x[i] - meanx;
    sx2 += t*t;
    f->b1 += t * y[i];
  }
  f->b1 /= sx2;
  f->b0 = (sy - sx*f->b1) / n;
 
  /* find the significances */
  f->sig0 = sqrt ((1.0 + sx*sx / (n * sx2)) / n);
  f->sig1 = sqrt (1.0 / sx2);
 
  /* find the sum of the squared residuals */
  f->chi2 = 0.0;
  for (i=0; i < n; i++) {
    t = y[i] - (f->b0 + f->b1 * x[i]);
    f->chi2 += t*t;
  }
  t = sqrt (f->chi2 / (n-2));
  f->sig0 *= t;
  f->sig1 *= t;

  /* r2 measures the size of the residuals relative to the original
     variance in y */
  f->r2 = 1 - f->chi2/sy2;
}

/* wfit: Same as fit except that the least squares fit is weighted */

void wfit (double x[], double y[], double w[], int n, Fit *f)
{
  int i;
  double ss = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double mx, t;
  double st2 = 0.0;
 
  for (i=0; i < n; i++) {
    ss += w[i];
    sx += x[i] * w[i];
    sy += y[i] * w[i];
  }
  mx = sx / ss;

  f->b1 = 0.0;
  for (i=0; i < n; i++) {
    t = x[i] - mx;
    st2 += t*t * w[i];
    f->b1 += t * y[i] * w[i];
  }
 
  f->b1 /= st2;
  f->b0 = (sy - sx*f->b1) / ss;

  f->sig0 = sqrt ((1.0 + sx*sx / (n * st2)) / ss);
  f->sig1 = sqrt (1.0 / st2);

  f->chi2 = 0.0; 
  for (i=0; i < n; i++) {
    t = y[i] - f->b0 - f->b1 * x[i];
    f->chi2 += t*t * w[i];
  }
}

/* find_weights: we assume that data far from the fitted line
   are less reliable, and give them weights accordingly */

void find_weights (double x[], double y[], double w[], int n, Fit *f)
{
  int i;
  double fit, res, err;
  double max = 0.0;

  /* find the datum with the largest abs(residual) */
  for (i=0; i<n; i++) {
    fit = x[i] * f->b1 + f->b0;
    res = y[i] - fit;
    err = fabs (res);
    w[i] = err;
    if (w[i] > max) max = w[i];
  }

  /* scale the weights so that the most distant outlier is
     ignored altogether, and other weights are proportional to error */
  for (i=0; i<n; i++) {
    w[i] = 1.0 - w[i]/max;
  }
}

/* print_fit: used for debugging.  prints the actual values (xs and ys),
   and then the actual xs with the fitted value for each s  */

void print_fit (double x[], double y[], int n, Fit *f)
{
  int i;
  double fit;

  printf ("\n");
  for (i=0; i<n; i++) {
    printf ("%lf\t%lf\n", x[i], y[i]);
  }

  printf ("\n");
  for (i=0; i<n; i++) {
    fit = x[i] * f->b1 + f->b0;
    printf ("%lf\t%lf\n", x[i], fit);
  }
}

/* fit_link: least squares fit for a link */

void fit_link (Link *link)
{
  int i;
  double x[MAX_SIZE];
  double y[MAX_SIZE];
  double w[MAX_SIZE];
  int n = 0;

  /* collect the message size and rtt for each sample */
  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      x[n] = (double) i;
      y[n] = link->sample[i]->min;
      n++;
    }
  }

  fit (x, y, n, link->fit);

  /* here's the code for iteratively-weighted least squares (IWLS)
  find_weights (x, y, w, n, link->fit);
  wfit (x, y, w, n, link->fit);
  */
}

/* fit_combined_links: merge the minima from the even and odd
   samples, and fit a curve through the collected data.

   This code is currently defunct because when I collect a datum
   now, it goes into either evens or odds, AND it goes into both,
   so there is no reason to merge the data later.  */

void fit_combined_links (Link *both, Link *evens, Link *odds)
{
  int i;
  double x[MAX_SIZE];
  double y[MAX_SIZE];
  double w[MAX_SIZE];
  int n = 0;

  /* collect the message size and rtt for each sample */

  /* IMPORTANT NOTE: we never actually assemble the data from the
     even and odd samples into the "both" Link.  Only the "fit"
     object gets filled in; no data. */

  for (i=0; i<MAX_SIZE; i++) {
    if (evens->sample[i] != NULL) {
      assert (odds->sample[i] != NULL);

      x[n] = (double) i;
      if (evens->sample[i]->min < odds->sample[i]->min) {
	y[n] = evens->sample[i]->min;
      } else {
	y[n] = odds->sample[i]->min;
      }

      both->sample[i] = make_sample (both->ttl, i, 0);
      both->sample[i]->min = y[n];

      n++;
    }
  }
  fit (x, y, n, both->fit);
}

/* total_probes: count the total number of probes for a given link */

int total_probes (Link *link)
{
  int i;
  int n = 0;

  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      n += link->sample[i]->n;
    }
  }
  return n;
}

/* calc_link_difference: subtract adjacent links and fill-in
   slope and inter */

void calc_link_difference (Link *link, Link *prev,
			   double *slope, double *latency)
{
  Fit *f0, *f1;
  Fit fake[1];
  int i = link->ttl;

  f1 = link->fit;

  /* if there is no predecessor, make a fake Fit with zeros */
  if (i==0 || prev == NULL) {
    fake->b0 = 0.0;
    fake->b1 = 0.0;
    f0 = fake;
  } else {
    f0 = prev->fit;
  }

  *latency = (f1->b0 - f0->b0) / 2;
  *slope = f1->b1 - f0->b1;
}

/* convergence: calculate the convergence criterion for this
   link, by estimating the bandwith using all four adjacent
   pairs, and then finding the ratio of the spread to the
   median.  bw is an array of four doubles that
   gets filled in with the four bandwidth estimates, sorted
   from low to high */

double convergence (Path *path, int i, double *bw)
{
  int j;
  double diff;
  double slope[4], inter[4];
  double bw1, bw2;
  double median;

  /* calculate four differences between adjacent links */
  calc_link_difference (path->evens[i], path->evens[i-1],
			&slope[0], &inter[0]);
  calc_link_difference (path->odds[i], path->odds[i-1],
			&slope[1], &inter[1]);
  calc_link_difference (path->evens[i], path->odds[i-1],
			&slope[2], &inter[2]);
  calc_link_difference (path->odds[i], path->evens[i-1],
			    &slope[3], &inter[3]);

  /* sort the estimated slopes and calc the bandwidths */
  sort_array (slope, 4);
  for (i=0; i<4; i++) {
    bw[i] = 1 / slope[i] / 1000 * 8;
  }
  sort_array (bw, 4);
  
  median = (slope[1] + slope[2]) / 2;
  diff = fabs ((slope[3] - slope[0]) / median);
  return diff;
}

/* print_link_mins: if the -M option is on, we generate a file
   for each link that contains the SORTT for each packet size

   the names of the file are host.1.mins, host.2.mins, etc */

void print_link_mins (Link *link)
{
  int i;
  FILE *fp;
  char filename[100];

  sprintf (filename, "%s.%d.mins", host, link->ttl);
  fp = fopen (filename, "w");
  if (fp == NULL) {
    printf ("Unable to open output file %s.", filename);
    exit (1);
  }

  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      fprintf (fp, "%d\t%lf\n", i, link->sample[i]->min);
    }
  }
  fclose (fp);
}

/* print_sample_delays: subtract the fitted value from each
   probe in the sample and print the residual */

void print_sample_delays (FILE *fp, Sample *s, Fit *fit)
{
  int i;
  double res;

  for (i=0; i < s->n; i++) {
    res = s->x[i] - (fit->b0 + fit->b1 * s->size);
    if (res >= 0.0) {
      fprintf (fp, "%lf\n", res);
    }
  }
}

/* print_link_delays: subtract the fitted value from each
   probe in this link and print the residual in a series of files
   with the names host.1.delays, host.2.delays, etc. */

void print_link_delays (Link *link)
{
  int i;
  FILE *fp;
  char filename[100];

  sprintf (filename, "%s.%d.delays", host, link->ttl);
  fp = fopen (filename, "w");
  if (fp == NULL) {
    printf ("Unable to open output file %s.", filename);
    exit (1);
  }

  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      print_sample_delays (fp, link->sample[i], link->fit);
    }
  }
  fclose (fp);
}

/* print_link: print a one-line summary of the estimated
   characteristics for the given link of the given path */

void print_link (Path *path, int i)
{
  double latency, slope, bandwidth;
  Fit *f0, *f1;
  Fit fake[1];
  Link *link = path->both[i];
  Link *prev = path->both[i-1];
  int probes;
  double diff;
  double bw[4];

  if (dump_minima) {
    print_link_mins (path->both[i]);
  }

  probes = total_probes (path->both[i]);

  calc_link_difference (link, prev, &slope, &latency);
  bandwidth = 1 / slope / 1000 * 8;
  diff = convergence (path, i, bw) * 100;

  /* print human-readable format */

  printf ("|\tn=%5d  lat= %.3lf ms\tbw= (%0.3lf, %0.3lf, %0.3lf) Mb/s\n",
     probes, latency, bw[0], bandwidth, bw[3]);

  /* print LaTeX table format */
  /*  printf ("%2d   &  %4d   &  %.2lf  &  %.2lf  \\\% \\\\\n",
      i, probes, bandwidth, diff); */
}

/* print_sockaddr: try to look up the name of the given IP address,
   and print the name and address, or just the address */

void print_sockaddr (int ttl, Sockaddr *addr)
{
  int stat = -1;
  char str[NI_MAXHOST];

  if (dns_resolve) {
    stat = getnameinfo (addr, salen, str, sizeof(str), NULL, 0, 0);
  }

  if (stat == 0) {
    printf ("%d %s (%s)\n", ttl, str, Sock_ntop_host (addr, salen));
  } else {
    printf ("%d %s\n", ttl, Sock_ntop_host (addr, salen));
  }
}

/* print_host: print information about a host or intermediate router */

void print_host (Path *path, int ttl)
{
  Sockaddr *addr = &(path->addr[ttl]);
  print_sockaddr (ttl, addr);
}

/* print_path: iterate through the given path printing information
   about the links and routers */

void print_path (Path *path)
{
  int i;

  for (i=0; i<MAX_LINKS-1; i++) {
    if (path->both[i] != NULL) {
      print_link (path, i);
      print_host (path, i);
    }
  }
}

/* sample_size: return the number of probes for this link at this size */

int sample_size (Link *link, int size)
{
  if (link->sample[size] == NULL) return 0;
  else return link->sample[size]->n;
}

/* add_to_sample: adds a new datum to the given sample */

void add_to_sample (double rtt, Sample *s)
{
  if (s->n == s->space) {
    s->space *= 2;
    s->x = (double *) realloc ((void *) s->x, s->space * sizeof (double));
    assert (s->x != NULL);
  }
  s->x[s->n] = rtt;
  s->n++;
}

/* add_to_link: adds a new datum to the given link */

void add_to_link (int size, double rtt, Link *link)
{
  if (size >= MAX_SIZE) {
    fprintf (stderr, "Ignoring sample for message size %d\n", size);
    return;
  }
  if (link->sample[size] == NULL) {
    link->sample[size] = make_sample (link->ttl, size, MIN_SAMPLE_SIZE);
  }
  add_to_sample (rtt, link->sample[size]);
}

/* add_to_path: adds a new datum to the given path */

void add_to_path (int ttl, int size, double rtt, Path *path)
{
  int diff;

  if (ttl >= MAX_LINKS) {
    fprintf (stderr, "Ignoring sample for link %d\n", ttl);
    return;
  }

  /* the first time we see a new link, make all three structures */
  if (path->both[ttl] == NULL) {
    path->both[ttl] = make_link (ttl);
    path->evens[ttl] = make_link (ttl);
    path->odds[ttl] = make_link (ttl);
  }

  /* if the "evens" and "odds" have different numbers of
     samples, even them out; otherwise choose one at random */

  add_to_link (size, rtt, path->both[ttl]);

  diff = sample_size (path->evens[ttl], size) -
         sample_size (path->odds[ttl], size);

  if (diff < 0) { 
    add_to_link (size, rtt, path->evens[ttl]);
  } else if (diff > 0) {
    add_to_link (size, rtt, path->odds[ttl]);
  } else {
    diff = random_int (0, 2);
    if (diff == 0) {
      add_to_link (size, rtt, path->evens[ttl]);
    } else {
      add_to_link (size, rtt, path->odds[ttl]);
    }
  }
}

/* get_sockaddr: go to the path with index p and extract the
   address of the router with the given ttl */

Sockaddr *get_sockaddr (int p, int ttl)
{
  if (paths[p] == NULL) return NULL;

  if (paths[p]->evens[ttl] == NULL && paths[p]->odds[ttl] == NULL) {
    return NULL;
  } else {
    return &(paths[p]->addr[ttl]);
  }
}

/* get_sockaddr: for the path with index p and the given ttl,
   set the router address to the given Sockaddr */

void set_sockaddr (int p, int ttl, Sockaddr *addr)
{
  memcpy (&(paths[p]->addr[ttl]), addr, salen);
}

/* add_datum: use the information in the Datum structure to
   add a new probe to the big, hairy data structure */

void add_datum (int ttl, int size, Datum *datum)
{
  int p;
  Sockaddr *addr;
  int cmp;

  for (p=0; p<MAX_ALT_PATH; p++) {
    addr = get_sockaddr (p, ttl);
    if (addr == NULL) {
      if (paths[p] == NULL) {
	paths[p] = make_path ();
      }
      set_sockaddr (p, ttl, datum->addr);
      add_to_path (ttl, size, datum->rtt, paths[p]);
      return;
    } else {
      cmp = sock_cmp_addr (addr, datum->addr, salen);
      if (cmp == 0) {
	add_to_path (ttl, size, datum->rtt, paths[p]);
	return;
      } else {
	/* if it doesn't match, increase p and try again */
      }
    }
  }
}

/* residual: the difference between the SORTT in this sample
   and the fitted line */

double residual (Sample *s, Fit *fit)
{
  return s->min - (fit->b0 + fit->b1 * s->size);
}

/* remove_sample_min: kill the lowest rtt in the sample by
   setting it to a very high value.

   Note that this operation changes the collected data, which
   will mess things up if we try to calculate a moment or a distribution
   of queue delays */

void remove_sample_min (Sample *s)
{
  int i;
  int mindex = 0;

  /* find the minimal value and set it to a super-high value */
  for (i=0; i<s->n; i++) {
    if (s->x[i] < s->x[mindex]) mindex = i;
  }
  s->x[mindex] = 1e37;

  /* find the lowest remaining value so we can update min */
  mindex = 0;
  for (i=0; i<s->n; i++) {
    if (s->x[i] < s->x[mindex]) mindex = i;
  }

  /* change the minimum only if it is still a legitimate datum.
     if we have clobbered all the real data, leave the min alone */
  if (s->x[mindex] < 1e37) {
    s->min = s->x[mindex];
  }
}

/* remove_link_min: find the packet size with the most negative
   residual and clobber it */

void remove_link_min (Link *link)
{
  int i, mindex = 0;
  double res, minres;

  minres = 1e37;

  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      res = residual (link->sample[i], link->fit);
      if (res < minres) {
	minres = res;
	mindex = i;
      }
    }
  }
  if (minres < 1e37) {
    remove_sample_min (link->sample[mindex]);
  }
}

/* send_probe_sample: send a probe with the ttl and size taken
   from the given sample in the hope that the new datum will
   supplement the sample.  Keep in mind that this may not
   succeed, since the new probe may not go to the "right" router.

   Returns -1 if the attempt does not affect this sample
            0 if the new datum does not alter the sample minimum
	    1 if the new datum is a new minimum     */

int send_probe_sample (Sample *s)
{
  int n = s->n;
  int x = send_one_probe (s->size, s->ttl);

  if (n == s->n) return -1;

  if (s->x[s->n-1] < s->min) {
    s->min = s->x[s->n-1];
    return 1;
  }
  return 0;
}

/* send_probes_until: repeatedly calls send_probe_sample
   we get a new minimum (return 1), or we
   try max times (return 0) */

int send_probes_until (Sample *s, int max)
{
  int i, x;

  for (i=0; i<max; i++) {
    x = send_probe_sample (s);
    if (x == 1) return 1;
  }
  return 0;
}

/* get_new_minimum: 
   a return value of 1 means the sample minimum changed;
   0 means we sent n probes, but they either went to the
   wrong router or failed to set a new SORTT */

int get_new_minimum (Sample *s)
{
  return send_probes_until (s, s->n);
}

/* resample_high_spots: traverse the samples in this link and tries
   to find a new minimum for any sample above the fitted line.
   Returns the number of samples whose minimum changed. */

int resample_high_spots (Link *link)
{
  int i;
  double res, maxres;
  int ret;
  int count = 0;

  maxres = -1e37;

  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      res = residual (link->sample[i], link->fit);
      if (res > 0) {
	ret = get_new_minimum (link->sample[i]);
	if (ret > 0) count++;
      }
    }
  }
  return count;
}

/* process_link: perform whatever computation we need to do for
   this link, typically a least-squares fit and computing some
   other stats */

void process_link (Link *link)
{
  int i;

  if (link == NULL) return;

  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) {
      calc_order_stats (link->sample[i]);
    }
  }

  /* you have to calc_order_stats before you can perform a fit */
  fit_link (link);

  /*
  print_link_delays (link);
  link->delays = collect_link_delays (link);
  calc_moments (link->delays);
  */
}

/* fiddle_link: perform whatever operation we perform to get
   a link to change.  For example: we might do directed data
   collection or selective data removal. */

int fiddle_link (Link *link)
{
  int i, ret;
  int n = 0;

  if (link == NULL) return;

  /* get more data for any size with a positive residual */
  process_link (link);
  ret = resample_high_spots (link);

  /* see how many sizes we have measured */
  for (i=0; i<MAX_SIZE; i++) {
    if (link->sample[i] != NULL) n++;
  }

  /* remove data from the quarter of the samples with the lowest residuals */
  for (i=0; i<0; i++) {
    remove_link_min (link);
  } 
  
  return ret;
}

/* converge_link: invoke fiddle_link repeatedly until the convergence
   criterion for the given link reaches the threshhold (currently 0.1) */

void converge_link (Path *path, int ttl)
{
  int n = 2;
  int avail;
  int f1, f2, f3, f4;
  double diff;
  double bw[4];

  do {
    f1 = fiddle_link (path->evens[ttl]);
    f2 = fiddle_link (path->odds[ttl]);
    f3 = fiddle_link (path->evens[ttl-1]);
    f4 = fiddle_link (path->odds[ttl-1]);
    if (f1 < 0 && f2 < 0 && f3 < 0 && f4 < 0) break;

    process_link (path->evens[ttl]);
    process_link (path->odds[ttl]);
    process_link (path->evens[ttl-1]);
    process_link (path->odds[ttl-1]);
    diff = convergence (path, ttl, bw);

  } while (diff > 0.1);
}

/* process_path: traverse the path from from to to and process each link */

void process_path (Path *path, int from, int to)
{
  int i;
  Link *prev;
  Sample *s;

  if (path == NULL) {
    err_quit ("Tried to process a path that contains no data.");
  }

  for (i=0; i<MAX_LINKS; i++) {
    if (path->both[i] != NULL) {
      /* converge_link (path, i); */
    }
  }

  for (i=from; i<=to; i++) {
    if (path->both[i] != NULL) {
      process_link (path->evens[i]);
      process_link (path->odds[i]);
      process_link (path->both[i]);
    }
  }
}

// rand_uniform: returns 0.0 through 1.0, including 1.0 but not 0.0

double rand_uniform ()
{
  return 1.0 - (double) random() / ((double) RAND_MAX + 1);
}

/* rand_exp: random number from an exponential distribution
   with mean mu */

double rand_exp (double mu)
{
  return -log(rand_uniform()) * mu;
}

/* send_one_probe: send a probe with the given size and ttl,
   and adds the resulting datum to the big, hairy data structure.
   prints information about the probe on the screen and in the
   dump file (if requested)

   returns a count of the number of probes that returned
   a PORT UNREACHABLE error (indicating that they reached the
   destination) */

int send_one_probe (int size, int ttl)
{
  int code;
  Datum datum[1];
  int usecs;
  int count = 0;

  code = send_probe (size, ttl, wait_time, datum);

  probe_count++;
  if (code == -1) count++;

  if (dump_file_fp != NULL) {
    fprintf (dump_file_fp, "%d\t%s\t%d\t%.3lf\n",
	     ttl, Sock_ntop_host (datum->addr, salen), size, datum->rtt);
  }

  if (verbose) {
    printf ("ttl=%d\t%s\tsize= %4d B\trtt= %.3lf ms\n",
	      ttl, Sock_ntop_host (datum->addr, salen), size, datum->rtt);
  } else {
    printf ("                                                        \r");
    printf ("n=%5d\t%s\tsize= %4d B\trtt= %.3lf ms\r",
       probe_count, Sock_ntop_host (datum->addr, salen), size, datum->rtt);
  }
  fflush(stdout);

  add_datum (ttl, size, datum);

  if (inter_sample_time == 0.0) {
    usecs = (int) (datum->rtt * 1000 * inter_sample_ratio);
  } else {
    usecs = (int) (inter_sample_time * 1000);
  }
  usecs = rand_exp (usecs);
  microsleep (usecs);
  return count;
}

/* send_probes: send one probe at each of the given sizes, in
   a random order.  After each probe, sleep for an amount of
   time equal to the rtt of the previous probe.  Returns the
   number of probes that yielded a PORT UNREACHABLE error,
   indicating that they reached the destination. */

int send_probes (Sizes *s, int ttl)
{
  int i;
  int index, size;
  int done_count = 0;

  permute_sizes (s);

  for (i=0; i<s->n; i++) {
    index = s->permute[i];
    size = s->size[index];
    done_count += send_one_probe (size, ttl);

  }
  return done_count;
}

/* read_data: fp is a file pointer to a dump file that is open
   for reading.  Read the data from the file and add each probe
   to the big, hairy data structure */

void read_data (FILE *fp)
{
  int ttl, size;
  char ip_addr[128];
  Datum datum[1];
  double rtt;
  Sockaddr_in *addr = (Sockaddr_in *) &(datum->addr);
  int n;

  while (1) {
    n = fscanf (fp, "%d %s %d %lf", &ttl, ip_addr, &size, &rtt);
    if (n != 4) break;

    n = convert_sockaddr (ip_addr, addr, salen);
    if (n <= 0) {
      err_quit ("Invalid IP address in dump file: %s\n", ip_addr);
    }
    datum->rtt = rtt;
    add_datum (ttl, size, datum);
  }
}

/* process_file: open the dump file, read the data, process
   the resulting path, and print the results. */

void process_file ()
{
  FILE *input_file_fp;

  input_file_fp = fopen (input_file, "r");
  if (input_file_fp == NULL) {
    err_quit ("Unable to open input file %s", input_file);
  }

  read_data (input_file_fp);
  process_path (paths[0], 0, MAX_LINKS);

  printf ("0 localhost\n");
  print_path (paths[0]);
}

/* count_probes: count the number of probes that have been sent
   to this link at this size. */

int count_probes (Link *link, int size)
{
  assert (link != NULL);
  if (link->sample[size] == NULL) return 0;
  return link->sample[size]->n;
}

/* check_link: check and see if any of the alternate paths with
   the given ttl have accumulated the requested number of probes
   (num_probes) at each of the desired packet sizes

   If so, return the index of the first path that passes the test.
   If not, return -1  */

int check_link (Sizes *s, int ttl) {
  int i, k, n;
  int size;
  int count;

  for (k=0; k<MAX_ALT_PATH; k++) {
    if (paths[k] == NULL) continue;
    count = 0;

    for (i=0; i<s->n; i++) {
      size = s->size[i];
      n = count_probes (paths[k]->both[ttl], size);
      if (n >= num_probes) count++;
    }    
    if (count == s->n) return k;
  }
  return -1;
}

/* fill_in_link: traverse each of the alternate paths with the
   given ttl, and for each one, traverse the sizes array in a
   random order looking for sizes that have not been adequately
   sampled (n < num_probes).

   If you find any, request a probe with the given ttl and size,
   keeping in mind that it may or may not follow the "right" path,
   and therefore may not actually fill in the hole that was
   discovered.

   If the requested probe succeeds in rendering one of the
   alternate paths complete, return the index of the winner.
   Otherwise, return -1 */

int fill_in_link (Sizes *s, int ttl) {
  int i, k, n;
  int size;
  int count = 0;
  int index;

  for (k=0; k<MAX_ALT_PATH; k++) {
    if (paths[k] == NULL) continue;
    permute_sizes (s);

    for (i=0; i<s->n; i++) {
      index = s->permute[i];
      size = s->size[index];
      n = count_probes (paths[k]->both[ttl], size);
      if (n < num_probes) {
	if (verbose) {
	  printf ("filling in path= %d, size= %d, n= %d\n", k, size, n);
	}
	send_one_probe (size, ttl);
      }
    }
    index = check_link (s, ttl);
    if (index != -1) return index;
  }
  return -1;
}

/* swap_link: swap all the data at the given ttl from two
   alternate paths (with indices k1 and k2) */

void swap_link (int ttl, int k1, int k2)
{
  Link *temp;
  Sockaddr sa[1];

  temp = paths[k1]->evens[ttl];
  paths[k1]->evens[ttl] = paths[k2]->evens[ttl];
  paths[k2]->evens[ttl] = temp;

  temp = paths[k1]->odds[ttl];
  paths[k1]->odds[ttl] = paths[k2]->odds[ttl];
  paths[k2]->odds[ttl] = temp;

  temp = paths[k1]->both[ttl];
  paths[k1]->both[ttl] = paths[k2]->both[ttl];
  paths[k2]->both[ttl] = temp;

  memcpy (sa, &paths[k1]->addr[ttl], sizeof (Sockaddr));
  memcpy (&paths[k1]->addr[ttl], &paths[k2]->addr[ttl], sizeof (Sockaddr));
  memcpy (&paths[k2]->addr[ttl], sa, sizeof (Sockaddr));
}

/* clink_loop: traverse the relevant values of ttl (from min_ttl
   to max_ttl) and send out the requested number of probes at
   each of the requested sizes.

   Check and see if all the probes went to the same router, in
   which case we're done, or whether we have to send more probes
   before one of the alternate paths gets "topped off".

   Once enough probes are sent, process the new data and print
   a summary of it */

void clink_loop ()
{
  int i, ttl, k;
  int index;
  int done_count;
  Sizes *sizes;

  if (kernel_timestamps) {
    tos = IPTOS_LOWDELAY | IPTOS_THROUGHPUT | IPTOS_MINCOST;
  }

  clink_init (host, tos);

  sizes = make_sizes (low, high, step);

  printf ("  %d probes at each of %d sizes (%d to %d by %d)\n",
	  num_probes, sizes->n, low, high, step);
  
  printf ("0 localhost\n");

  for (ttl=min_ttl; ttl<=max_ttl; ttl++) {
    probe_count = 0;                  /* how many probes we've sent */
    done_count = 0;                   /* how many PORT UNREACH errors */

    for (i=0; i<num_probes; i++) {
      done_count += send_probes (sizes, ttl);
    }

    while ( (index = check_link (sizes, ttl)) == -1) {
      index = fill_in_link (sizes, ttl);
      if (index != -1) break;
    }

    if (index != 0) {
      swap_link (ttl, 0, index);
    }

    printf ("                                                       \r");

    process_path (paths[0], ttl, ttl);
    print_link (paths[0], ttl);
    print_host (paths[0], ttl);

    if (done_count > 0) break;
  }
}

/* main: process the options and call the appropriate procedures */

main (int argc, char *argv[])
{
  int i, c;

  opterr = 0;
  while ( (c = getopt (argc, argv, "knvMD:I:f:h:i:l:m:q:r:s:t:w:")) != -1) {
    switch (c) {
    case 'k':
      kernel_timestamps = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'n':
      dns_resolve = 0;
      break;
    case 'D':
      dump_file = optarg;
      dump_file_fp = fopen (dump_file, "w");
      if (dump_file_fp == NULL) {
	printf ("Unable to open output file %s. Ignoring -D option.",
		dump_file);
      }
      break;
    case 'I':
      input_file = optarg;
      break;
    case 'M':
      dump_minima = 1;
      break;
    case 'f':
      if ( (min_ttl = atoi(optarg)) < 1) {
	err_quit ("invalid -f value (min_ttl)");
      }
      break;
    case 'm':
      if ( (max_ttl = atoi(optarg)) < 1) {
	err_quit ("invalid -m value (max_ttl)");
      }
      break;
    case 'h':
      if ( (high = atoi(optarg)) < 1) {
	err_quit ("invalid -h value (high)");
      }
      break;
    case 'l':
      if ( (low = atoi(optarg)) < 1) {
	err_quit ("invalid -l value (low)");
      }
      break;
    case 'i':
      if ( (inter_sample_time = atof(optarg)) < 0.0) {
	err_quit ("invalid -i value (inter_sample_time)");
      }
      break;
    case 'r':
      if ( (inter_sample_ratio = atof(optarg)) < 0.0) {
	err_quit ("invalid -r value (inter_sample_ratio)");
      }
      break;
    case 's':
      if ( (step = atoi(optarg)) < 1) {
	err_quit ("invalid -s value (step)");
      }
      break;
    case 't':
      if ( (tos = atoi(optarg)) < 1) {
	err_quit ("invalid -t value (tos)");
      }
      if (tos > 31) {
	err_quit ("invalid -t value (tos)");
      }
      break;
    case 'q':
      if ( (num_probes = atoi(optarg)) < 2) {
	err_quit ("invalid -q value (num_probes); must be at least 2");
      }
      break;
    case 'w':
      if ( (wait_time = atoi(optarg)) < 1) {
	err_quit ("invalid -w value (wait_time)");
      }
      break;
    case '?':
      err_quit ("unrecognized option");
    }
  }

  for (i=0; i<MAX_ALT_PATH; i++) {
    paths[i] = NULL;
  }

  if (input_file != NULL) {
    host = input_file;
    process_file ();
    exit (0);
  }

  if (optind != argc - 1) {
    err_quit ("usage: clink [options] <hostname> or clink -Iinput_file");
  }

  host = argv[optind];

  clink_loop ();

  return 0;
}
