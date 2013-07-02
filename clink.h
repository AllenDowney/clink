#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <time.h>
#include <netinet/in.h>
#include <netinet/in_systm.h>
#include <netinet/ip.h>
#include <netinet/ip_icmp.h>
#include <netinet/udp.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/un.h>
#include <stdarg.h>
#include <syslog.h>

/* maximum message size */
#define MAX_SIZE 3000

/* default number of probes initially allocated for a sample 
   (doesn't matter much, since samples expand dynamically) */
#define MIN_SAMPLE_SIZE 8

/* maximum number of links in a path */
#define MAX_LINKS 30

typedef struct timeval Timeval;
typedef struct sockaddr Sockaddr;
typedef struct sockaddr_in Sockaddr_in;

typedef struct datum {
  double rtt, krtt;
  Sockaddr addr[1];
} Datum;

typedef struct sizes {
  int *size;
  int *permute;
  int n;
} Sizes;

typedef struct {
  double *x;                       /* data */
  int n, space;                    /* number of items and amount of space */
  int ttl, size;
  double min, lep;
  double zeta, sig2;
  double m1p, m2p;
  double m1, m2, m3;               /* moments about the mean */
} Sample;

typedef struct {
  double b0, b1;
  double sig0, sig1;
  double chi2, r2;
} Fit;

typedef struct {
  int ttl;
  Sample *sample[MAX_SIZE];
  Fit fit[1];
  Sample *delays;
} Link;

typedef struct {
  Link *evens[MAX_LINKS];
  Link *odds[MAX_LINKS];
  Link *both[MAX_LINKS];
  Sockaddr addr[MAX_LINKS];
} Path;

/* the following are a few definitions from Stevens' unp.h */

typedef	void Sigfunc(int);        /* for signal handlers */

#define max(a,b) ((a) > (b) ? (a) : (b))

/* the following are prototypes for the Stevens utilities in util.c */

char *Sock_ntop_host(const struct sockaddr *sa, socklen_t salen);
void sock_set_port(struct sockaddr *sa, socklen_t salen, int port);
int sock_cmp_addr(const struct sockaddr *sa1,
		  const struct sockaddr *sa2, socklen_t salen);
void tv_sub (struct timeval *out, struct timeval *in);
char *icmpcode_v4(int code);
Sigfunc *Signal(int signo, Sigfunc *func);
void *Calloc(size_t n, size_t size);
void Gettimeofday(struct timeval *tv, void *foo);
void Pipe(int *fds);
void Bind(int fd, const struct sockaddr *sa, socklen_t salen);
void Setsockopt(int fd, int level, int optname, const void *optval,
		socklen_t optlen);
void Sendto(int fd, const void *ptr, size_t nbytes, int flags,
	    const struct sockaddr *sa, socklen_t salen);
struct addrinfo *Host_serv(const char *host, const char *serv,
                           int family, int socktype);
ssize_t Read(int fd, void *ptr, size_t nbytes);
void Write(int fd, void *ptr, size_t nbytes);
ssize_t Recvfrom(int fd, void *ptr, size_t nbytes, int flags,
		 struct sockaddr *sa, socklen_t *salenptr);
int Socket(int family, int type, int protocol);
void err_sys (char *fmt, ...);
void err_quit (char *fmt, ...);
void err_msg (char *fmt, ...);
void microsleep (int usecs);
int convert_sockaddr (char *ip_addr, Sockaddr_in *addr, socklen_t salen);



