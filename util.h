#include "clink.h"
#include <unistd.h>

#define MAXLINE 4096

void err_doit (int errnoflag, int level, const char *fmt, va_list ap);;
void err_sys (char *fmt, ...);;
void err_quit (char *fmt, ...);
void err_msg (char *fmt, ...);
char *sock_ntop_host(const struct sockaddr *sa, socklen_t salen);
char *Sock_ntop_host(const struct sockaddr *sa, socklen_t salen);
void sock_set_port(struct sockaddr *sa, socklen_t salen, int port);
int sock_cmp_addr(const struct sockaddr *sa1, const struct sockaddr *sa2,socklen_t salen);
void tv_sub (struct timeval *out, struct timeval *in);
char *icmpcode_v4(int code);
Sigfunc *signal(int signo, Sigfunc *func);
Sigfunc *Signal(int signo, Sigfunc *func)	/* for our signal() function */;
void *Malloc(size_t size);
void *Calloc(size_t n, size_t size);
void Gettimeofday(struct timeval *tv, void *foo);
void Pipe(int *fds);
void Bind(int fd, const struct sockaddr *sa, socklen_t salen);
void Setsockopt(int fd, int level, int optname, const void *optval,		socklen_t optlen);
struct addrinfo *host_serv (char *host, char *serv, int family, int socktype);
struct addrinfo *Host_serv(const char *host, const char *serv, int family, int socktype);
ssize_t Read(int fd, void *ptr, size_t nbytes);
void Write(int fd, void *ptr, size_t nbytes);
ssize_t Recvfrom(int fd, void *ptr, size_t nbytes, int flags,		 struct sockaddr *sa, socklen_t *salenptr);
int Socket(int family, int type, int protocol);
void Sendto(int fd, const void *ptr, size_t nbytes, int flags,	   const struct sockaddr *sa, socklen_t salen);
int convert_sockaddr (char *ip_addr, Sockaddr_in *addr, socklen_t salen);
void microsleep (int usecs);
