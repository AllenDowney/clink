#include "clink.h"

extern int kernel_timestamps;
extern int dns_resolve;
extern int verbose;
extern socklen_t salen;

/* other global variables */

char recvbuf[MAX_SIZE];
char sendbuf[MAX_SIZE];

int sendfd, recvfd;
int pipefd[2];              /* the pipe for the alarm handler */

Sockaddr *sasend;    /* socket addresses for various purposes */
Sockaddr *sarecv;
Sockaddr *sabind;

#define MIN_DPORT (32768 + 668)     /* 668 = the neighbor of the beast */
#define MAX_DPORT 65535

u_short dport = MIN_DPORT;          /* destination port cycles through
				       the high numbers, might get
				       unlucky occasionally */
u_short sport;                      /* source UDP port # (set in main) */

Timeval sendtv[1];
Timeval recvtv[1];
Timeval difftv[1];
Timeval kdifftv[1];
Timeval senddiff[1];
Timeval recvdiff[1];

Timeval *ksendtv;
Timeval *krecvtv;

/* NOTE: system calls beginning with a capital letter are Stevens's
   wrapper functions.  Each one invokes the method and checks the
   return value.  If the call fails, it invokes err_sys, which prints
   the error message and quits.

   Types that begin with a capital letter are usually typedefs
   that I added because (1) I hate seeing the word struct all over
   the place, and (2) it lets me pretend I am writing Java. */

/* sig_alrm: alarm handler sends a message to the process through
   a pipe, causing select to return */

void sig_alrm (int signo)
{
  Write (pipefd[1], "", 1);  /* write 1 null byte to pipe */
  return;
}

/* extract_timeval: in the kernel-hacked version, the returning
   timestamp is hidden in the src and dest fields of the IP header
   of the enclosed packet.  */

void extract_timeval (struct ip *ip)
{
  krecvtv = (Timeval *) &ip->ip_src;
}

/* process_ip: extracts all the info from an incoming ICMP packet

     Just for kicks, I changed the BSD-style names of the ICMP
     errors to Linux-style names, mostly so that they will be
     consistent with the changes I made in the kernel and so my
     head won't explode. */

int process_ip (struct ip *ip, int len)
{
  int hlen1, hlen2, icmplen;
  struct icmp *icmp;
  struct ip *hip;
  struct udphdr *udp;

  hlen1 = ip->ip_hl << 2;                        /* length of IP header */
  icmp = (struct icmp *) (recvbuf + hlen1);
  icmplen = len - hlen1;

  if (verbose) {
    printf ("received ICMP type= %d  code= %d\n",
	    icmp->icmp_type, icmp->icmp_code);
  }

  if (icmplen < 8 + 20 + 8) return 0;

  if (icmp->icmp_type != ICMP_TIME_EXCEEDED &&
      icmp->icmp_type != ICMP_DEST_UNREACH)
    return 0;

  /* hip is the header of the enclosed IP packets, supposedly
     the header of the packet that caused the error */

  hip = (struct ip *) (recvbuf + hlen1 + 8);
  if (hip->ip_p != IPPROTO_UDP) return 0;

  hlen2 = hip->ip_hl << 2;
  udp = (struct udphdr *) (recvbuf + hlen1 + 8 + hlen2);

  if (udp->source != htons (sport)) return 0;
  if (udp->dest != htons (dport)) return 0;

  /* now we know it's an ICMP packet caused by a UDP
     datagram sent by us and sent to the port we happen to
     be sending to.  It's probably one of ours. */

  extract_timeval (hip);

  if (icmp->icmp_type == ICMP_TIME_EXCEEDED) {
    if (icmp->icmp_code == ICMP_EXC_TTL) {
      return -2;
    } else {
      return 0;
    }
  }

  if (icmp->icmp_type == ICMP_DEST_UNREACH) {
    if (icmp->icmp_code == ICMP_PORT_UNREACH) {
      return -1;
    } else {
      return 0;
    }
  }
}

/* recv_dgram: reads all incoming datagrams and checks for
   returning ICMP packets.
   returns -3 on timeout,
           -2 on ICMP time exceeded in transit (we reached a router)
	   -1 on ICMP port unreachable (we reached the destination)
	    0 on ICMP that has nothing to do with us  */

  /* as Stevens points out in Section 18.5 of Unix Network Programming,
     many programs with alarms have a race condition, which is that
     the alarm might go off before we get to the recvfrom, in which
     case it does nothing and the recvfrom might wait indefinitely.

     In earlier versions of this code, this problem seemed to pop
     up occasionally (although I am not positive about that).

     The use of select here solves that problem.  When the alarm
     goes off, the alarm handler sends a message through the pipe,
     which is one of the things select waits for.

     When select returns, we know that we have received a datagram
     OR the alarm has gone off OR both.  We then use rset to find
     out which, and deal with it.

     According to the specification of select, it should not be possible
     to get to the recvfrom unless there is a datagram waiting, and
     therefore the recvfrom should never block.  Nevertheless, it sometimes
     does, which is why, when we opened it, we set the NONBLOCK flag
     and why, if it fails (errno = EAGAIN) we just go on. */

int recv_dgram (int timeout)
{
  int err;
  socklen_t len;
  ssize_t n;
  struct ip *ip;
  int maxfdp1 = max (recvfd, pipefd[0]) + 1;
  fd_set rset[1];  
  FD_ZERO (rset);

  alarm(timeout);    /* set the timeout alarm to handle dropped packets */

  while (1) {
    FD_SET (recvfd, rset);
    FD_SET (pipefd[0], rset);

    n = select (maxfdp1, rset, NULL, NULL, NULL);
    if (n < 0 && errno != EINTR) {
      err_sys ("select error");
    }

    if (FD_ISSET (recvfd, rset)) {
      len = salen;
      n = recvfrom (recvfd, recvbuf, sizeof(recvbuf), 0, sarecv, &len);
      err = errno;
      Gettimeofday (recvtv, NULL);   /* get time of packet arrival */
      if (n < 0 && err != EAGAIN) {
	err_sys ("recvfrom error");
      }
    }

    if (FD_ISSET (pipefd[0], rset)) {
      Read (pipefd[0], &n, 1);

      if (verbose) {
	printf ("timeout\n");
      }
      return -3;                 /* timeout */
    }

    ip = (struct ip *) recvbuf;
    return process_ip (ip, n);
  }
}

/* sub_tv: subtract minus from plus and put the result in res */

void sub_tv (Timeval *plus, Timeval *minus, Timeval *res)
{
  res->tv_sec = plus->tv_sec - minus->tv_sec;
  res->tv_usec = plus->tv_usec - minus->tv_usec;

  if (res->tv_usec < 0) {
    res->tv_sec--;
    res->tv_usec += 1000000;
  }
}

/* time_to_double: convert a Timeval to a double.  This only
   works with Timevals that are small (like the difference between
   two real Timevals) */

double time_to_double (Timeval *time)
{
  return time->tv_sec * 1000.0 + time->tv_usec / 1000.0;
}

void fill_in_datum (Datum *datum)
{
  /* stat = sock_cmp_addr (sarecv, salast, salen); */

  /* copy the address of the router into the return structure */
  memcpy (datum->addr, sarecv, salen);

  /* calculate the round trip time using user-level timestamps */
  sub_tv (recvtv, sendtv, difftv);
  datum->rtt = time_to_double (difftv);

  if (kernel_timestamps) {
    /* check and see if the kernel-level timestamps make sense */
    sub_tv (ksendtv, sendtv, senddiff);
    sub_tv (recvtv, krecvtv, recvdiff);

    if (senddiff->tv_sec < 0 || senddiff->tv_sec > 2) {
      err_sys ("\nI don't think the outgoing kernel timestamps are working\n");
    }

    if (recvdiff->tv_sec < 0 || recvdiff->tv_sec > 2) {
      err_sys ("\nI don't think the incoming kernel timestamps are working\n");
    }

    sub_tv (krecvtv, ksendtv, kdifftv);
    datum->krtt = time_to_double (kdifftv);
  } else {
    datum->krtt = 0.0;
  }
}

/* send_dgram: generate an outgoing UDP packet

   we add the sequence number to the destination port so
   that when we get the ICMP packet back, we can make sure
   the sequence number is right

   the second effort send is a kludge to handle a funny
   thing, which is that the last host seems to refuse the
   second or third connection consistently, which might
   might mean that something breaks when we get the
   ICMP_DEST_UNREACH error.  The second attempt seems
   to succeed consistently. */

void send_dgram (int datalen, int ttl)
{
  int n;

  dport++;
  if (dport >= MAX_DPORT) dport = MIN_DPORT;
  if (dport < MIN_DPORT) dport = MIN_DPORT;
  sock_set_port (sasend, salen, htons(dport));

  if (verbose) {
    printf ("sending datagram to port %d\n", dport);
  }

  //strcpy (sendbuf, "<clink!>");
  Gettimeofday (sendtv, NULL);
  Sendto (sendfd, sendbuf, datalen, 0, sasend, salen);
}

/* send_probe: sends a probes with the given size and ttl and
   then waits for the reply.  

   If we want the total packet size to be _size_, we have to
   make the datalen = size - 28 (IP header = 20 bytes, UDP
   header = 8 */

int send_probe (int size, int ttl, int timeout, Datum *datum) 
{
  int code;
  int datalen = size - 28;
  int flag = 1;

  if (datalen < 0) {
    err_sys ("packet size too small, has to be at least 28\n");
  }

  if (datalen > MAX_SIZE) {
    err_sys ("packet size too big, MAX_SIZE = %d\n", MAX_SIZE);
  }

  Setsockopt (sendfd, IPPROTO_IP, IP_TTL, &ttl, sizeof(int));

  /* the following option is only necessary on Linux machines
     because they have the unusual behavior of returning some ICMP
     errors to the process that generated them, which means that
     when we get a PORT UNREACHABLE error, it causes the next
     sendto to fail and return ECONNREFUSED.  Switching to
     BSD_COMPAT mode solves that problem.

     See Stevens, UNPv12e Section 8.9 Server Not Running for
     an explanation of weird UDP asynchronous errors.
 */

  Setsockopt (sendfd, SOL_SOCKET, SO_BSDCOMPAT, &flag, sizeof(int));

  do {
    send_dgram (datalen, ttl);

    do {
      code = recv_dgram (timeout);
    } while (code == 0);

    /* if we get an irrelevant ICMP, go back and get another */

  } while (code == -3);

  /* if we get a timeout, send another datagram */

  fill_in_datum (datum);
  return code;
}

/* clink_init :

   The weird TOS options are there
   as a signal to the kernel to identify clink packets so it can
   fill in the timestamps.  I am assuming that they don't have
   any actual effect.  */

void clink_init (char *host, int tos)
{
  int i, n;
  struct addrinfo *ai = NULL;

  if (dns_resolve) {
    ai = Host_serv (host, NULL, 0, 0);
    printf ("clink to %s (%s)\n",
	    ai->ai_canonname,
	    Sock_ntop_host (ai->ai_addr, ai->ai_addrlen));

    if (ai->ai_family != AF_INET) {
      err_quit ("unknown address family %d", ai->ai_family);
    }

    sasend = ai->ai_addr;
  } else {
    sasend = Calloc (1, salen);
    n = convert_sockaddr (host, (Sockaddr_in *) sasend, salen);
    if (n <= 0) {
      err_msg ("clink is unable to resolve this address: %s", host);
      err_msg ("The address may be ill-formed, or it may be failing");
      err_msg ("because you are using the -n option, which prevents");
      err_quit ("clink from using DNS to resolve this address.");
    }
  }

  sarecv = Calloc (1, salen);
  sabind = Calloc (1, salen);

  Pipe (pipefd);     /* the pipe for the alarm handler */

  recvfd = socket (sasend->sa_family, SOCK_RAW, IPPROTO_ICMP);
  if (recvfd == -1) {
    if (errno == EPERM) {
      printf ("\nclink was unable to open a raw socket.  The most\n");
      printf ("likely cause is that you are not running it as root.\n");
      exit (1);
    } else {
      err_sys ("could not open raw socket");
    }
  }

  n = fcntl (recvfd, F_SETFL, O_NONBLOCK);
  if (n != 0) {
    err_sys ("fcntl could not make socket nonblocking");
  }

  /* if we are running setuid, we can stop now */
  setuid (getuid ());

  /* create the outgoing socket and set the funky TOS */
  sendfd = Socket (sasend->sa_family, SOCK_DGRAM, 0);
  Setsockopt (sendfd, IPPROTO_IP, IP_TOS, &tos, sizeof(int));

  ksendtv = (Timeval *) sendbuf;

  sabind->sa_family = sasend->sa_family;
  sport = (getpid() & 0xffff) | 0x8000;       /* source UDP port # */
  sock_set_port (sabind, salen, htons(sport));
  Bind (sendfd, sabind, salen);

  Signal (SIGALRM, sig_alrm);

  /* there is some possibility that the contents of the message
     will affect its transmission speed on links that do some kind
     of zero-balancing */
  for (i=0; i<MAX_SIZE; i++) {
    sendbuf[i] = 'U';
  }
}




