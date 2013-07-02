CLINK -- Characterize LINKs
---------------------------

Clink is a utility that estimates the latency and bandwidth of
Internet links by sending UDP packets from a single source and
measuring round-trip times.  The basic mechanism is similar to
ping and traceroute, except that clink generally has to send
many more packets.

Clink was written by Allen Downey (downey@colby.edu), based on
pathchar, a similar program written by Van Jaconbson.  The interface
of clink is based on the interface of pathchar, and the underlying
mechanism is based on Jacobson's description of pathchar.  No pathchar
source code is included in clink.

Clink is based on trout, which is a simple version of traceroute
written by Allen Downey, but which is heavily based on the version of
traceroute Richard Stevens presents in his book UNIX Network
Programming: Volume 1, Second Edition.

For Copyright information, see the file named COPYRIGHT that
should have been included with this distribution or contact:

	Allen Downey
	Colby College
	5850 Mayflower Hill
	Waterville, ME 04901
	downey@colby.edu
	http://www.cs.colby.edu/~downey



Information about pathchar and clink:

There is a web page at CAIDA that describes pathchar:

	http://www.caida.org/Pathchar/

At this time (19 July 1999), an alpha version of pathchar
is available in binary form for FreeBSD, Linix, Solaris and
some other platforms:

	ftp://ftp.ee.lbl.gov/pathchar/

There is also a set of slides there that Jacobson used to
present pathchar at a conference.

In Summer 1998, I wrote the data processing part of clink,
which became process.c.  I used the alpha version of pathchar
to collect data, and wrote a paper describing my experiences,
and proposing some improvements to pathchar's data processing.
That paper will appear at SIGCOMM in September, 1999.  A
Postscript version of the paper is available from

	http://www.cs.colby.edu/~downey/pathchar

In Summer 1999, I wrote the data collection part of clink,
called collect.c, and assembled the pieces into this distribution.
I am currently working on (1) implementing the improvements I
suggested in my SIGCOMM paper, and (2) testing additional
improvements.

So far, the primary differences between pathchar and clink are:

1) clink uses the even-odd technique described in the SIGCOMM
	paper to generate interval estimates for bandwidth.

2) when clink encounters a routing instability, it collects data
	for ALL the paths it encounters, until one of the paths
	generates enough data to yield an estimate.

The next major addition to clink will be adaptive data collection
as described in the SIGCOMM paper.



HOW TO COMPILE CLINK
--------------------

To unpack the distribution:

	tar -xzf clink.1.0.tar.gz    (if the tar file is compressed)
or
	tar -xf clink.1.0.tar        (if the tar file is not compressed)

Move into the newly created directory and compile clink:
	
	cd clink.1.0
	make

Clink will probably only compile and run on Linux machines, but
it should be possible to port it to other UNIX platforms without
an enormous amount of work.

Because clink has to open a raw socket, you must either run it as root
or make the executable setuid root.  The latter is preferable; here's
how:

	su root
	(type the root password)
	chmod +s clink
	exit

If you try to run clink as a mere mortal, you will get a message like:

	clink was unable to open a raw socket.  The most
	likely cause is that you are not running it as root.


HOW TO USE CLINK
----------------

Basic invocation
----------------

To characterize the links along a path between you and another host,
type

	clink host_name

If clink succeeds, the output will look something like:

clink to host-06.colby.edu (137.146.210.39)
  8 probes at each of 93 sizes (28 to 1500 by 16)
0 localhost
|       n=  744  lat= 0.234 ms  bw= (8.493, 8.506, 8.506) Mb/s
1 cb3500-02.switches.network.colby.edu (137.146.194.17)
|       n=  744  lat= 0.182 ms  bw= (9.034, 9.134, 9.552) Mb/s
2 port-0.router-1.network.colby.edu (137.146.238.209)
|       n=  744  lat= 1.020 ms  bw= (4.325, 4.389, 4.458) Mb/s
3 host-06.colby.edu (137.146.210.39)

n is the number of probes that were used to characterize each
link.  In this case, clink makes 8 measurements at each of 93
sizes, for a total of 744 links.  If clink encounters a routing
instability, it may have to send more probes before it gets a
complete set of probes at each size.  If you encounter an alternating
link, you might want to use the -D option to generate a dump
file, and then examine the dump file for more information about
the instability.

lat indicates latency, in milliseconds.  bw indicates bandwidth, in
megabits per second.  Three values are given for bandwidth: a low
estimate, a high estimate, and (in the middle) a best estimate.  The
distance between the high and low estimates gives some indication of
how reliable the estimate is.  For reasons that are explained in the
SIGCOMM paper, the "best" estimate does not necessarily fall between
the high and low values.


Controlling the starting and ending links
-----------------------------------------

If you want to start at something other than the first link,
use the -f option:

	clink -f2 host_name

In this case, clink starts with the ttl set to 2, meaning that the
packet will go two hops before returning.  Of course, if the starting
ttl is greater than 1, the estimated characteristics for the
first "link" will really be the aggregate characteristics of the first
few links.

To stop before you get to the destination, use the -m option:

	clink -m3 host_name

This example would characterize only the first three links
of the path to host_name.


Controlling the number of probes
--------------------------------

To set the number of probes sent over each link at each size,
use the -q option:

	clink -q16 host_name

This example would send 16 probes instead of the default (8).

The -l, -h, and -s option control the range of packet sizes clink
uses, by specifying the low size, the high size, and the step,
respectively.  The default values are 28, 1500 and 16, indicating
sizes from 28 to 1500 by steps of 16 (yielding 93 different values).

These sizes include IP and UDP headers, hence the smallest possible
value is 28 (a UDP packet with no data).  1500 is a common MTU
(maximum transfer unit) for Ethernets.  If the range of sizes spans
the MTU of one of the links, it messes up clink's estimates.

The packet sizes do not include Ethernet headers, which means
that the actual size while the packets are on Ethernet links is
a little bigger (8 bytes bigger, I think).  This omission does
not affect the bandwidth estimates, but it does perturb the latency
estimates of some links just a little.


Controlling the rate of probing
-------------------------------

clink works by sending out a lot of UDP packets that contain
deliberate errors (like a ttl too small to get to the destination or a
port number that is unlikely to be in use at the destination).  These
probes impose a load on the network and on the destination machine.
In fact, sending a barrage of probes to an unsuspecting machine can be
interpreted as an act of hostility.

Users of clink should use caution to avoid overloading networks or
harassing machines at other sites.  By default, clink only sends out
one probe at a time, and it waits a while between probes.
Specifically, it waits for a period equal to 10 times the round trip
time of the previous packet.  Thus, in the worst case, clink-
generated traffic uses 1/11th of the network's capacity.

You can control the wait time between packers with the -r and -i
options.  The -r option changes the ratio between the wait time and
the round trip time.

	clink -r12 hostname

This example waits for a period 12 times the round trip time of 
each packet before sending the next.

The -i option specifies the intersample time in milliseconds,
independent of the round trip time.

	clink -i100 hostname

This example waits 100ms between samples.


Generating dump files
---------------------

The -D option generates a file that contains the raw data
clink collects.  You can specify the name of the dump file.

	clink -Dhost.dump host_name

The file host.dump will look like this:

1       137.146.194.17  876     1.303
1       137.146.194.17  1420    1.812
1       137.146.194.17  1340    1.743
1       137.146.194.17  1100    1.505

The first column is the ttl (how many links the packet traversed).
The second column is IP address of the router or host that generated
the ICMP packet.  The third column is the packet size in bytes.  clink
sends the various sizes out in random order.  The fourth column is the
round trip time in milliseconds.

clink also provides an option that prints the SORTT (shortest
observed round trip time) for each packet size sent to each link.
The -M option causes clink to write one file for each link,
with the names hostname.0.mins, hostname.1.mins, etc.  The format
of these file is

28      0.515000
44      0.515000
60      0.529000
76      0.545000
92      0.557000

The first column is the packet size, the second column is the
SORTT in milliseconds.


Reading from a dump file
------------------------

To generate link estimates from a dump file, use the -I option:

	clink -Ihost.dump

The output is identical to the output that was generated when
the data was collected.  This feature is probably not useful for
clink users, but it is very useful for trying out variations on
clink's estimation methods.


Avoiding DNS resolution
-----------------------

Normally clinks uses DNS to find the IP address the destination
machine and the names of the intermediate routers.  Because DNS
resolution can be slow (or unavailable), it is sometimes useful
to avoid it.  The -n option forces clink to avoid using DNS.  Of
course, this means that you have to provide the destination machine's
IP address instead of a name.

	clink 137.146.210.39

If you are reading data from a dump file, it should be possible
to run clink even if you are not connected to a network.


Verbose mode
------------

The -v flag generates verbose output that looks like this:

clink to host-06.colby.edu (137.146.210.39)
  8 probes at each of 93 sizes (28 to 1500 by 16)
0 localhost
sending datagram to port 33437
received ICMP type= 11  code= 0
ttl=1   137.146.194.17  size=  188 B    rtt= 0.832 ms
sending datagram to port 33438
received ICMP type= 11  code= 0
ttl=1   137.146.194.17  size=  780 B    rtt= 1.742 ms

Notice that the destination port gets incremented each time we
send a message.  That makes it possible to check incoming ICMP
messages and make sure they match up with what we sent.  When
the destination port number reaches the maximum (65535) it wraps
back around to 33437.

ICMP type 11, code 0 indicates a time-exceeded error, which
really means that we reached an intermediate router.

ICMP type 3, code 3 indicates a port unreachable error, which
means we have reached the destination.  That's how clink knows
when to stop.

In verbose mode, clink prints messages about all incoming ICMP
packets, including ones that have nothing to do with clink.


kernel-level timestamps
-----------------------

The -k option only works if you have hacked your kernel to put
timestamps in outgoing UDP packets and incoming ICMP packets.  I
did this to my kernel in order to measure the amount of variance
being added by the user-kernel switch, and the effect of that
variance on clink's estimates.

If you use the -k option with an unmodified kernel, you should
get a message like:

	I don't think the outgoing kernel timestamps are working

In all likelihood, the incoming stamps aren't working either,
but clink didn't get that far.

If you are interested in the hack I used, contact me (info above).



