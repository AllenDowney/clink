#pragma once
void clink_init (char *host, int tos);
int send_probe (int size, int ttl, int timeout, Datum *datum);
