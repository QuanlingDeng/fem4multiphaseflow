/**************************************************************************/
/* File:   osockstream.cpp                                                */
/* Author: Tony, Veselin                                                  */
/* Date:   02/27/01                                                       */
/**************************************************************************/

#include "osockstream.hpp"
#include <iostream>
#include <stdio.h>
#include <netdb.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <netinet/in.h>


osockstream::osockstream(int port, const char *hostname)
   : ostrstream()
{
   int n;

   sent = 0;
   portnum = port;
   n = strlen(hostname);
   if (n >= 128)
   {
      strncpy(machine, hostname, 127);
      machine[127] = '\0';
   }
   else
      strcpy(machine, hostname);
}

int osockstream::call_socket()
{
   struct sockaddr_in  sa;
   struct hostent     *hp;
   int s;

   if ( (hp= gethostbyname(machine)) == NULL )
      return(-1);

   memset(&sa,0,sizeof(sa));
   memcpy((char *)&sa.sin_addr,hp->h_addr,hp->h_length); /* set address */
   sa.sin_family= hp->h_addrtype;
   sa.sin_port= htons((u_short)portnum);

   if ((s= socket(hp->h_addrtype,SOCK_STREAM,0)) < 0)   /* get socket */
      return(-2);
   if (connect(s,(const sockaddr*)&sa,sizeof sa) < 0) { /* connect */
      close(s);
      return(-3);
   }
   return(s);
}

int osockstream::send_through_socket(int s, char *Buf, int size)
{
   int bcount = 0; /* counts bytes read */
   int br     = 0; /* bytes read this pass */
   char *b    = Buf;

   char length[32];
   sprintf(length,"%d", size);

   ::write(s, length, 32);

   while (bcount < size)
      if ((br= ::write(s,b,size-bcount)) > 0) {
         bcount += br;
         b += br;
      }
      else 
         if (br < 0)
            return (-1);
   return 1;
}

int osockstream::send()
{
   int sock, size, ierr = 0;
   char *Buf;

   if (sent)
      return (-1); // Data already sent
   sent = 1;
   switch ( sock = call_socket() )
   {
      case -1:
         cerr << "Unknown host: " << machine << endl;
         ierr = -2; break;
      case -2:
         cerr << "Unable to establish socket on the local machine" << endl;
         ierr = -3; break;
      case -3:
         cerr << "Unable to connect to port " << portnum << " on "
              << machine << endl;
         ierr = -4; break;
      default:
         Buf = str();
         size = pcount();
         switch ( send_through_socket(sock, Buf, size) )
         {
            case -1:
               cerr << "Error sending data to port " << portnum << " on "
                    << machine << endl;
               ierr = -5; break;
            default:
               break;
         }
         close(sock);
         delete [] Buf;
         break;
   }
   return ierr;
}

osockstream::~osockstream()
{
   if (!sent)
      send();
}

