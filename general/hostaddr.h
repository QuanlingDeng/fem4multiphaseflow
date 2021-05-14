#include <unistd.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <string.h>

int getfullhostname(char *name, int maxlen)
{
   struct hostent *he;

   if (gethostname (name, maxlen))
      return (-1);
   he = gethostbyname (name);
   if (!he)
      return (-2);
   strncpy(name, he->h_name, maxlen);
   name[maxlen-1] = '\0';
   return 0;
}

int gethostip(char *ip, int maxlen)
{
   struct hostent *he;
   unsigned char *addr;

   if (gethostname (ip, maxlen))
      return (-1);
   he = gethostbyname (ip);
   if (!he)
      return (-2);
   addr = (unsigned char *) ( he->h_addr_list[0] );
   switch (he->h_length)
   {
      case 4:  snprintf(ip, maxlen, "%u.%u.%u.%u", addr[0], addr[1],
                                                   addr[2], addr[3]);
               break;
      default:
               return (-3);
   }
   return 0;
}
