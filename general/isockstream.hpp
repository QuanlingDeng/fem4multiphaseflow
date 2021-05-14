/**************************************************************************/
/* File:   isockstream.hpp                                                */
/* Author: Tony                                                           */
/* Date:   06/29/01                                                       */
/**************************************************************************/
#include <string.h>
#include <stdlib.h>
#include <strstream>
using namespace std;
#include <errno.h>
#include <signal.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>


/** Data type for input socket stream class. The class is used as server
    to receive data from a client on specified port number. The user gets
    data from the stream as from any othr input stream.                  **/

class isockstream
{
   private:
      int portnum, portID, socketID;
      char *Buf;

      int establish();
      int read_data(int socketid, char *buf, int size); 

   public:

      /** The constructor takes as input the portnumber port on which
          it establishes a server.                                       **/
      isockstream(int port);

      //  Start waiting for data and return it in an input stream.
      void receive(istrstream **in);

      /** Virtual destructor. If the data hasn't been sent it sends it.  **/
      ~isockstream();
};

//=============================================================================

isockstream::isockstream(int port){
   portnum = port;

   if ( (portID = establish()) < 0)	
      cout << "Server couldn't be established on port " 
           << portnum << endl;
   Buf = NULL;
}

//=============================================================================

int isockstream::establish(){
   char   myname[129];
   int    port;
   struct sockaddr_in sa;
   struct hostent *hp;

   memset(&sa, 0, sizeof(struct sockaddr_in));
   gethostname(myname, 128);
   hp= gethostbyname(myname);

   if (hp == NULL)
     return(-1);

   sa.sin_family= hp->h_addrtype;
   sa.sin_port= htons(portnum);

   if ((port = socket(AF_INET, SOCK_STREAM, 0)) < 0)
     return(-1);

   int on=1;
   setsockopt(port, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));

   if (bind(port,(const sockaddr*)&sa,sizeof(struct sockaddr_in)) < 0) {
     close(port);
     return(-1);
   }
 
   listen(port, 4);
   return(port);
}

//=============================================================================

int isockstream::read_data(int s, char *buf, int n){ 
  int bcount;                      // counts bytes read
  int br;                          // bytes read this pass

  bcount= 0;
  br= 0;
  while (bcount < n) {             // loop until full buffer
    if ((br= read(s,buf,n-bcount)) > 0) {
      bcount += br;                // increment byte counter
      buf += br;                   // move buffer ptr for next read
    }
    else if (br < 0)               // signal an error to the caller
      return(-1);
  }
  return(bcount);
}

//=============================================================================

void isockstream::receive(istrstream **in){
   int size;
   char length[32];

   if (portID >= 0)
     if ((socketID = accept(portID, NULL, NULL)) < 0)
        if (socketID != EINTR)
	  cout << "Server failed to accept connection." << endl;

   read(socketID, length, 32);
   size = atoi(length);
    
   if (Buf != NULL)
     delete [] Buf;
   Buf = new char[size+1];
   if (size != read_data(socketID, Buf, size))
     cout << "Not all the data has been read" << endl;
   else
#ifdef DEBUG
     cout << "Reading " << size << " bytes is successful" << endl;
#endif
   Buf[size] = '\0';
   
   close(socketID); 
   if ((*in) != NULL)
     delete (*in);
   (*in) = new istrstream(Buf);
}

//=============================================================================

isockstream::~isockstream(){
   if (Buf != NULL)
      delete [] Buf;
   close(portID);
}

//=============================================================================
