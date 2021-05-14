/**************************************************************************/
/* File:   osockstream.hpp                                                */
/* Author: Tony, Veselin                                                  */
/* Date:   02/27/01                                                       */
/**************************************************************************/

#include <strstream>
using namespace std;

/** Data type for output socket stream class. The class is used as client
    to send data to a server on a specified port number. One object of the 
    class can be used for one time send of data to the server. The user
    writes in the stream, as in any other output stream and when the data
    is ready to be send function send() has to be executed. Otherwise (if
    not executed) the destructor will send the data.                     **/
class osockstream : public ostrstream
{
   private:
      int portnum, sent;
      char machine[128];

      // return socket number;
      // (-1) for unknown machine;
      // (-2) for socket fail;
      // (-3) for connect fail.
      int call_socket();
      int send_through_socket(int s, char *Buf, int size);

   public:

      /** The constructor takes as input the name of the server and the
          port number through which the communication will take place.   **/
      osockstream(int port, const char *hostname);

      /** Send the current in the stream data to the server specified by
          name "hostname" (in the constructor) on port number "port".    
          Return -1 if data has already been sent or 0 for success.      **/
      int send();

      /** Virtual destructor. If the data hasn't been sent it sends it.  **/
      ~osockstream();
};
