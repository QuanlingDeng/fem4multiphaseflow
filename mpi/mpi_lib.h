#include <mpi.h>

//start the main thread, call comm_rank and comm_size
//init performs necissary initializations
//size = number of processors, P
//rank = process identification \in (0, P-1) 
void mpi_Init(int argc, char **argv, int &totalnodes, int &mynode, 
	      MPI_Comm comm = MPI_COMM_WORLD);

int mpi_Rank(MPI_Comm comm = MPI_COMM_WORLD);

int mpi_Size(MPI_Comm comm = MPI_COMM_WORLD);

//ends the mpi implimentation. 
//note: the program runs entirely in parrelel, 
//so this does not stop "end" the parallel implementation
void mpi_Finalize();

//Performs a blocking send
int mpi_Send(void *message, int count, MPI_Datatype datatype,
	     int dest, int tag, MPI_Comm comm = MPI_COMM_WORLD);

//Performs a blocking recv 
int mpi_Recv(void* message, int count, MPI_Datatype datatype, 
	     int source, int tag, MPI_Comm comm = MPI_COMM_WORLD);

/*
MPI_Send and MPI_Recv are blocking communications, which means 
that they will not return until it is safe to modify or use the
contents of the send/recv buffer respectively. 

note: send and recieve are coupled, from one proccesor you send, 
and from another you recieve.
*/ 

//Begins a nonblocking send 
int mpi_Isend(void* message, int count, MPI_Datatype datatype, 
	      int dest, int tag, MPI_Request *request, MPI_Comm comm = MPI_COMM_WORLD);
 
//Begins a nonblocking receive 
int mpi_Irecv(void* message, int count, MPI_Datatype datatype, 
	      int source, int tag, MPI_Request *request, MPI_Comm comm = MPI_COMM_WORLD);

//wait till you have completed the isend and irecv command
int mpi_Wait(MPI_Request* request, MPI_Status* status);


/*
MPI_Send and MPI_Recv are nonblocking communications, which means 
that they allow a process to post that it wants to send to or recieve 
from a process, and then later allows it to call a function MPI_Wait 
to complete the send/recieving. These functions are useful in that 
they allow the programer to appropriately stagger computation and 
communication to minimize the total waiting time resulting from 
communication.

note: the I stands for immediate!
*/ 

//Sends and receives a message 
int mpi_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
		 int dest, int sendtag, void* recvbuf, int recvcount, 
		 MPI_Datatype recvtype, int source, int recvtag, MPI_Status *status,
		 MPI_Comm comm = MPI_COMM_WORLD);

//Sends and receives using a single buffer 
int mpi_Sendrecv_replace(void* buffer, int count, MPI_Datatype datatype, int dest, 
			 int sendtag, int source, int recvtag, MPI_Status *status,
			 MPI_Comm comm = MPI_COMM_WORLD);

//Blocks until all processes in the communicator have reached this routine. 
int mpi_Barrier(MPI_Comm comm = MPI_COMM_WORLD);


//Returns an elapsed time on the calling processor 
double mpi_Wtime();

//Returns the resolution of MPI_Wtime 
double mpi_Wtick();

//Reduces values on all processes to a single value
int mpi_Reduce(void* operand, void* result, int count, MPI_Datatype datatype, 
	       MPI_Op operate, int root, MPI_Comm comm  = MPI_COMM_WORLD);

//Combines values from all processes and distributes the result back to all processes
int mpi_Allreduce(void* operand, void* result, int count, MPI_Datatype datatype, 
		  MPI_Op operate, MPI_Comm comm = MPI_COMM_WORLD);

//Gathers together values from a group of processes 
int mpi_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, 
	       int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm = MPI_COMM_WORLD); 

// Gathers data from all tasks and distribute the combined data to all tasks 
int mpi_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, 
		  int recvcount, MPI_Datatype recvtype, MPI_Comm comm = MPI_COMM_WORLD);

//Gathers data from all tasks and distribute the combined data to all tasks 
int mpi_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, 
		   int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm = MPI_COMM_WORLD);

//Broadcasts a message from the process with rank "root" to all other processes of the communicator 
int mpi_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm = MPI_COMM_WORLD);

// Sends data from one process to all other processes in a communicator 
int mpi_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, 
		int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm = MPI_COMM_WORLD);

//Sends data from all to all processes 
int mpi_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
		 int recvcount, MPI_Datatype recvtype, MPI_Comm comm = MPI_COMM_WORLD);
