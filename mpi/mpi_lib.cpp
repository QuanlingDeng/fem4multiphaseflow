#include "mpi_header.h"

using namespace std;

void mpi_Init( int argc, char **argv, int &totalnodes, int &mynode, MPI_Comm comm) 
{
    MPI_Init( &argc, &argv);
    totalnodes = mpi_Size(comm);
    mynode = mpi_Rank(comm);
}

int mpi_Rank(MPI_Comm comm) {

    int mynode;

    MPI_Comm_rank(comm, &mynode);
    return mynode;
}

int mpi_Size(MPI_Comm comm) {

    int totalnodes;

    MPI_Comm_size(comm, &totalnodes);
    return totalnodes;

}

void mpi_Finalize()
{
  MPI_Finalize();
  exit(0);
}

int mpi_Send(void *message, int count, MPI_Datatype datatype, int dest, 
	      int tag, MPI_Comm comm)
{
    return MPI_Send(message, count, datatype, dest, tag, comm);
}

int mpi_Recv(void* message, int count, MPI_Datatype datatype, int source, 
	     int tag, MPI_Comm comm) 
{
    MPI_Status status;
    return MPI_Recv(message, count, datatype, source, tag, comm, &status);
}

int mpi_Isend(void* message, int count, MPI_Datatype datatype, int dest, 
	      int tag, MPI_Request *request, MPI_Comm comm) 
{

    return MPI_Isend(message, count, datatype, dest, tag, comm, request);
}

int mpi_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, 
		 int sendtag, void* recvbuf, int recvcount, MPI_Datatype recvtype,
		 int source, int recvtag, MPI_Status *status, MPI_Comm comm) 
{
  return MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, 
		      recvcount, recvtype, source, recvtag, comm, status);
}

int mpi_Sendrecv_replace(void* buffer, int count, MPI_Datatype datatype, 
			 int dest, int sendtag, int source, int recvtag, 
			 MPI_Status *status, MPI_Comm comm) 
{
  return MPI_Sendrecv_replace(buffer, count, datatype, dest, sendtag, 
			      source, recvtag, comm, status);
}

int mpi_Irecv(void* message, int count, MPI_Datatype datatype, int source, 
	      int tag, MPI_Request *request, MPI_Comm comm) 
{

  return MPI_Irecv( message, count, datatype, source, tag, comm, request);
}

int mpi_Wait(MPI_Request* request, MPI_Status* status)
{
  return MPI_Wait( request, status);
}


int mpi_Barrier(MPI_Comm comm) 
{
    return MPI_Barrier(comm);
}

double mpi_Wtime()
{
  return MPI_Wtime();
}

double mpi_Wtick()
{
  return MPI_Wtick();
}

int mpi_Reduce(void* operand, void* result, int count,
	       MPI_Datatype datatype, MPI_Op operate, int root, MPI_Comm comm) 
{
  return MPI_Reduce(operand, result, count, datatype, operate, root, comm);
}

int mpi_Allreduce(void* operand, void* result, int count, 
		  MPI_Datatype datatype, MPI_Op operate, MPI_Comm comm) 
{
  return MPI_Allreduce(operand, result, count, datatype, operate, comm);
}

int mpi_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
	       void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) 
{
  return mpi_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
}

int mpi_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, 
		  int recvcount, MPI_Datatype recvtype, MPI_Comm comm) 
{
  return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}

int mpi_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		   void *recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm)
{
  return MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
}

int mpi_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
  return MPI_Bcast(buffer, count, datatype, root, comm);
}

int mpi_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, 
		int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) 
{
  return mpi_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
}

int mpi_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, 
		 int recvcount, MPI_Datatype recvtype, MPI_Comm comm) 
{
  return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
