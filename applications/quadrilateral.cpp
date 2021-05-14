#include "../fem/fem_header.h"
#include "../mesh/mesh_header.h"


#define Nx 3
#define Ny 3

double k(double x, double y)
{
    return 1;
}

double exactu(double x, double y) { return -x*x*x/6.0 + x*x/4.0 -1.0/24.0; }

double f(double x, double y)
{
    return ( x-0.5 );
}


int main()
{
    enum Type { POINT, SEGMENT, TRIANGLE, QUADRILATERAL, TETRAHEDRAL };
    
    Array<int> n(2);
    n[0] = Nx;
    n[1] = Ny;
    
    Array<double> L(4);
    L[0] = 0.0;
    L[1] = 1.0;
    L[2] = 0.0;
    L[3] = 1.0;
    
    Mesh *rec_mesh = new TMesh<2>(n, L, Element::QUADRILATERAL);
    
    int TotalNV = (Nx+1)*(Ny+1);
    
    Array<double *> linellipticdata(1);
    Array<double *> linforcedata(1);
    
    linellipticdata[0] = new double[TotalNV];
    linforcedata[0] = new double[TotalNV];
    
    double dx = 1.0/Nx;
    double dy = 1.0/Ny;
    
    for(int i=0; i<=Ny; i++)
    {
        for(int j=0; j<=Nx; j++)
        {
                linellipticdata[0][i*(Nx+1)+j] = k(j*dx, i*dy );
                linforcedata[0][i*(Nx+1)+j] = f(j*dx, i*dy);
        }
    }
    
    
    
    Data *data;
    data = new Data();
    data->SetEllipticData(linellipticdata);
    data->SetForceData(linforcedata);
    
    
    FEM *fun = new BilinearFEM2D(rec_mesh, data);
    fun->Assemble();

    
    Vector sol();
    fun->Solve(sol, 10000, 1);
    
    
    
    ofstream femsol("femsol.txt");
    
    int index = 0;
    for(int j=0; j<=Ny; j++)
    {
        double coord = 0;
        for(int i=0; i<=Nx; i++)
        {
            femsol << exactu(i*dx, j*dy) << " " << sol(index) <<  endl;
            index++;
            
        }
    }
    
    femsol.close();
    
    delete data;
    delete []linellipticdata[0];
    delete []linforcedata[0];
    delete fun;


}
