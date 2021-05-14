#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main()
{
  double dt = 5;
  for(int i=0; i<=7200; i+=3600)
    {
      double time = i;
      char gle[256];
      sprintf(gle, "saturation_%d.gle", (int) time);
      ofstream file(gle);

      char satname[256];
      sprintf(satname, "saturation_%d.dat", (int) time);
      ifstream inputfile;
      inputfile.open(satname);

      file << "include \"color.gle\"" << endl;
      file << "include \"shape.gle\"" << endl;
      file << endl;

      file << "sub red z" << endl;
      file << " local r = 0" << endl; 
      file << "if (z <= 0.25)                 then r = 0.5 - z*4" << endl;
      file <<  "if (z > 0.25) and (z <= 0.50)  then r = 0" << endl;
      file << "if (z > 0.50) and (z <= 0.75)  then r = 0" << endl;
      file << "if (z > 0.75)                  then r = 0" << endl;
      file << "return r" << endl;
      file << "end sub" << endl;
      file << endl;

      file << "sub green z" << endl;
      file << "local g = 0" << endl;
      file << "if (z <= 0.25)                 then g = z*4" << endl;
      file << "if (z > 0.25) and (z <= 0.50)  then g = 1 - (z-0.25)*4" << endl;
      file << "if (z > 0.50) and (z <= 0.75)  then g = 0 " << endl;
      file << "if (z > 0.75)                  then g = 0" << endl;
      file << "return g" << endl;
      file << "end sub" << endl;
      file << endl;

      file << "sub blue z" << endl;
      file << "local b = 0" << endl;
      file << "if (z <= 0.25)                 then b = z*4" << endl;
      file << "if (z >  0.25) and (z <= 0.50) then b = 1" << endl;
      file << "if (z >  0.50) and (z <= 0.75) then b = 1" << endl;
      file << "if (z >  0.75)                 then b = 1 - (z-0.75)*2" << endl;
      file << "return b" << endl;
      file << "end sub" << endl;
      file << endl;

      file << "size 5.0 1.0" << endl;

      file << "set lwidth 0.001" << endl;
      file << "set cap square" << endl;
      file << "set lstyle 0" << endl;
      file << endl;

      int ncv;
      inputfile >> ncv;
      for(int k=0; k<ncv; k++)
	{
	  int numvert;
	  double satval;
	  inputfile >> numvert >> satval;
	  double vertices[numvert];
	  for(int j=0; j<numvert; j++){ inputfile >> vertices[j]; }  

	  
	  file << "r = red(" << satval << ")" << endl;
	  file << "g = green(" << satval << ")" << endl;
	  file << "b = blue(" << satval << ")" << endl;

	  file << "set color rgb(r, g, b)" << endl;	  	  
	  file << "begin path stroke fill rgb(r, g, b)" << endl;
	  file << "amove " << vertices[0] << " " << vertices[1] << endl; 
	  for(int j=2; j<numvert; j+=2)
	    {
	      file << "aline " << vertices[j] << " " << vertices[j+1] << endl;
	    }  
	  file << "closepath" << endl;
	  file << "end path" << endl;
	  file << endl;

	}





      //file << "amove xg(xgmax)+0.3 yg(ygmin)" << endl;
      //file << "color_range_vertical zmin 0 zmax 1 zstep 0.1 format \"fix 1\" pixels 500 palette palette_custom" << endl;

      file.close();
      inputfile.close();
    }

}

