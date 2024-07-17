
#include "output.hpp"
#include "lbm.hpp"
#include "cantera.hpp"
#include "units.hpp"
#include <stdio.h>
#include <iomanip>

void OutputVTK(int &nout, LBM *lbm)
{
	LBM lb = *lbm;
    int Nx = lb.get_Nx();
    int Ny = lb.get_Ny();
    int Nz = lb.get_Nz();
	int dx = lb.get_dx();
	int nSpecies = lb.get_nSpecies();
	std::vector<std::string> speciesName = lb.get_speciesName();
	int		i,j,k;
	char	filename[128];
	FILE	*fp;
	unsigned int array_size;
	unsigned long int offset=0;
	// short  num16; // Int16 2byte
	float  val32; // Float32 4byte

	sprintf(filename,"./field%06d.vtr",nout);
	fp=fopen(filename,"wb");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
	fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",Nx,Ny,Nz);
	fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Nx,Ny,Nz);
	fprintf(fp,"    <PointData>\n");
	fprintf(fp,"    </PointData>\n");
	fprintf(fp,"    <CellData>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Temperature\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Pressure\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);
	for (int a = 0; a < nSpecies ; ++a) 
		{fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Mole Fraction %s\" format=\"appended\" offset=\"%ld\" />\n", speciesName[a].c_str(), offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);}

	for (int a = 0; a < nSpecies ; ++a) 
		{fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Species Velocity %s\" format=\"appended\" offset=\"%ld\" />\n", speciesName[a].c_str(), offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);}
	fprintf(fp,"    </CellData>\n");
	fprintf(fp,"    <Coordinates>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Ny+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nz+1);
	fprintf(fp,"    </Coordinates>\n");
	fprintf(fp,"  </Piece>\n");
	fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"  <AppendedData encoding=\"raw\">");
	fprintf(fp,"_");

    // Density (cell)
	array_size=1*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				#ifdef OUTPUT_SI
					val32=(float)units.si_rho(lb.mixture[i][j][k].rho); fwrite(&val32,sizeof(float),1,fp);
				#else
					val32=(float)lb.mixture[i][j][k].rho; fwrite(&val32,sizeof(float),1,fp);
				#endif
			}
		}
	}
    // Velocity (cell)
	array_size=3*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				#ifdef OUTPUT_SI
					val32=(float)units.si_u(lb.mixture[i][j][k].u); fwrite(&val32,sizeof(float),1,fp);
					val32=(float)units.si_u(lb.mixture[i][j][k].v); fwrite(&val32,sizeof(float),1,fp);
					val32=(float)units.si_u(lb.mixture[i][j][k].w); fwrite(&val32,sizeof(float),1,fp);
				#else
					val32=(float)lb.mixture[i][j][k].u; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.mixture[i][j][k].v; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)lb.mixture[i][j][k].w; fwrite(&val32,sizeof(float),1,fp);
				#endif
			}
		}
	}

    // CellType (cell)
	array_size=1*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				val32=(int)lb.mixture[i][j][k].type; fwrite(&val32,sizeof(int),1,fp);
			}
		}
	}

	// Temperature (cell)
	array_size=1*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				#ifdef OUTPUT_SI
					val32=(float)units.si_temp(lb.mixture[i][j][k].temp); fwrite(&val32,sizeof(float),1,fp);
				#else
					val32=(float)lb.mixture[i][j][k].temp; fwrite(&val32,sizeof(float),1,fp);
					// val32=(float)lb.mixture[i][j][k].internalEnergy; fwrite(&val32,sizeof(float),1,fp);
				#endif
			}
		}
	}

	// Pressure (cell)
	array_size=1*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				#ifdef OUTPUT_SI
					val32=(float)units.si_p(lb.mixture[i][j][k].p); fwrite(&val32,sizeof(float),1,fp);
				#else
					val32=(float)lb.mixture[i][j][k].p; fwrite(&val32,sizeof(float),1,fp);
				#endif
			}
		}
	}

	#ifdef MULTICOMP
	for(int a = 0; a < nSpecies; ++a)
	{
		array_size=1*4*(Nx)*(Ny)*(Nz);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=0;k<Nz;k++){
			for(j=0;j<Ny;j++){
				for(i=0;i<Nx;i++)
				{
					// val32=(float)( units.si_rho(lb.species[a][i][j][k].rho) ) ; fwrite(&val32,sizeof(float),1,fp);
					val32=(float)( lb.species[a][i][j][k].X ) ; fwrite(&val32,sizeof(float),1,fp);
				}
			}
		}
	}

	for(int a = 0; a < nSpecies; ++a)
	{
		array_size=1*4*(Nx)*(Ny)*(Nz);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=0;k<Nz;k++){
			for(j=0;j<Ny;j++){
				for(i=0;i<Nx;i++)
				{
					#ifdef OUTPUT_SI
						val32=(float) units.si_u(lb.species[a][i][j][k].Vdiff_x) ; fwrite(&val32,sizeof(float),1,fp);
					#else
						val32=(float) lb.species[a][i][j][k].u-lb.mixture[i][j][k].u ; fwrite(&val32,sizeof(float),1,fp);
					#endif					
				}
			}
		}
	}

	
	#endif

	// Coordinates (vertices)
	array_size=1*4*(Nx+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(i=0;i<Nx+1;i++){ val32=(float)(i*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Ny+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(j=0;j<Ny+1;j++){ val32=(float)(j*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Nz+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz+1;k++){ val32=(float)(k*dx); fwrite(&val32,sizeof(float),1,fp); }

	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

void OutputKeEns(int &step, LBM *lbm)
{	
	LBM lb = *lbm;
	int Nx = lb.get_Nx();
    int Ny = lb.get_Ny();
    int Nz = lb.get_Nz();
	double e_kinetic = 0;
	double enstro = 0;
	double uu = 0 ;
	double vort = 0;
	double toten = 0.0; // Total Energy

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				if (lb.mixture[i][j][k].type == TYPE_F)
				{ 
					uu = lb.mixture[i][j][k].u*lb.mixture[i][j][k].u + lb.mixture[i][j][k].v*lb.mixture[i][j][k].v + lb.mixture[i][j][k].w*lb.mixture[i][j][k].w;
					e_kinetic = e_kinetic + lb.mixture[i][j][k].rho * uu;
					
					int i_a = i + 1;
					int i_b = i - 1;
					int j_a = j + 1;
					int j_b = j - 1;
					int k_a = k + 1;
					int k_b = k - 1;

					if (i_b < 1) i_b = Nx-2;
					else if(i_a > Nx-2) i_a = 1;

					if (j_b < 1) j_b = Ny-2;
					else if(j_a > Ny-2) j_a = 1;

					if (k_b < 1) k_b = Nz-2;
					else if(k_a > Nz-2) k_a = 1;
					

					double uy = (lb.mixture[i][j_a][k].u - lb.mixture[i][j_b][k].u) / 2;
					double uz = (lb.mixture[i][j][k_a].u - lb.mixture[i][j][k_b].u) / 2;
					
					double vx = (lb.mixture[i_a][j][k].v - lb.mixture[i_b][j][k].v) / 2;
					double vz = (lb.mixture[i][j][k_a].v - lb.mixture[i][j][k_b].v) / 2;

					double wx = (lb.mixture[i_a][j][k].w - lb.mixture[i_b][j][k].w) / 2;
					double wy = (lb.mixture[i][j_a][k].w - lb.mixture[i][j_b][k].w) / 2;

					vort = (wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy);
					enstro = enstro + lb.mixture[i][j][k].rho * vort;

					#ifndef ISOTHERM
					toten = toten + lb.mixture[i][j][k].rhoe;
					#endif
				}
			}
		}
	}
	
	std::cout << step << ",";
	std::cout << e_kinetic << ",";
	std::cout << enstro << ",";
	std::cout << std::setprecision(15) << std::fixed << toten << std::endl;
	

}

void calcError(int &t,LBM &lb)
{
	int Nx = lb.get_Nx();
    int Ny = lb.get_Ny();
    int Nz = lb.get_Nz();
	int NX = lb.get_NX();
    int NY = lb.get_NY();
	int NU = lb.get_nu();
	
	double kx = 2.0*M_PI/NX;
	double ky = 2.0*M_PI/NY;

	double ua = 0.0;
	double va = 0.0;

	double sumue2 = 0.0;
	double sumve2 = 0.0;

	double sumua2 = 0.0;
	double sumva2 = 0.0;

	double td = 1.0/(NU*(kx*kx+ky*ky));

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				if (lb.mixture[i][j][k].type == TYPE_F)
				{
					double X = i ;
					double Y = j ;
					ua =-0.04 * cos(kx*X) * sin(ky*Y) *exp(-1.0*t/td);
					va = 0.04 * sin(kx*X) * cos(ky*Y) *exp(-1.0*t/td);

					sumue2 += (lb.mixture[i][j][k].u/lb.mixture[i][j][k].rho - ua) * (lb.mixture[i][j][k].u/lb.mixture[i][j][k].rho - ua);
					sumve2 += (lb.mixture[i][j][k].v/lb.mixture[i][j][k].rho - va) * (lb.mixture[i][j][k].v/lb.mixture[i][j][k].rho - va);
					
					sumua2 += ua * ua;
					sumva2 += va * va;
				}
			}
		}
	}

	std::cout << "u error		= " << sqrt(sumue2/sumua2) << std::endl;
	std::cout << "v error		= " << sqrt(sumve2/sumva2) << std::endl; 

}

void printLogo()
{
}