//A
//  Created by G Santambrogio on 15/8/14.
//  Copyright (c) 2014 ambrogio. All rights reserved.
//

// compile with: gcc -o MolsGen MolGen.c -O3 -lm -lgsl -lgslcblas -lconfig 

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <libconfig.h>


static const double fwhm2sigma = 2.35482; // 2 * sqrt (2*ln(2) )= 2.35482  FWHM=2.35482*sigma

int main (int argc, char *argv[]){
    int i;
    //double t;
    //y = {v_x, v_y, v_z, x, y, z}
    //f = {a_x, a_y, a_z, v_x, v_y, v_z}
    double y[6];
    
    if(argc < 2){
        fprintf(stderr,  " The program has different usages.\n Case 1:\n It calculates a distribution of molecules\n generated at position <z nozzle>, excited at a position <z laser>\n by a laser with spotsize of <Radius Laser>,\n and arriving at the entrance of the chip <z chip>.\n Only molecules within +/- halfHeight in the y direction are taken into account.\n All numbers in micron.\n Usage: MolsGen 1 <number of molecules> <z nozzle> <z laser> <Radius Laser> <z Chip entrance> \n \n Case 2:\n It returns an homogeneous and isotropic distribution in x, y, z, vx, vy, and vz centered at {0,0,0,0,0,0} (all units in um and us)\n Usage: MolsGen 2 <number of molecules> <x width> <y width> <z width> <vx width> <vy width> <vz width> <vz>\n\n Case 3:\n It calculates a distribution of molecules\n generated at position <z nozzle>, excited at a position <z laser>\n by a laser with spotsize of <Radius Laser>.\n The laser fires at the right time to capture the middle of the molecular cloud.\n\n");
        return 1;
    }
    int gmod = atoi(argv[1]);
    
    if(argc != 7 && gmod==1){
        fprintf(stderr, " For Case 1:\n Usage: MolsGen3D 1 <number of molecules> <z nozzle> <z laser> <Radius Laser> <z Chip entrance> \n");
        return 1;
    }
    if(argc != 10 && gmod==2){
        fprintf(stderr, " For Case 2:\n Usage: MolsGen 2 <number of molecules> <x width> <y width> <z width> <vx width> <vy width> <vz width> <vz>\n");
        return 1;
    }
    if(argc != 6 && gmod==3){
        fprintf(stderr, " For Case 3:\n Usage: MolsGen 3 <number of molecules> <z nozzle> <z laser> <Radius Laser>  \n");
        return 1;
    }
    
    config_t cfg;

    config_init(&cfg);
    if(! config_read_file(&cfg, "config.cfg"))
      {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
              config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    }

    switch(gmod){
        case 1:
        {
            int imax = atoi(argv[2]);
            double znozzle = atof(argv[3]);
            double zlas = atof(argv[4]);
            double lasradius = atof(argv[5]);
            double zChipbegin = atof(argv[6]);
            double SourceRadius=config_setting_get_float(config_lookup(&cfg, "molSource.ValveDia"));
	    double velocity=config_setting_get_float(config_lookup(&cfg, "molSource.vz"));
	    double velFWHM=config_setting_get_float(config_lookup(&cfg, "molSource.vz_width"));
	    double valveopentime=config_setting_get_float(config_lookup(&cfg, "molSource.ValvePulseDuration"));
	    double halfHeight=0.5*config_setting_get_float(config_lookup(&cfg, "microchip.entranceSlitHeight"));
	    double halfWidth=0.5*config_setting_get_float(config_lookup(&cfg, "microchip.arrayWidth"));	    	    
            gsl_rng * r=gsl_rng_alloc (gsl_rng_mt19937);
            gsl_rng_set(r, 0);
            
            //  goodcout<<"#time"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"vx"<<"\t"<<"vy"<<endl;
	    printf("name \t #time \t x \t y \t z \t vx \t vy \t vz \n");

            i=0;
            while(i<imax){
                double rho = SourceRadius*sqrt(gsl_rng_uniform(r));
                double phi = 2.0*M_PI*gsl_rng_uniform(r);
                y[3] = rho*cos(phi);
                y[4] = rho*sin(phi);
                
                //cout<<i<<"\t"<<y[3]<<"\t"<<y[4]<<endl;
                //first approximation: angle = sin = tan: @ 3.5 degrees the error is 7e-5
                y[2] = velocity+gsl_ran_gaussian(r, velFWHM/fwhm2sigma); // gsl_ran_gaussian(double x, double sigma)
                double angleY = gsl_ran_gaussian (r, 0.05); //Rotation around Y axis, i.e. in the xz plane
                double angleX = gsl_ran_gaussian (r, 0.05);// Rotation around X axis, i.e. in the yz plane
                //now we have some 40 cm of flight and a distribution with a FWHM of about 4 cm, so sigma is about 2 cm
                //So sigma of the angle is atan(2/40)
                y[1] = y[2]*angleX;
                y[0] = y[2]*angleY;
                
                y[5] = znozzle+velocity*valveopentime*(gsl_rng_uniform(r)-0.5);
                
                double tlas = (zlas - znozzle)/velocity; // this can be changed      cout<<tlas<<endl;
                
                // totcout<<t<<"\t"<<y[2]<<"\t"<<y[3]<<"\t"<<y[0]<<"\t"<<y[1]<<endl;
                // The molecules arrive at the center of the laser spot at a time (zlas-y[5])/y[2]
                // so the following condition is on time for the z direction, whether the molecules see the laser,
                // and on space for the x,y directions, assuming that every molecule within the laser spot is excited.
                if( (((zlas-lasradius-y[5])/y[2])<tlas)
                   && (tlas<((zlas+lasradius-y[5])/y[2]))
                   && (fabs(y[3]+(zlas-y[5])*angleY) < lasradius)
                   && (fabs(y[4]+(zlas-y[5])*angleX) < lasradius)){
                    //    goodcout<<(xlas-y[2])/y[0]<<"\t"<<y[2]<<"\t"<<y[3]+(xlas-y[2])*angle<<"\t"<<y[0]<<"\t"<<y[1]<<endl;
                    double xPosAtChipBegin = y[3]+(zChipbegin - y[5]) * angleY;
                    double yPosAtChipBegin = y[4]+(zChipbegin - y[5]) * angleX;
                    if( (fabs(xPosAtChipBegin) < halfWidth) && ( fabs(yPosAtChipBegin) < halfHeight) ) {
                        //goodcout<<(xslit-y[2])/y[0]<<"\t"<<xslit<<"\t"<<yPosAtSlit<<"\t"<<y[0]<<"\t"<<y[1]<<endl;
                        printf("%d %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", i, (zChipbegin - y[5])/y[2], xPosAtChipBegin, yPosAtChipBegin, zChipbegin, y[0], y[1], y[2]);
                        i++;
                    }
                    
                    //        cout<<(xslit-y[2])/y[0]<<"\t"<<xslit<<"\t"<<yPosAtSlit<<"\t"<<y[0]<<"\t"<<y[1]<<endl;
                }
                
            }
        }
            break;
            
        case 2:
        {
            int imax = atoi(argv[2]);
            double xwidth = atof(argv[3]);
            double ywidth = atof(argv[4]);
            double zwidth = atof(argv[5]);
            double vxwidth = atof(argv[6]);
            double vywidth = atof(argv[7]);
            double vzwidth = atof(argv[8]);
	    double vz = atof(argv[9]);
            
            gsl_rng * r=gsl_rng_alloc (gsl_rng_mt19937);
            gsl_rng_set(r, 0);
            printf("name \t #time \t x \t y \t z \t vx \t vy \t vz \n");
            for(i=0;i<imax;i++){
                y[3]=xwidth*(gsl_rng_uniform(r)-0.5);
                y[4]=ywidth*(gsl_rng_uniform(r)-0.5);
                y[5]=zwidth*(gsl_rng_uniform(r)-0.5);
                y[0]=vxwidth*(gsl_rng_uniform(r)-0.5);
                y[1]=vywidth*(gsl_rng_uniform(r)-0.5);
                y[2]=vzwidth*(gsl_rng_uniform(r)-0.5)+vz;
                printf("%d %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", i, 0.0, y[3], y[4], y[5], y[0], y[1], y[2]);
            }
        }
            break;
        case 3:
        {
            int imax = atoi(argv[2]);
            double znozzle = atof(argv[3]);
            double zlas = atof(argv[4]);
            double lasradius = atof(argv[5]);
            double SourceRadius=0.5*config_setting_get_float(config_lookup(&cfg, "molSource.ValveDia"));
	    double velocity=config_setting_get_float(config_lookup(&cfg, "molSource.vz"));
	    double velxFWHM=config_setting_get_float(config_lookup(&cfg, "molSource.vx_width"));
	    double velyFWHM=config_setting_get_float(config_lookup(&cfg, "molSource.vy_width"));
	    double velzFWHM=config_setting_get_float(config_lookup(&cfg, "molSource.vz_width"));	    
	    double valveopentime=config_setting_get_float(config_lookup(&cfg, "molSource.ValvePulseDuration"));	    
            gsl_rng * r=gsl_rng_alloc (gsl_rng_mt19937);
            gsl_rng_set(r, 0);
	    printf("name \t #time \t x \t y \t z \t vx \t vy \t vz \n");
            
            i=0;
            while(i<imax){
                double rho = SourceRadius*sqrt(gsl_rng_uniform(r));
                double phi = 2.0*M_PI*gsl_rng_uniform(r);
                y[3] = rho*cos(phi);
                y[4] = rho*sin(phi);

		y[0] = gsl_ran_gaussian(r, velxFWHM/fwhm2sigma);
		y[1] = gsl_ran_gaussian(r, velyFWHM/fwhm2sigma); 		
                y[2] = velocity+gsl_ran_gaussian(r, velzFWHM/fwhm2sigma); // gsl_ran_gaussian(double x, double sigma)
		//double angleY = gsl_ran_gaussian (r, 0.05); //Rotation around Y axis, i.e. in the xz plane
                //double angleX = gsl_ran_gaussian (r, 0.05);// Rotation around X axis, i.e. in the yz plane
                //now we have some 40 cm of flight and a distribution with a FWHM of about 4 cm, so sigma is about 2 cm
                //So sigma of the angle is atan(2/40)
		//NB This is true for Kripton. Think it again for Neon!!!! #############
                //y[1] = y[2]*angleX;
                //y[0] = y[2]*angleY;
                
                y[5] = znozzle+velocity*valveopentime*(gsl_rng_uniform(r)-0.5);
                
                double tlas = (zlas - znozzle)/velocity; // this can be changed 
                

                // The molecules arrive at the center of the laser spot at a time (zlas-y[5])/y[2]
                // so the following condition is on time for the z direction, whether the molecules see the laser,
                // and on space for the x,y directions, assuming that every molecule within the laser spot is excited.
		double xPosAtLas = y[3]+tlas*y[0];
		double yPosAtLas = y[4]+tlas*y[1];
                if( (((zlas-lasradius-y[5])/y[2])<tlas)
                   && (tlas<((zlas+lasradius-y[5])/y[2]))
		    && ((xPosAtLas*xPosAtLas+yPosAtLas*yPosAtLas) < (lasradius*lasradius))){

                        printf("%d %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", i, (zlas - y[5])/y[2], xPosAtLas, yPosAtLas, zlas, y[0], y[1], y[2]);
                        i++;
		}
                
            }
        }
            break;

    }
    
    
    return 0;
}
