/*
 * Elena Sellentin
 * University of Leiden
 * Sterrewacht Leiden
 * 
 * Gravitational wave filtering exercise
 * Example program accompanying my lecture
 * "Modern Astrostatistics"
 * 
 * Compilation:
 * g++ -o GravWaves GWFiltering.cpp -lgsl -lgslcblas
 * 
 * Execution:
 * ./GravWaves
 */
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;


double f(double x, double Msolar, double tc)  //functions are named as on the exercise sheet
{
    return 16.6*pow(Msolar,-5./8)*pow(tc-x,-3./8.);
}


double phi(double x, double Msolar, double tc)
{
    return 2*3.1418*16.6*pow(Msolar,-5./8) *(-3./8.) * pow(tc-x,5./8.);
}


double h(double x, double Amplitude, double Msolar, double tc, double phi_c)  //the strain
{
    if(Amplitude > 1.0)
    {
        /*
          the >> and << are like unix pipes. All it says is that the message is then printed to the terminal,
           if you make the mistake of entering an Amplitude greater than 1.0
         */
        cout << "The amplitude should lie between 0 and 1, you entered " << Amplitude << endl;  
    }
    double r = fabs(Amplitude)*pow(f(x,Msolar,tc) ,2./3.) * cos( phi(x,Msolar,tc) + phi_c )  ;
    
    if(x > tc)
    {
        r=0;
    }
    //This is fixing a house-made problem because I am not using full General-Relativity templates, but Newtonian ones instead. The Newtonian ones make 
   // and incorrect spike at the end, which I here cut away
    if( fabs(r) > 0.5 )
    {
        r = 0.0; //Disallowing chi-by-eye due to our lacking ringdown
    }
    return r;
}



double Filter_output(vector<double> datastream, vector<double> datatime, double standard_dev, double Amplitude, double Msolar, double tc, double phi_c )
{
    double SN = 0;
    
    double norm = 0;
    
    int d = datastream.size();
    
    for(int i = 0; i < d; i++)
    {
        double signal = h(datatime[i], Amplitude,Msolar,tc,phi_c);
        SN += pow( datastream[i]*signal,2);
        norm += pow(signal/standard_dev,2);
    }
    
    if(standard_dev > 0.0)
    {
        SN = (SN/norm);
    }
    
    return SN;
}



//C++ programs start at the main function
int main()
{
 
    /*
     * Generating 3 signals;
     */
    int points = 2*1500;
    vector<double> time(points);   //C's version of np.zeros 
    vector<double> GW1(points);
    vector<double> GW2(points);
    vector<double> GW3(points);
    
    double tmin = -1500;  //arbitrary zero point for my time
    for(int i = 0; i < points;i++)
    {
        time[i] = tmin;       //These are the correct values for the GWs; if your filtering works out, the filters should spike for these values
        GW1[i] = h(tmin,1.0,90.0,300.1,0);
        GW2[i] = h(tmin,0.8,31.0,1000.5,0);
        GW3[i] = h(tmin, 0.1,17.0,100.2,3.1415);
        tmin++;
    }
    

    
    ofstream All;  //np.savetxt replacement
    All.open("All3Waves.dat");   //filename where it shall be saved
    for(int i = 0; i < points; i++)
    {
        All << time[i] << " " << GW1[i] << " " << GW2[i] << " " << GW3[i] << endl;  //pipes: again shuffle the data in the wanted formating into output stream "All". "All" creates the textfile "All3Waves.dat"
    }
    All.close(); //close the textfile
    
    /*
     * By now we have noise-free gravitational waves.
     * Now we add Gaussian noise.
     */
    
    
    
    gsl_rng * r = gsl_rng_alloc (gsl_rng_default);  //replacement for np.random
    double standard_deviation = 0.2;
    
    for(int i = 0; i < points;i++)
    {
        GW1[i] += gsl_ran_gaussian (r, standard_deviation); //replacement of np.random.normal
        GW2[i] += gsl_ran_gaussian (r, standard_deviation);
        GW3[i] += gsl_ran_gaussian (r, standard_deviation);
    }
    
    
    ofstream All2;  //write the noisy gravitational waves into a new file
    All2.open("AllWithNoise.dat");
    for(int i = 0; i < points; i++)
    {
        All2 << time[i] << " " << GW1[i] << " " << GW2[i] << " " << GW3[i] << endl;
    }
    All2.close();
    
    gsl_rng_free(r); //releasing memory (it is a C-thing, ignore it, if you come from python, it has nothing to do with statistics)
    
    
    
    /*
     * Now we search for the gravitational waves
     * Here I put in all the right numbers apart from collapse time.
     * I then output the SN amplitude
     */
    
    ofstream SN_stream;
    SN_stream.open("SN_collapsetime.dat");
    for(int t = 5; t < 1200; t++)
    {
      SN_stream << t << " " <<  Filter_output(GW1, time, standard_deviation, 1.0,90.0,t,0) << " " << Filter_output(GW2,time, standard_deviation, 0.8,31.0,t,0) << " " << Filter_output(GW3,time, standard_deviation, 0.1,17.0, t ,3.1415) << endl;
        
    }
    
    SN_stream.close();
    
    

}

