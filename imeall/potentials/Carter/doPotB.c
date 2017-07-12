#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define A0 2.8553

#define NFePair 13
const double aFePair[NFePair] = 
  {
    -27.444805994228,
    15.738054058489,
    2.2077118733936,
    -2.4989799053251,
    4.2099676494795,
    -0.77361294129713,
    0.80656414937789,
    -2.3194358924605,
    2.6577406128280,
    -1.0260416933564,
    0.35018615891957,
    -0.058531821042271,
    -0.0030458824556234
  };

const double rFePair[NFePair] =
  { 
    2.2, 2.3, 2.4, 2.5,
    2.6, 2.7, 2.8, 3.0, 
    3.3, 3.7, 4.2, 4.7, 
    5.3
  };


#define NFeRho 3
const double aFeRho[NFeRho] =
  {
    11.686859407970,
    -0.014710740098830,
    0.47193527075943
  };

const double rFeRho[NFeRho] = 
  {
    2.4, 3.2, 4.2
  };


const double aFeEmb2 = -6.7314115586063e-4;
const double aFeEmb4 =  7.6514905604792e-8;


double Fe_embedding (double rho)
{
  return (
          -sqrt(rho) + aFeEmb2*rho*rho
          + aFeEmb4*rho*rho*rho*rho
          );
  
};



double FeFe_electron_density (double r)
{
  double sum = 0.0;
  int i;
  
  for(i=NFeRho-1; i>=0; i--){
    if (r> rFeRho[i])
      break;
    else
      sum += aFeRho[i]*(rFeRho[i]-r)*(rFeRho[i]-r)*(rFeRho[i]-r);
  }
  
  return(sum);
}



#define r1  1.0
#define r2  2.05
#define B0  7.4122709384068
#define B1  -0.64180690713367
#define B2  -2.6043547961722
#define B3  0.62625393931230


double FePair1 (double r)
{
  double phi;
  if(r<1.e-8)
    phi = 0.0;
  
  else{
    phi = (9.7342365892908E+03/r)*
      (0.1818*exp(-28.616724320005*r)
       +0.5099*exp(-8.4267310396064*r)
       +0.2802*exp(-3.6030244464156*r)
       +0.02817*exp(-1.8028536321603*r)
       );
  }
  
  return(phi);
}



double FePair2 (double r)
{
  return(exp(B0 + r*(B1 + r*(B2 + B3*r))));
}



double Fe_two_body (double r)
{
  double sum = 0.0;
  int i;
  
  if(r<r1)
    sum = FePair1(r);
  
  else if (r1<=r && r<= r2)
    sum = FePair2(r);
  
  else{
    for(i=NFePair-1; i>=0; i--){
      if (r> rFePair[i])
	break;
      else
	sum += aFePair[i]*(rFePair[i]-r)*(rFePair[i]-r)*(rFePair[i]-r);
    }
  }
  
  return(sum);
}


// End of Mendelev's Fe


// Start new H and FeH
const double rc_H = 2.4;
const double rc_FeH = 4.2;


#define NFeH 7
const double rFeH[NFeH] =
  {
    1.6, 1.7, 1.8, 2.0, 
    2.5, 3.2, 4.2
  };

double aFeH[NFeH] =
  {
    14.0786236766230779,
    -4.4526835638887965,
    5.5025349784052979,
    -1.0687331741292405,
    -0.3461226670484926,
    -0.0064991313802717,
    -0.0357322844877736
  };


#define NFeHRho 6
const double rFeHRho[NFeHRho] =
  {
    1.6, 1.8, 2.0, 
    2.4, 3.2, 4.2
  };

double aFeHRho[NFeHRho] =
  {
    10.0073629218346891,
    32.4862873850836635,
    -0.9494211670931015,
    11.6683860903729624,
    -0.0147079871493827,
    0.4945807618408609
  };


#define NHFeRho 5
const double rHFeRho[NHFeRho] =
  {
    1.5, 2.0, 2.5, 
    3.0, 4.2
  };


double aHFeRho[NHFeRho] =
  {
    11.1667357634216433,
    -3.0351469477486712,
    3.6092404272928578,
    0.0212508491354509,
    0.0303904795842773
  };


#define NHemb 6


double aHemb[NHemb] =
  {
    -0.0581047132616673,
    0.0022873205657864,
    -0.0000313966169286,
    0.0000013788174098,
    -0.0000000253074673,
    0.0000000001487789
  };


double rswitch = 0.9;
double C1_H =  0.0;
double C2_H =  0.0;
double C3_H =  0.25;


double H_fc (double r)
{
  double fc;
  fc = expl(1/(r-rc_H));

  return(fc);
}


double FeH_electron_density (double r)
{
  double rho = 0.0;
  int i;
  
  for(i=NFeHRho-1; i>=0; i--){
    if (r> rFeHRho[i])
      break;
    else
      rho += aFeHRho[i]*(rFeHRho[i]-r)*(rFeHRho[i]-r)*(rFeHRho[i]-r);
  }
  
  return (rho);
}



double HFe_electron_density (double r)
{
  double rho = 0.0;
  int i;
  
  for(i=NHFeRho-1; i>=0; i--){
    if (r> rHFeRho[i])
      break;
    else
      rho += aHFeRho[i]*(rHFeRho[i]-r)*(rHFeRho[i]-r)*(rHFeRho[i]-r);
  }
  
  return (rho);
}


double HH_electron_density (double r)
{
  double rho;
  
  const double alpha =  1800.0;

  if(r <= rc_H){
    rho = alpha*r*r*expl(-2.*r/0.53)*H_fc(r);
  }
  else
    rho = 0.0;
  
  return (rho);
}

double H_embedding (double rho)
{
  double emb = rho* (aHemb[0] + 
		     rho*(aHemb[1] + 
			  rho*(aHemb[2] + 
			       rho*(aHemb[3] + 
				    rho*(aHemb[4] + 
					 rho*aHemb[5]
					 )
				    )
			       )
			  )
		     );
  return(emb);
}



// H-H two body
const double Eb_H =   2.37;
const double r0_H =   0.74;
const double lambda = 0.4899;

double H_two_body (double r)
{
  double phi;

  if(r<1.e-6 || r>= rc_H)
    phi = 0.0;

  else{
    double s, Emol, a, rho;
    s = 0.5*(1-tanh(25*(r-rswitch)));
    a = (r-r0_H)/(r0_H*lambda);
    Emol = -2.*Eb_H*(1+a)*expl(-1.*a);
    rho = HH_electron_density(r);
    phi = s*(Emol-2*H_embedding(rho))
      + (1-s)*(C1_H*H_fc(r) + C2_H*rho);
    //phi = (Emol-2*H_embedding(rho));

  }

  return(phi);
}




double FeHPair1 (double r)
{
  double phi;
  if(r<1.e-8)
    phi = 0.0;
  
  else{ 
    phi = (374.39371497272307692307/r)*
      (0.1818*exp(-21.356399424569375*r)
       +0.5099*exp(-6.2887922430536625*r)
       +0.2802*exp(-2.688904165049687*r)
       +0.02817*exp(-1.3454531637478704*r)
       );
  }
  
  return(phi);
}






double FeHPair2 (double r)
{
#define C0   1242.154614241987
#define C1  -6013.4610429013765
#define C2  12339.275191444543
#define C3 -12959.339514470237
#define C4   6817.662603221567
#define C5  -1422.130403271231

  return (C0 + r*(C1 + r*(C2 + r*(C3 + r*(C4 + C5*r)))));
}


double FeH_two_body (double r)
{
  double phi = 0.0;
  int i;
  
  if(r<0.6)
    phi = FeHPair1 (r);

  else if(r>=0.6 && r<=1.2)
    phi = FeHPair2(r);
  
  else{
    for(i=NFeH-1; i>=0; i--){
      if (r> rFeH[i])
	break;
      else
	phi += aFeH[i]*(rFeH[i]-r)*(rFeH[i]-r)*(rFeH[i]-r);
    }
  }
  
  return(phi);
}




#define MAX(a,b) (a>b ? a : b)


int main()
{  
  int nrho=10000, nr=10000;
  double rc = MAX(rFePair[NFePair-1], MAX(rc_FeH, rc_H));
  double dr = rc/(double)nr;
  double rho_max = 300.0;
  double drho = rho_max/(double)nrho;
  
  double emb_Fe[nrho], rho_FeFe[nr], rho_FeH[nr],
    emb_H[nrho], rho_HFe[nr], rho_HH[nr],
    phi_Fe[nr], phi_H[nr], phi_FeH[nr];
  
  int i;
  double x;
  
  // Generate embedding functions
  x = 0.0;
  for(i=0; i<nrho; i++)
    {
      emb_Fe[i] = Fe_embedding (x);
      emb_H[i]  = H_embedding (x);
      x += drho;
    }
  
  
  // Generate electron density & two-body interactions      
  x = 0.0;
  for(i=0; i<nr; i++)
    {
      rho_FeFe[i] = FeFe_electron_density(x);
      rho_FeH[i] = FeH_electron_density(x);
      phi_Fe[i] = x*Fe_two_body(x);
      
      // deal with H
      rho_HFe[i] = HFe_electron_density (x);
      rho_HH[i] = HH_electron_density (x);
      phi_H[i] = x*H_two_body (x);
      
      // Fe-H term
      phi_FeH[i] = x*FeH_two_body(x);
      
      x += dr;
    }
  
  
  // Write potential file
  FILE *fp;
  fp = fopen("PotentialB.fs", "w");
  
  // comments on first 3 lines
  fprintf(fp, "Fe-H EAM potential\n");
  fprintf(fp, "Mendelev Fe + splines for H\n");
  time_t timer;
  timer=time(NULL);
  fprintf(fp, "Tabulated by Ashwin, %s", asctime(localtime(&timer)));
  
  //elements on line 4
  fprintf(fp, "    2   Fe   H\n");
  
  //nrho, drho, nr, dr on line 6
  fprintf(fp, "%5d  % 0.16e  %5d  % 0.16e  %0.16e\n", nrho, drho, nr, dr, rc);
  
  //Fe data on line 7
  fprintf(fp, "    26    55.847   %f   BCC\n", A0);

  // Dump Fe_embedding
  for(i=0; i<nrho/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", emb_Fe[5*i],
	    emb_Fe[5*i+1], emb_Fe[5*i+2], emb_Fe[5*i+3], emb_Fe[5*i+4]);
  
  // Dump FeFe_electron_density
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", rho_FeFe[5*i],
	    rho_FeFe[5*i+1], rho_FeFe[5*i+2], rho_FeFe[5*i+3], rho_FeFe[5*i+4]);

  // Dump FeH_electron_density
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", rho_FeH[5*i],
	    rho_FeH[5*i+1], rho_FeH[5*i+2], rho_FeH[5*i+3], rho_FeH[5*i+4]);
  
  
  // H data on line ...
  fprintf(fp, "     1     1.008   1.8   BCC\n");

  // Dump H_embedding
  for(i=0; i<nrho/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", emb_H[5*i],
	    emb_H[5*i+1], emb_H[5*i+2], emb_H[5*i+3], emb_H[5*i+4]);
  
  // Dump HFe_electron_density
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", rho_HFe[5*i],
	    rho_HFe[5*i+1], rho_HFe[5*i+2], rho_HFe[5*i+3], rho_HFe[5*i+4]);
  
  // Dump HH_electron_density
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", rho_HH[5*i],
	    rho_HH[5*i+1], rho_HH[5*i+2], rho_HH[5*i+3], rho_HH[5*i+4]);
  
  // Dump Fe-Fe pair potential
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", phi_Fe[5*i],
	    phi_Fe[5*i+1], phi_Fe[5*i+2], phi_Fe[5*i+3], phi_Fe[5*i+4]);


  // Dump H-Fe pair potential
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", phi_FeH[5*i],
	    phi_FeH[5*i+1], phi_FeH[5*i+2], phi_FeH[5*i+3], phi_FeH[5*i+4]);

  // Dump H-H pair potential
  for(i=0; i<nr/5; i++)
    fprintf(fp, "% 0.16e  % 0.16e  % 0.16e  % 0.16e  % 0.16e\n", phi_H[5*i],
	    phi_H[5*i+1], phi_H[5*i+2], phi_H[5*i+3], phi_H[5*i+4]);

  fclose(fp);

  
  // Write files for plotting
  
  fp = fopen("pairfuncsMFePs.dat", "w");
  fprintf(fp, "# r  phi_Fe  phi_H  phi_FeH\n");
  x = 0.0;
  for(i=0; i<nr; i++){
    fprintf(fp, "% le % le % le % le\n", x, Fe_two_body(x), H_two_body(x), 
	    FeH_two_body(x) );
    x += dr;
  }
  fclose(fp);
  
  fp = fopen("rhofuncsMFePs.dat", "w");
  fprintf(fp, "# r  rho_Fe  rho_H\n");
  x = 0.0;
  for(i=0; i<nr; i++){
    fprintf(fp, "% le % le % le % le % le\n", x, FeFe_electron_density(x), 
	    FeH_electron_density(x), HFe_electron_density(x), 
	    HH_electron_density(x) );
    x += dr;
  }
  fclose(fp);
  
  fp = fopen("embfuncsMFePs.dat", "w");
  fprintf(fp, "# r  emb_Fe  emb_H\n");
  x = 0.0;
  for(i=0; i<nr; i++){
    fprintf(fp, "% le % le % le\n", x, Fe_embedding(x), H_embedding(x));
    x += drho;
  }
  fclose(fp);
}
