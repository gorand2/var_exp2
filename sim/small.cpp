#include <cmath>
#include <cstdio>

class small_cyl {
	public:
  double f;
  double V;
  double T;
  double T_env;
  double P_env;
  double n;
  double V_big;
  double P_big;
  double T_big;
  double n_big;
  double fill(double V_frac);
  double search_vfrac(bool &err,double mul);
  double push(bool &err);
  double large_compress(bool &err);
};


double small_cyl::fill(double V_frac) {
  n=V_frac*V*P_env/(f*0.5*8.41*T_env);
  double inv=V_frac*pow(T_env,f/2.0);
  T=pow(inv,2.0/f);
  return f/2.0*n*8.41*(T_env-T)+P_env*V*V_frac;
}

double small_cyl::push(bool &err) {
  err=false;
  /*assuming temperature equalises upon valve opening */
  double T_new_start=(n*T+n_big*T_big)/(n+n_big);
  double inv=(V+V_big)*pow(T_new_start,f/2.0);
  double T_new=pow(inv/V_big,2.0/f);
  double P_new=f/2*(n+n_big)*8.41*T_new/V_big;
  if (P_new>P_env) {err=true; return 0.0;};
  P_big=P_new;
  T_big=T_new;
  n_big+=n;
  return f/2.0*n*8.41*(T_new_start-T_new);
}

double small_cyl::search_vfrac(bool &err,double mul) {
  err=false;
  double inv=P_env*1.0;
  double ex=f/(f+2.0);
  double V_new=pow(inv/(P_big+mul),ex);
  if ((P_big+mul)>=P_env) { err=true; return 0.0; }
  return 1.0/V_new;
}

double small_cyl::large_compress(bool &err) {
  err=false;
  double inv=P_big*pow(V_big,(f+2)/f);
  double V_new_1=pow(inv/P_env,f/(f+2));
  double inv2=V_big*pow(T_big,f/2.0);
  double T_new_1=pow(inv2/V_new_1,2.0/f);
  if (P_env<P_big) { err=true; return 0.0; }
  return -(T_new_1-T_big)*n_big*f*0.5*8.41-V_new_1*P_env;
}

double sim(double V_fr_first,double Vol_rat, double mul) {
  small_cyl cyl;
  cyl.f=5.6;
  cyl.V=1.0;
  cyl.T_env=300.0;
  cyl.P_env=1e5;
  cyl.V_big=Vol_rat;
  cyl.n_big=0.0;
  cyl.P_big=0.0;
  cyl.T_big=0.0;

  double E=0.0,de;
  bool err;
  de=cyl.fill(V_fr_first);
  //printf("fill 1, dE=%lf, E=%lf\n",de,de);
  E+=de;
  de=cyl.push(err);
  //printf("push 1, dE=%lf, E=%lf\n",de,E+de);
  E+=de;
  int n=2;
  do {
      small_cyl cyl1=cyl;
      double fr=cyl.search_vfrac(err,mul);
      double E1=E;
      de=cyl.fill(fr);
      //printf("fill %i, \tdE=%lf, \tE=%lf\n",n,de,E1+de);
      E1+=de;
      de=cyl.push(err);
      //printf("push %i, \tdE=%lf, \tE=%lf\n",n,de,E1+de);
      E1+=de;
      if (err) {
	  cyl=cyl1;
	  //printf("disregard last\n");
      } else {
	  E=E1;
      }
      n++;
  } while(!err);
  de=cyl.large_compress(err);
  //printf("push long, dE=%lf, E=%lf\n",de,E+de);
  return E/-de;
}



int main(int argc, char *argv[]) {
  small_cyl cyl;
  cyl.f=5.6;
  cyl.V=1.0;
  cyl.T_env=300.0;
  cyl.P_env=1e5;
  cyl.V_big=40.0;
  cyl.n_big=0.0;
  cyl.P_big=0.0;
  cyl.T_big=0.0;

  double fill=0.6;
  double mul=700.0;
  if (argc==2) {
      double V_fr=0.01;
      double V_rat=3.0;
      int i,j;
      double score=0;
      double sc;
      for(i=0;i<90;i++) {
	  V_rat=3.0;
          for(j=0;j<190;j++) {
              sc=sim(V_fr,V_rat,mul);
	      if (sc>score) {
		  fill=V_fr;
		  cyl.V_big=V_rat;
		  score=sc;
	      }
	      V_rat+=0.3;
	  }
	  V_fr+=0.01;
      }
  }

  if (argc==3) {
      sscanf(argv[1],"%lf",&fill);
      sscanf(argv[2],"%lf",&cyl.V_big);
  }

  printf("initial fill ratio %lf, volume ratio %lf\n",fill,cyl.V_big); 
  double E=0.0,de;
  bool err;
  de=cyl.fill(fill);
  printf("fill 1, dE=%lf, E=%lf\n",de,de);
  E+=de;
  de=cyl.push(err);
  printf("push 1, dE=%lf, E=%lf\n",de,E+de);
  E+=de;
  int n=2;
  do {
      small_cyl cyl1=cyl;
      double fr=cyl.search_vfrac(err,mul);
      double E1=E;
      de=cyl.fill(fr);
      printf("fill %i, \tdE=%lf, \tE=%lf\n",n,de,E1+de);
      E1+=de;
      de=cyl.push(err);
      printf("push %i, \tdE=%lf, \tE=%lf\n",n,de,E1+de);
      E1+=de;
      if (err) {
	  cyl=cyl1;
	  printf("disregard last\n");
      } else {
	  E=E1;
      }
      n++;
  } while(!err);
  de=cyl.large_compress(err);
  printf("push long, dE=%lf, E=%lf\n",de,E+de);
}
