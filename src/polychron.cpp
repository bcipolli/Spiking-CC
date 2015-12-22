#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h>
#include <time.h>

#include <sys/stat.h> 

#include "defs.h"
#include "const.cpp"



/********************************************/

double	a[N], d[N];				// parameters for izhikevich neuron
int		post[N][M];				// neuron #s of post-synaptic connections
double	s[N][M], sd[N][M];	  	// synaptic strength and synaptic change
short	delays[N][D][M];		// synaptic transmission delay, indexed by source neuron and delay
short	delays_length[N][D];  	//
int		N_pre[N], I_pre[N][3*M], D_pre[N][3*M];	   //
double	*s_pre[N][3*M], *sd_pre[N][3*M];
double	LTP[N][1001+D], LTD[N];		  //
double	v[N], u[N];	  //
int		N_firings;
int		firings[N_firings_max][2];


bool FileExists(char* strFilename) {
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
	
	// Attempt to get the file attributes
	intStat = stat(strFilename,&stFileInfo);
	if(intStat == 0) {
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	} else {
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}
	
	return(blnReturn);
}

void initialize()
{	int i,j,k,jj,dd, exists, r;

//	a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
	for (i=0;i<Ne;i++) a[i]=0.02;
	for (i=Ne;i<N;i++) a[i]=0.1;

//	d=[   8*ones(Ne,1);      2*ones(Ni,1)];
	for (i=0;i<Ne;i++) d[i]=8;
	for (i=Ne;i<N;i++) d[i]=2;


//	post=ceil([Ne+Ni*rand(Ne,M*Ni/N),Ne*rand(Ne,M*Ne/N);Ne*rand(Ni,M)]); 
	//LH excitatory
	for (i=0;i<Ne/2;i++)
	{
		//intrahemispheric connections
		for (j=0;j<M_INTRA;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(N*rand01));
				if (!in_lh(r)) exists = 1; // intrahemispheric cxn
				if (r == i)    exists = 1;  // cannot connect to self
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1; //do not duplicate connection
			}
			while (exists == 1);
			post[i][j]=r;
		}
		//interhemispheric connections
		for (j=M_INTRA;j<M;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(N*rand01));
				if (!in_rh(r)) exists = 1; // intrahemispheric cxn
				if (r == i)    exists;  // cannot connect to self
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1; //do not duplicate connection
				
			}
			while (exists == 1);
			post[i][j]=r;
		}
	}
	//LH inhibitory: only intrahemispheric
	for (i=Ne;i<(Ne+Ni/2);i++)
	{
		for (j=0;j<M;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(Ne/2*rand01));  
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1;
			}
			while (exists == 1);
			post[i][j]=r; 
		}
	}

	//RH excitatory
	for (i=Ne/2;i<Ne;i++)
	{
		//intrahemispheric connections
		for (j=0;j<M_INTRA;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(N*rand01));
				if (!in_rh(r)) exists = 1; // intrahemispheric cxn
				if (r == i)    exists = 1;    // cannot connect to self
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1; //do not duplicate connection
			}
			while (exists == 1);
			post[i][j]=r;
		}
		//interhemispheric connections
		for (j=M_INTRA;j<M;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(N*rand01));
				if (!in_lh(r)) exists = 1; // intrahemispheric cxn
				if (r == i)    exists = 1;    // cannot connect to self
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1; //do not duplicate connection
			}
			while (exists == 1);
			post[i][j]=r;
		}
	}
	//RH inhibitory: only intrahemispheric
	for (i=(Ne+Ni/2);i<N;i++)
	{
		for (j=0;j<M;j++) 
		{
			do
			{
				exists = 0;
				r = Ne/2+int(floor(Ne/2*rand01));  
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1;
			}
			while (exists == 1);
			post[i][j]=r; 
		}
	}




//	s=[6*ones(Ne,M);-5*ones(Ni,M)];     % initial synaptic weights
	for (i=0;i<Ne;i++) for (j=0;j<M;j++) s[i][j]=6;
	for (i=Ne;i<N;i++) for (j=0;j<M;j++) s[i][j]=-5;
	
//	sd=zeros(N,M);                      % derivatives of synaptic weights
  	for (i=0;i<N;i++)for (j=0;j<M;j++) sd[i][j]=0;
	

//	for i=1:N
  	for (i=0;i<N;i++)
	{
		short ind=0;

//		if i<=Ne
		if (i<Ne)
		{

//			for j=1:D
//            delays{i,j}=M/D*(j-1)+(1:M/D);
//			end;
			for (j=0;j<D;j++) delays_length[i][j] = 0;
			
			for (j=0;j<D_INTRA;j++) 
				delays_length[i][j]+=M_INTRA/D_INTRA;
			for (j=(D_INTER-1);j<D_INTER;j++)
				delays_length[i][j] += M_INTER;

			for (j=0;j<D;j++)
				for (k=0;k<delays_length[i][j];k++)
					delays[i][j][k]=ind++;

			// Make sure the interhemispheric delays match up with
			//   the very last indices.
			if (D_INTER < D)
				cout << "Error: interhemispheric delay < intrahemispheric delay NYI";
				
		}
//		else
//        delays{i,1}=1:M;
//		end;
		else
		{
			for (j=0;j<D;j++) delays_length[i][j]=0;
			delays_length[i][0]=M;
			for (k=0;k<delays_length[i][0];k++)
					delays[i][0][k]=ind++;
		}
	}
	
//		pre{i} = find(post==i & s>0);     % Indeces of pre excitatory neurons
//		aux{i} = N*(D-1-ceil(ceil(pre{i}/N)/(M/D))) + 1+mod(pre{i}-1,N);
  	for (i=0;i<N;i++)
	{
		N_pre[i]=0;
		for (j=0;j<Ne;j++)
		for (k=0;k<M;k++)
		if (post[j][k] == i) 
		{
			I_pre[i][N_pre[i]]=j;
			for (dd=0;dd<D;dd++)
				for (jj=0;jj<delays_length[j][dd];jj++)
					if (post[j][delays[j][dd][jj]]==i) D_pre[i][N_pre[i]]=dd;
			s_pre[i][N_pre[i]]=&s[j][k];
			sd_pre[i][N_pre[i]++]=&sd[j][k];
		}

//	end;
	}

//	LTP = zeros(N,1001+D);
	for (i=0;i<N;i++)
		for (j=0;j<1+D;j++)
			LTP[i][j]=0;

//	LTD = zeros(N,1);
	for (i=0;i<N;i++)	LTD[i]=0;

//	v = -65+10*rand(N,1);               % initial values for v
	for (i=0;i<N;i++)	v[i]=-65+10*rand01;

//	u = 0.2.*v;                         % initial values for u
	for (i=0;i<N;i++)	u[i]=0.2*v[i];

//	firings=[-D 0];                     % spike timings
	N_firings=1;
	firings[0][0]=-D;
	firings[0][1]=0;
}



void load_all(char fname[30])
{
	fprintf(stdout, "Loading data...");
	
	int		i,j, k, Np;
	float  x;
	int		dd;
	
	FILE	*stream;
	stream = fopen( fname, "r" );
    if( stream == NULL )
		cout << " \n\tError: The file " << fname << " cannot be opened" << "\n";
    else
    {
	  /* Set pointer to beginning of file: */
      fseek( stream, 0L, SEEK_SET );
	  for (i=0; i < N; ++i)
	  {
		fscanf( stream, "%f", &x);
		v[i]=x;
		fscanf( stream, "%f", &x);
		u[i]=x;

		for (j=0; j < M; ++j)
		{
			fscanf( stream, "%d", &dd);
			post[i][j]=dd;
			fscanf( stream, "%f", &x);
			s[i][j]=x;
			fscanf( stream, "%f", &x);
			sd[i][j]=x;
		}
		for (k=0; k < D; ++k)
		{
			fscanf( stream, "%d", &dd);
			delays_length[i][k]=dd;
			for (j=0; j < delays_length[i][k]; ++j)
			{
				fscanf( stream, "%d", &dd);
				delays[i][k][j]=dd;
			}
		}

		fscanf( stream, "%d", &dd);
		N_pre[i] = dd;
	    for (j=0; j < N_pre[i]; ++j)
		{
		  fscanf( stream, "%d", &dd);
		  I_pre[i][j]=dd;
		  fscanf( stream, "%d", &dd);
		  D_pre[i][j]=dd;
		}

		fscanf( stream, "%f", &x);
		LTD[i]=x;
		for (j=0; j < D+1; ++j)
		{
			fscanf( stream, "%f", &x);
			LTP[i][j]=x;
		}
	  }
	  
	  fscanf( stream, "%d", &dd);
	  N_firings=dd;
	  for (i=0; i < N_firings; ++i)
	  {
			fscanf( stream, "%d", &dd);
			firings[i][0]=dd;
			fscanf( stream, "%d", &dd);
			firings[i][1]=dd;
	  }
	
	  fclose( stream );

	  for (i=0; i < N; ++i)
	  {
		for (Np=0;Np<N_pre[i];Np++)
		{
			j = I_pre[i][Np];
			for (k=0;k<M;k++)
			if (post[j][k] == i) 
			{
				s_pre[i][Np]=&s[j][k];
				sd_pre[i][Np++]=&sd[j][k];
			}
		}
	  }
	  
	}
	
	fprintf(stdout, "done.\n");
}






//--------------------------------------------------------------
int			N_polychronous;


double		C_rel = 0.95*C_max;

const int	polylenmax = N;

int			N_postspikes[polylenmax], I_postspikes[polylenmax][N], J_postspikes[polylenmax][N], D_postspikes[polylenmax][N], L_postspikes[polylenmax][N];
double		C_postspikes[polylenmax][N];
int			N_links, links[2*W*polylenmax][4];
int			group[polylenmax], t_fired[polylenmax], layer[polylenmax];
int			gr3[W], tf3[W];
int			I_my_pre[3*M], D_my_pre[3*M], N_my_pre;
int			N_fired;


FILE		*fpoly_summary;
FILE		*fpoly_spikes;
FILE		*fpoly_links;





// The new (default) algorithm to find polychronous groups

const	int	max_latency = D; // maximum latency 


//--------------------------------------------------------------
void	polychronous(int nnum)
{
	int	i,j, t, p, k;
	int npre[W];
	int dd;
	int	t_last, timing;
	int	Dmax, L_max; 
	int	used[W], discard;

	double v[N],u[N],I[N];
	
	cout << "Neuron # " << nnum << ":" << endl;
	//fprintf(stdout, "a\n"); fflush(stdout);
	
	N_my_pre = 0;
	for (i=0;i<N_pre[nnum];i++) {
		//fprintf(stdout, "%d/%d (%d); M=%d: %d\n", i, N_pre[nnum], N_my_pre, M, s_pre[nnum][i]); fflush(stdout);
		
		if (*s_pre[nnum][i] > C_rel) {
			I_my_pre[N_my_pre]=I_pre[nnum][i];
			D_my_pre[N_my_pre]=D_pre[nnum][i];
			N_my_pre++;
		}
	}
	//fprintf(stdout, "b.a\n"); fflush(stdout);
	if (N_my_pre<W) return;

	//fprintf(stdout, "b\n"); fflush(stdout);
	
	for (i=0;i<W;i++)	npre[i]=i;

	// Find anatomical groups
	while (0==0) 
	{
		Dmax=0;
		for (i=0;i<W;i++) if (Dmax < D_my_pre[npre[i]]) Dmax=D_my_pre[npre[i]];
		
		for (i=0;i<W;i++)
		{
			group[i]=I_my_pre[npre[i]];
			t_fired[i]= Dmax-D_my_pre[npre[i]];
			layer[i]=1;
			
			for (dd=0; dd<D; dd++)		 
			for (j=0; j<delays_length[group[i]][dd]; j++)
			{
				p = post[group[i]][delays[group[i]][dd][j]];
				if ((s[group[i]][delays[group[i]][dd][j]] > C_rel) & (dd>=D_my_pre[npre[i]]))
				{
					timing = t_fired[i]+dd+1;
					J_postspikes[timing][N_postspikes[timing]]=group[i];				// presynaptic
					D_postspikes[timing][N_postspikes[timing]]=dd;					// delay index
					C_postspikes[timing][N_postspikes[timing]]=s[group[i]][delays[group[i]][dd][j]];	// syn weight
					I_postspikes[timing][N_postspikes[timing]++]=p;						// index of post target	
				}
			}
		}
	
		for (i=0;i<N;i++) {v[i]=-70; u[i]=0.2*v[i]; I[i]=0;};

		N_links = 0;
		N_fired=W;
		t_last = D+D+max_latency+1; //?
		t=-1;
		while ((++t<t_last) & (N_fired < polylenmax))
		{
			for (p=0;p<N_postspikes[t];p++) 
			  I[I_postspikes[t][p]]+=C_postspikes[t][p]; 
 
		  	for (i=0;i<N;i++)
			{
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				u[i]+=a[i]*(0.2*v[i]-u[i]);
				I[i]=0;
			}

			for (i=0;i<N;i++) 
			{
				if (v[i]<30) continue;

				v[i] = -65;
				u[i]+=d[i];

				if (N_fired >= polylenmax) break;
				
				t_fired[N_fired]= t;
				group[N_fired++]=i;
				for (dd=0; dd<D; dd++)
				{
					for (j=0; j<delays_length[i][dd]; j++)
					{
						if ((s[i][delays[i][dd][j]] <= C_rel) && (i<Ne)) continue;

						timing = t+dd+1;
						if (timing >= polylenmax)
						{	cout << "timing > polylenmax" << endl; continue; }
						
						J_postspikes[timing][N_postspikes[timing]]=i;				// presynaptic
						D_postspikes[timing][N_postspikes[timing]]=dd;				// delay+1
//							L_postspikes[timing][N_postspikes[timing]]=NL+1;			// layer
						C_postspikes[timing][N_postspikes[timing]]=s[i][delays[i][dd][j]];	   // syn weight
						I_postspikes[timing][N_postspikes[timing]++]=post[i][delays[i][dd][j]];// index of post target	
					}

					if (t_last < timing+1) 
					{
						t_last = timing+1;
						if (t_last > polylenmax-D-1) t_last = polylenmax-D-1;
					}
				} // for dd
			} // for i
		}
		
		if (N_fired>2*W)
		{
			N_links=0;
			L_max=0;
			for (i=W;i<N_fired;i++)
			{
				layer[i]=0;
				for (p=t_fired[i]; (p>t_fired[i]-max_latency) & (p>=0); p--)
				{
					for (j=0;j<N_postspikes[p];j++)
					{
						if ((I_postspikes[p][j]!=group[i]) | (J_postspikes[p][j]>=Ne)) continue; 
						if (N_links >= 2*W*polylenmax) continue;

						for (k=0;k<i;k++)
							if ((group[k]==J_postspikes[p][j]) & (layer[k]+1>layer[i])) layer[i]=layer[k]+1;
							
						links[N_links][0]=J_postspikes[p][j];  //presynaptic neuron (output column #3)
						links[N_links][1]=I_postspikes[p][j];  //index of post-synaptic (output column #4)
						links[N_links][2]=D_postspikes[p][j];  // delay INDEX
						links[N_links++][3]=layer[i];          //current length, as counted to determine max_length
					   	if (L_max < layer[i]) L_max = layer[i]; 
					}
				}
			}

			discard = 0;
			for (i=0;i<W;i++)
			{
				used[i]=0;
				for (j=0;j<N_links;j++) if ((links[j][0] == group[i]) & (links[j][1] < Ne)) used[i]++;
				if (used[i] == 1) discard = 1;
			}

//			if ((discard == 0) & (t_fired[N_fired-1] > min_group_time) )  // (L_max >= min_group_path))
			if ((discard == 0) & (L_max >= min_group_path))
			{

				for (i=0;i<W;i++) {gr3[i]=group[i]; tf3[i]=t_fired[i];};

				N_polychronous++;
				
				cout << "\tN_polychronous= " << N_polychronous+POLY_START_NUMBER-1;
				cout << ",\tN_fired = " << N_fired;
				cout << ",\tL_max = " << L_max;
				cout << ",\tT=" << t_fired[N_fired-1];
				cout << endl;
				
				if (fpoly_summary)
				{	fprintf(fpoly_summary, "%d\t%d\t%d\t%d\n", nnum, N_fired, L_max, t_fired[N_fired-1]); 
					fflush(fpoly_summary); }
					
				if (fpoly_spikes)
				{	for (i=0; i<N_fired; i++)
						fprintf(fpoly_spikes, "%d\t%d\t%d\t%d\n", nnum, N_polychronous+POLY_START_NUMBER-1, group[i], t_fired[i]);
					fflush(fpoly_spikes); }
					
				if (fpoly_links)
				{	for (j=0;j<N_links;j++)
				   		fprintf(fpoly_links, "%d\t%d\t%d\t%d\t%d\t%d\n", 
				   		nnum, 
				   		N_polychronous+POLY_START_NUMBER-1, 
				   		links[j][0], 
				   		links[j][1], 
				   		links[j][2], 
				   		links[j][3]);  
					fflush(fpoly_links); }
			}
		}

  		for (dd=Dmax;dd<t_last;dd++) N_postspikes[dd]=0;
		if (t_last == polylenmax-D) for (dd=t_last;dd<polylenmax;dd++) N_postspikes[dd]=0;

		i=1;
		while (++npre[W-i] > N_my_pre-i) if (++i > W) return; 
		while (i>1) {npre[W-i+1]=npre[W-i]+1; i--;}
	}	
}


 
//--------------------------------------------------------------
void	all_polychronous()
{
	int	i;
	N_polychronous=0;
   	for (i=0;i<polylenmax;i++) N_postspikes[i]=0;

	for (i=POLY_START_NEURON;i<POLY_END_NEURON;i++) polychronous(i);
fprintf(stdout,"xx\n"); fflush(stdout);

	cout << "\nN_polychronous=" << (N_polychronous+POLY_START_NUMBER-1) << "\n";
}




void shuffle()
{
	int i, j, ri, rj;
	double x,y;
	cout << "***** scrambling ****";
	for (i=0;i<Ne;i++)
	for (j=0;j<M;j++)
	if (post[i][j] < Ne)
	{
		ri = int(floor(rand01*Ne));
		do 
		{
			rj = int(floor(rand01*M));
		}
		while (post[ri][rj] >= Ne);	 
		x=s[ri][rj];
		y=sd[ri][rj];
		s[i][j]=s[ri][rj];
		sd[i][j]=sd[ri][rj];
		s[ri][rj]=x;
		sd[ri][rj]=y;
	}
}


// --------------------------------------------------------------------------
int main()
{
	int k;

	srand(RSEED); 
	initialize();

	// Load from a previous run
	if (FileExists((char*)poly_results_fn))
		load_all((char*)poly_results_fn);

	// Cannot initialize to spnet results
	else if (!FileExists((char*)sp_full_fn)) {
		fprintf(stdout, "Error: %s does not exist.  Run spnet first.\n", sp_full_fn);
		return 1; 
	}
	
	// Initialize to spnet results and find poly groups
	else {
		load_all((char*)sp_full_fn);
	}
	
	// Pull out polychronous groups (saving to files)
//	if (FileExists((char*)poly_summary_fn)) remove((char*)poly_summary_fn);
//	if (FileExists((char*)poly_spikes_fn))  remove((char*)poly_spikes_fn);
//	if (FileExists((char*)poly_links_fn))   remove((char*)poly_links_fn);
	fpoly_summary = fopen(poly_summary_fn,"a");
	fpoly_spikes  = fopen(poly_spikes_fn,"a");
	fpoly_links   = fopen(poly_links_fn, "a");
	
	all_polychronous(); k=N_polychronous+POLY_START_NUMBER-1;
	
	fclose(fpoly_summary); fclose(fpoly_spikes); fclose(fpoly_links);
	
	// Get the # of polychronous groups in the shuffled case
	//shuffle();
	//all_polychronous(); 
	//cout << "ratio (true/shuffled) = " << double(k)/(N_polychronous+1) << "\n";

	return 0;
}