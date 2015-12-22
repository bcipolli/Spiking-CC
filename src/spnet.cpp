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

	// Validate parameters
	if (M_INTRA%D_INTRA != 0) {
		fprintf(stdout, "M_INTRA(%d) not divisible by D_INTRA(%d).\n", M_INTRA, D_INTRA);
		exit(1);
	}

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



// ----------------------------------------------------------------------------------

void save_all(char* fname)
{
	int		i,j,k;
	FILE	*fce;
	fce = fopen(fname,"w");
	cout << "\nsaving data \n";
	
	for (i=0; i < N; ++i)
	{
		fprintf(fce, "%5.3f %5.3f \n", v[i], u[i]);
		for (j=0; j < M; ++j)
			fprintf(fce, "%d %5.3f %6.5f \n", post[i][j], s[i][j], sd[i][j]);
		
		
		for (k=0; k < D; ++k)
		{
			fprintf(fce, "%d  ", delays_length[i][k]);
			for (j=0; j < delays_length[i][k]; ++j)
				fprintf(fce, "%d ", delays[i][k][j]);
			fprintf(fce, "\n");
		}
		
		fprintf(fce, "%d  ", N_pre[i]);
		for (j=0; j < N_pre[i]; ++j)
			fprintf(fce, "%d %d ", I_pre[i][j], D_pre[i][j]);
		
		fprintf(fce, "\n %5.4f ", LTD[i]);
		for (j=0; j < D+1; ++j)
			fprintf(fce, "%5.4f ", LTP[i][j]);
		fprintf(fce, "\n");
	}

	fprintf(fce, " %d", N_firings);
	for (i=0; i < N_firings; ++i)
		fprintf(fce, "%d %d ", firings[i][0],firings[i][1]);

	fclose(fce);
}

void run_all(char* fname)
{
	cout << "Running simulation for " << NHOURS << "hours" << "\n";
	
	int		i,j,k;
	int		sec, t;
	double	I[N];
	FILE	*fs, *fx, *fd;
	
	//	for sec=1:60*60*5
	for (sec=0; sec<60*60*NHOURS; sec++)
	{
		
		
		
		//		for t=1:1000                  % simulation of 1 sec
		for (t=0;t<1000;t++)
		{
			
			for (i=0;i<N;i++) I[i] = 0;
			I[int(floor(N*rand01))]=20;
			
			
			for (i=0;i<N;i++) 
				//			fired = find(v>=30);          % indices of fired neurons
				if (v[i]>=30)
				{
					
					
					//			    v(fired)=-65;
					v[i] = -65;
					
					//	            u(fired)=u(fired)+d(fired);
					u[i]+=d[i];
					
					//	            LTP(fired,t+D)=0.1;
					LTP[i][t+D]= 0.1;
					
					//	            LTD(fired)=0.12;
					LTD[i]=0.12;
					
					//				for k=1:length(fired)
					//					sd(pre{fired(k)}) = sd(pre{fired(k)})+LTP(N*t+aux{fired(k)});
					//				end;
					for (j=0;j<N_pre[i];j++) *sd_pre[i][j]+=LTP[I_pre[i][j]][t+D-D_pre[i][j]-1];
					
					//	            firings=[firings; t+zeros(length(fired),1), fired];
					firings[N_firings  ][0]=t;
					firings[N_firings++][1]=i;
					if (N_firings == N_firings_max)
					{
						cout << "*** Too many spikes, t=" << t << "*** (ignoring)";
						N_firings=1;
					}
				}
			
			//	        k=size(firings,1);
			k=N_firings-1;
			
			//	        while t-firings(k,1)<D
			while (t-firings[k][0] <D)
			{
				
				//	            del=delays{firings(k,2),t-firings(k,1)+1};
				for (j=0; j< delays_length[firings[k][1]][t-firings[k][0]]; j++)
				{
					//					ind = post(firings(k,2),del);
					i=post[firings[k][1]][delays[firings[k][1]][t-firings[k][0]][j]]; 
					
					//					I(ind)=I(ind)+s(firings(k,2), del)';
					I[i]+=s[firings[k][1]][delays[firings[k][1]][t-firings[k][0]][j]];
					
					//					if firings(k,2) <=Ne
					if (firings[k][1] <Ne)
						
						//						sd(firings(k,2),del)=sd(firings(k,2),del)-LTD(ind)';
						sd[firings[k][1]][delays[firings[k][1]][t-firings[k][0]][j]]-=LTD[i];
					//					end;
				}
				
				//				k=k-1;
				k--;
				
				//		    end;
			}
			
			for (i=0;i<N;i++)
			{
				//		        v = v + 0.5*((0.04*v+5).*v+140-u+I);    % for numerical stability
				//			    v = v + 0.5*((0.04*v+5).*v+140-u+I);    % time step is 0.5 ms
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				
				//				u = u + a.*(0.2*v-u);
				u[i]+=a[i]*(0.2*v[i]-u[i]);
				
				//				LTP(:,t+D+1)=0.95*LTP(:,t+D); % tau = 20 ms
				LTP[i][t+D+1]=0.95*LTP[i][t+D];
				
				//				LTD=0.95*LTD;                 % tau = 20 ms
				LTD[i]*=0.95;
			}
			
			//		end;
		}
		
		//		frate(end+1)=sum(firings(:,2)<=Ne)/Ne;
		double	frate=0;
		for (i=1;i<N_firings;i++)
			if ((firings[i][0] >=0) && (firings[i][1] <Ne)) frate++;
		frate = frate/Ne;
		
		//		str(end+1) = sum(sum(s(find(post<=Ne)) > 6.3))/Ne;
		double	str=0;
		for (i=0;i<Ne;i++)
			for (j=0;j<M;j++)
				if ((s[i][j] > 0.9*C_max) && (post[i][j] <Ne)) str++;
		str=100*str/Ne/M;
		
		//		sec, [frate(end), str(end), sum(firings(:,2)>Ne)/Ni]
		double	ifrate=0;
		for (i=1;i<N_firings;i++)
			if ((firings[i][0] >=0) && (firings[i][1] >=Ne)) ifrate++;
		ifrate = ifrate/Ni;
		cout << "sec=" << sec << ",\texc. frate=" << frate << ",\texc->exc str=" << str << ",\tinh. frate=" << ifrate << ".\n";
		fx = fopen(sp_summary_fn,"a");
		fprintf(fx, "%d\t%2.2f\t%2.2f\t%2.2f\n", sec, frate, str, ifrate);
		fflush(fx);
		fclose(fx);
		
		
		//		plot(firings(:,1),firings(:,2),'.');
		//		axis([0 1000 0 N]); drawnow;
   		/*
		fs = fopen(sp_spikes_fn,"a");
		for (i=1;i<N_firings;i++)
			if (firings[i][0] >=0)
				fprintf(fs, "%d\t%d\n", sec*000+firings[i][0], firings[i][1]);
		fclose(fs);
		*/
		//remove("spikes.dat"); 
		//rename( "spikes.tmp.dat", "spikes.dat" );
		
		
		//		LTP(:,1:D+1)=LTP(:,1001:1001+D);
		for (i=0;i<N;i++)
			for (j=0;j<D+1;j++)
				LTP[i][j]=LTP[i][1000+j];
		
		//	    ind = find(1001-firings(:,1) < D);
		k=N_firings-1;
		while (1000-firings[k][0]<D) k--;
		
		//		firings=[-D 0; firings(ind,1)-1000, firings(ind,2)];
		for (i=1;i<N_firings-k;i++)
		{
			firings[i][0]=firings[k+i][0]-1000;
			firings[i][1]=firings[k+i][1];
		}
		N_firings = N_firings-k;
		
		//      sd=0.9*sd;                         % tau = 250 ms
		//      s(1:Ne,:)=max(0,min(7, 0.01+s(1:Ne,:)+sd(1:Ne,:)));
		for (i=0;i<Ne;i++)
			for (j=0;j<M;j++)
			{
				sd[i][j]*=0.9;
				s[i][j]+=0.01+sd[i][j];
				if (s[i][j]>C_max) s[i][j]=C_max;
				if (s[i][j]<0) s[i][j]=0;
			}
		
		//    if mod(sec,10)==0, 
		//        save all; 
		//    end;
		if ( (sec%100==0) & (sec > 0)) 
		{
			save_all(fname);
			
			fs = fopen(sp_conns_fn, "w");
			for (i=0; i<Ne; i++)
				for (j=0;j<M; j++)
					fprintf(fs, "%d\t%3.3f\n", post[i][j], s[i][j]);
			fclose(fs);
		}
		
			
		//	end;
		
	}
}


// --------------------------------------------------------------------------
int main()
{
	int k;
	
	srand(RSEED);
	initialize();
	
	run_all((char*)sp_full_fn);

	return 0;
}