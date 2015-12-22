#define rand01 (0.9999999*double(rand())/RAND_MAX) 
#define getrandom(max1) (((rand())%(max1))) // random integer between 0 and max-1

#define in_lh(idx) (((idx)<Ne/2) || (Ne<=(idx) && (idx)<(Ne+Ni/2)))
#define in_rh(idx) ((Ne/2<=(idx) && (idx)<Ne) || ((Ne+Ni/2)<=(idx)))

