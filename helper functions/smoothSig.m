function sig=smoothSig(sig)

siz=3;

h=ones(1,siz)/siz;

sig=filter(h,1,sig);
