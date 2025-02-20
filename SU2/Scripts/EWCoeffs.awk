BEGIN{FS="|"}{if(NF==7){CD=$5;CL=$6}}END{print alpha,M,CD,CL}
