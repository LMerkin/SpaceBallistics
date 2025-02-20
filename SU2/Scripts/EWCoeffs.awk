BEGIN{FS="|"}{if(NF==8){CD=$6;CL=$7}}END{print alpha,M,CD,CL}
