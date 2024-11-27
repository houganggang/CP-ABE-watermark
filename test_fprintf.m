clear
close all;
for i = 1:10
    itear = i;
    fid = fopen('aaaa.txt','a+' );%'at');
    fprintf(fid, 'itear=%.0f:\r\n',itear);
    fclose(fid);
end