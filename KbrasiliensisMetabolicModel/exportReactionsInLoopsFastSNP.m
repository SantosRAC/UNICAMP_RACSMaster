% Function that iterate through loop laws and print reactions in each loop
fid=fopen('outReactionsInLoops03062017_3.txt','w');
for i = 1:355
    a = transpose(model_Kbr.rxns(find(Nsnp([1:end],i))));
    b = transpose(find(Nsnp([1:end],i)));
    for j = 1:length(a)
        fprintf(fid, '%d\t%s\t%d\n', i, a{j}, b(1,j));
    end
end
fclose(fid);
