function Save2file(path, data)
    fid=fopen(path,'wt'); %写入的文件，各函数后面有说明?
    [m,n]=size(data); 
    for i=1:1:m 
        for j=1:1:n 
            if j==n 
                fprintf(fid,'%g\n',data(i,j)); 
            else 
                fprintf(fid,'%g\t',data(i,j)); 
            end
        end
    end
    fclose(fid);
end