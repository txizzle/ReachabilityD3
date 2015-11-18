function v2contours(data)
    mkdir('data_contours');
    fileID = fopen('data_contours/twoislands.txt','w');
    len_x = length(data(:,1,1,1));
    len_y = length(data(1,:,1,1));
    len_e = length(data(1,1,:,1));
    len_t = length(data(1,1,1,:));
    for e=1:len_e
        disp(strcat('Currently processing e=', int2str(e))) 
        for t=1:len_t
            xx = 1:1:len_x;
            yy = 1:1:len_y;
            [C,~] = contourf(xx,yy,data(:,:,e,t), [0 0]);
            Csplit={};
            breakpts = find(C(1,:)==0);
            breakpts = [breakpts, size(C,2)+1];
            for j=1:length(breakpts)-1
                Csplit{j} = C(:,breakpts(j)+1:breakpts(j+1)-1);
            end
            figure, hold on;
            for j=1:length(Csplit)
                plot(Csplit{j}(1,:),Csplit{j}(2,:));
            end
            fn = strcat('t',int2str(t), 'z', int2str(e));
            full_fn = strcat('data_contours/', fn);
            if (length(breakpts)-1 > 1)
                fprintf(fileID,strcat(fn, '\n'));
            end
            d3ify(gcf, full_fn);
            close all;
        end
    end
end