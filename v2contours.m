function v2contours(data)
    mkdir('data_contours_pt2');
    dimensions = size(size(data));
    is4d = (dimensions(2) == 5);
    
    fileID = fopen('data_contours/twoislands.txt','w');
    len_x = length(data(:,1,1,1,1));
    len_y = length(data(1,:,1,1,1));
    len_3 = length(data(1,1,:,1,1));
    len_4 = length(data(1,1,1,:,1));
    if (is4d)
        len_5 = length(data(1,1,1,1,:));
    end
    
    for var1=23:len_3
        disp(strcat('Currently processing e=', int2str(var1))) 
        for var2=1:len_4
            if (is4d)
                for var3=1:len_5
                    xx = 1:1:len_x;
                    yy = 1:1:len_y;
                    [C,~] = contourf(xx,yy,data(:,:,var1,var2,var3), [0 0]);
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
                    fn = strcat('t',int2str(var3), 'd', int2str(var2), 'z', int2str(var1));
                    full_fn = strcat('data_contours/', fn);
                    if (length(breakpts)-1 > 1)
                        fprintf(fileID,strcat(fn, '\n'));
                    end
                    d3ify(gcf, full_fn);
                    close all;
                end
            else
                xx = 1:1:len_x;
                yy = 1:1:len_y;
                
                [C,~] = contourf(xx,yy,data(:,:,var1,var2), [0 0]);
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
                fn = strcat('t',int2str(var2), 'z', int2str(var1));
                full_fn = strcat('data_contours/', fn);
                if (length(breakpts)-1 > 1)
                    fprintf(fileID,strcat(fn, '\n'));
                end
                d3ify(gcf, full_fn);
                close all;
            end
        end
    end
end