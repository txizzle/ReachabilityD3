for e=1:86
    for t=1:127
        xx = 1:1:31;
        yy = 1:1:31;
        [C,~] = contourf(xx,yy,value(:,:,e,t), [0 0]);
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
        d3ify(gcf, fn);
    end
end