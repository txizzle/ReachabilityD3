for t=1:127
    for e=1:86
        fn = strcat('t',int2str(t), 'z', int2str(e), '.csv'); 
        csvwrite(fn, value(:,:,e,t));
    end
end