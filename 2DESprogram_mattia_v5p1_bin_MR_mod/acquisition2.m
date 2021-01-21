function acquisition2(hObject, eventdata, handles)
Data = [];
PhInt = [];
set(hObject,'BackgroundColor','red');
pat = get(handles.path,'String');

if ispc
    a = regexp(pat,'\');
else
    a = regexp(pat,'/');
end

ind = pat(1:a(end));
files = dir(pat);
n = length(files);

cont = 1;
 for ii = 1:n
    if (isempty(findstr(files(cont).name,'reference')) == 0 ||...
            isempty(findstr(files(cont).name,'clean')) == 0) % (isempty(findstr(files(cont).name,'refl')) == 0 ||...
        
            files(cont) = [];
    else
            cont = cont+1;
    end
cont;
 end   

n = length(files);
T2 = zeros(n,1);
NumScans = zeros(n,1);

%%%%%%Ordering filename T2%%%%%%%%%
for ii = 1:n
    fprintf(['Number ' num2str(ii) '\n'])
    pati = [ind,files(ii).name];
    c = strfind(pati,'_');
    d = strfind(pati,'fs');
    T2(ii) = str2num(pati(c(end)+1:d(end)-1));
    %%%%%%%%%%%%%
    %     binfile = fopen(pati,'r');
    %         dass = fread(binfile);
    dass = read_bin(pati);
    fprintf(['t_2 = ' num2str(T2(ii)) ', dass = ' num2str(size(dass)) '\n'])
    NumScans(ii) = dass(2,1);
    Data{ii} = dass(4:end,3:end);
    PhInt{ii} = dass(4:end,2);
    
    if ii == 1
        wpr = dass(1,3:end);
        mStep = dass(4:end,1);
    end

    %%%%%%%%%%%
end

[T2, index] = sort(T2);
%files = files(index);
%NumScans = zeros(n,1);
NumScans = NumScans(index);
Data = Data(index);
PhInt = PhInt(index);
textLabel = sprintf('%d', n);
set(handles.Nfiles, 'String', textLabel);
%%%%%Acquisition data and creation of file.mat%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% for ii = 1:n
%    
%     pati = [ind,files(ii).name];
%     dass = dlmread(pati);
%     NumScans(ii) = dass(2,1);
%     Data{ii} = dass(4:end,3:end);
%     PhInt{ii} = dass(4:end,2);
%     if ii==1
%         wpr = dass(1,3:end);
%         mStep = dass(4:end,1);
%     end
%     
% end
% toc
%%%%%%%%%%%%%%%%%%%%%%%
set(hObject,'BackgroundColor','white')
filemat = [pati(1:(c(end)-1)) '.mat'];

save(filemat,'Data','T2','PhInt','wpr','mStep','NumScans','-v7.3')
set(handles.pathmat,'String',filemat);
end

