function acquisition(hObject, eventdata, handles)
Data=[];
PhInt=[];
set(hObject,'BackgroundColor','red');
pat=get(handles.path,'String');
a=regexp(pat,'\');
ind=pat(1:a(end));
files =dir(pat);
n=length(files);
cont=1;
for i=1:n
    if (isempty(findstr(files(cont).name,'refl'))==0 || isempty(findstr(files(cont).name,'reference'))==0 || isempty(findstr(files(cont).name,'clean'))==0)
            files(cont)=[];
    else
            cont=cont+1;
    end
cont;
end
    
n=length(files);
T2=zeros(n,1);
NumScans=zeros(n,1);
%%%%%%Ordering filename T2%%%%%%%%%
for i=1:n
    pati=[ind,files(i).name];
    c=strfind(pati,'_');
    d=strfind(pati,'fs');
    T2(i)=str2num(pati(c(end)+1:d-1));
    %%%%%%%%%%%%%
        dass=dlmread(pati);
     NumScans(i)=dass(2,1);
     Data{i}=dass(4:end,3:end);
     PhInt{i}=dass(4:end,2);
     if i==1
         wpr=dass(1,3:end);
         mStep=dass(4:end,1);
     end
     %%%%%%%%%%%
end
[T2, index]=sort(T2);
%files=files(index);
%NumScans=zeros(n,1);
NumScans=NumScans(index);
Data=Data(index);
 PhInt=PhInt(index);
textLabel = sprintf('%d', n);
set(handles.Nfiles, 'String', textLabel);
%%%%%Acquisition data and creation of file.mat%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% for i=1:n
%    
%     pati=[ind,files(i).name];
%     dass=dlmread(pati);
%     NumScans(i)=dass(2,1);
%     Data{i}=dass(4:end,3:end);
%     PhInt{i}=dass(4:end,2);
%     if i==1
%         wpr=dass(1,3:end);
%         mStep=dass(4:end,1);
%     end
%     
% end
% toc
%%%%%%%%%%%%%%%%%%%%%%%
set(hObject,'BackgroundColor','white')
filemat=[pati(1:(c(end)-1)) '.mat'];

save(filemat,'Data','T2','PhInt','wpr','mStep','NumScans')
set(handles.pathmat,'String',filemat);