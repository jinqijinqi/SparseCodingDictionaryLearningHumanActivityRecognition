function  Joint= read_joint( fileName)
    fid = fopen(fileName);
    tline = num2str(fgetl(fid));
    found = false;
    Joint=[];
    while ischar(tline)&~strcmp(tline,'END'),
        [fnum,tline] = strtok(tline,',');
        fnum=str2num(fnum);
        found = false;
%         if fnum == frameNum,
            fprintf('found frame num %d\n', fnum);
            data = zeros(11,14);
            data_pos = zeros(4,4);
            for i=1:11,
                for j=1:14,
                    [element,tline] = strtok(tline,',');
                    data(i,j) = str2num(element);
                end
            end
            for i=1:4,
                for j=1:4,
                    [element,tline] = strtok(tline,',');
                    data_pos(i,j) = str2num(element);
                end
            end
            Joint(1:11,:,fnum)=data(1:11,11:13);
            Joint(12:15,:,fnum)=data_pos(:,1:3);
            [element,tline] = strtok(tline,',');

            if strcmp(element, '') ~=1,
                error('ERROR! more data exist... parsing error..');
            end

%             showSkel(data, data_pos, figureTitle);
%             found = true;
%         end

%         if found,
%             break;
%         end

        % read next line
        tline = fgetl(fid);
%         if fnum==1351
%             pause;
%         end
    end
    
%     if ~found,
%         error('NO SUCH FRAME NUMBER EXIST IN THE FILE!');
%     end

end