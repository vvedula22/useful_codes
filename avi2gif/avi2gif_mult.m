% This program converts AVI video file into an animated GIF file.
% The GIF format have advantages especially in Power-Point presentation,
% and in internet browsers.
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% September 2010 (C)

clear all

prompt = {'Enter number of videos to convert'};
name = 'Input AVI';
defans = {'1'};
options.Interpreter = 'tex';
ans = inputdlg(prompt,name,[1 40],defans,options);
if (isempty(ans{:}))
    return;
end
nvidz = str2num(ans{:});

ntot = 0;
for ivid = 1:nvidz
    [file_name file_path]=uigetfile({'*.avi','AVI video files'},'Select Video file');
    avi_info(ivid)=aviinfo([file_path,file_name]);
    ntot = ntot + avi_info(ivid).NumFrames;
    vidObj(ivid) = VideoReader([file_path,file_name]);
end
    
[file_name file_path]=uiputfile('*.gif','Save as animated GIF',[file_path,file_name(1:end-3)]);

lps=questdlg('How many loops?','Loops','Forever','None','Other','Forever');
switch lps
    case 'Forever'
        loops=65535;
    case 'None'
        loops=1;
    case 'Other'
        loops=inputdlg('Enter number of loops? (must be an integer between 1-65535)        .','Loops');
        loops=str2num(loops{1});
end
fps=avi_info(1).FramesPerSecond;
delay=inputdlg('What is the delay time? (in seconds)        .','Delay',1,{num2str(1/fps)});
delay=str2num(delay{1});
dly=questdlg('Different delay for the first image?','Delay','Yes','No','No');
if strcmp(dly,'Yes')
    delay1=inputdlg('What is the delay time for the first image? (in seconds)        .','Delay');
    delay1=str2num(delay1{1});
else
    delay1=delay;
end
dly=questdlg('Different delay for the last image?','Delay','Yes','No','No');
if strcmp(dly,'Yes')
    delay2=inputdlg('What is the delay time for the last image? (in seconds)        .','Delay');
    delay2=str2num(delay2{1});
else
    delay2=delay;
end

h = waitbar(0,['0% done'],'name','Progress') ;
cnt = 0;
for ivid = 1:nvidz
    nframes = avi_info(ivid).NumFrames;
    for i=1:nframes
        cnt = cnt+1;
        if strcmpi(avi_info(ivid).ImageType,'truecolor')
            %a=avi_file(i).cdata;
            a = readFrame(vidObj(ivid));
            [M  c_map]= rgb2ind(a,256);
        else
            %M=avi_file(i).cdata;
            M = readFrame(vidObj(ivid));
            error('Unable to decode AVI color map');
            c_map=avi_file(i).colormap;
        end
        if (i==1 && ivid==1)
            imwrite(M,c_map,[file_path,file_name],'gif','LoopCount',loops,'DelayTime',delay1)
        elseif (i==nframes && ivid==1)
            imwrite(M,c_map,[file_path,file_name],'gif','WriteMode','append','DelayTime',delay2)
        else
            imwrite(M,c_map,[file_path,file_name],'gif','WriteMode','append','DelayTime',delay)
        end
        waitbar(cnt/ntot,h,[num2str(round(100*cnt/ntot)),'% done']) ;
    end
end
close(h);
msgbox('Finished Successfully!')