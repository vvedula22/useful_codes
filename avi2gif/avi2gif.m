% This program converts AVI video file into an animated GIF file.
% The GIF format have advantages especially in Power-Point presentation,
% and in internet browsers.
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% September 2010 (C)

clear all
[file_name file_path]=uigetfile({'*.avi','AVI video files'},'Select Video file');
[file_name2 file_path2]=uiputfile('*.gif','Save as animated GIF',[file_path,file_name(1:end-3)]);
avi_info=aviinfo([file_path,file_name]);
vidObj = VideoReader([file_path,file_name]);
%avi_file=aviread([file_path,file_name]);

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
fps=avi_info.FramesPerSecond;
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
%nframes = length(avi_file);
nframes = avi_info.NumFrames;
for i=1:nframes
    if strcmpi(avi_info.ImageType,'truecolor')
        %a=avi_file(i).cdata;
        a = readFrame(vidObj);
        [M  c_map]= rgb2ind(a,256);
    else
        %M=avi_file(i).cdata;
        M = readFrame(vidObj);
        c_map=avi_file(i).colormap;
    end
    if i==1
        imwrite(M,c_map,[file_path2,file_name2],'gif','LoopCount',loops,'DelayTime',delay1)
    elseif i==nframes
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay2)
    else
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
    end
    waitbar(i/nframes,h,[num2str(round(100*i/nframes)),'% done']) ;
end
close(h);
msgbox('Finished Successfully!')