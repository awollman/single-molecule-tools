% A function recreating the output of Isabels extract_image_sequence_data
% but only for multipage tifs and outputing an array rather than a
% structure, modifed then from Adam's code of extract_image_sequence_dataAWarray
%If all=1, whole image is loaded, else loads from optional startFrame to endFrame
%Also loads folders of tif image sequences
%Now loads ascii files but currently only 512x512
function [numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = ExtractImageSequence3(image_label, all, startFrame, endFrame)
%image_label is the number before the file extension of the tif
filePath=dir(strcat('*',image_label,'*'));
if isempty(strfind(filePath.name,'tif'))==0
    tifImagePath0 = dir(strcat(image_label,'.tif'));
    
    % Error control if no .tif file can be found with that label.
    if isempty(tifImagePath0) % If there is no .tif image sequence file for such image_label, show error and exit function:
        error('Check you are in the correct directory and run again. No .tif file found for that image_label.');
    end
    
    image_path = tifImagePath0.name;
    InfoImage=imfinfo(image_path);
    frame_Ysize=InfoImage(1).Width;
    frame_Xsize=InfoImage(1).Height;
    numFrames=length(InfoImage);
    
    if all==0
        image_data=zeros(frame_Xsize,frame_Ysize,endFrame-startFrame+1,'uint16');
        for i=1:endFrame-startFrame+1
            image_data(:,:,i)=imread(image_path,i+startFrame-1);
        end;
    else
        image_data=zeros(frame_Xsize,frame_Ysize,numFrames,'uint16');
        for i=1:numFrames
            
            image_data(:,:,i)=imread(image_path,i);
        end;
    end
elseif isempty(strfind(filePath.name,'asc'))==0
    M=dlmread(filePath.name);
    frame_Xsize=length(unique(M(:,1)));
    frame_Ysize=size(M,2)-1;
    numFrames=size(M,1)/frame_Xsize;
      image_path=pwd;
    image_data=zeros(frame_Xsize,frame_Ysize,numFrames);
    for i=1:numFrames
       image_data(:,:,i)= M(1+(i-1)*512:i*512,2:end);
    end
    image_data(:,513,:)=[];
else
    if filePath.isdir==1
    
        cd(filePath.name)
        TifNames=dir('*.tif');
        InfoImage=imfinfo(TifNames(1).name);
        image_path=pwd;
        frame_Ysize=InfoImage(1).Width;
        frame_Xsize=InfoImage(1).Height;
        numFrames=length(TifNames(startFrame:endFrame));
        image_data=zeros(frame_Xsize,frame_Ysize,numFrames,'uint16');
        for k=1:numFrames
            image_data(:,:,k)=imread(TifNames(startFrame+k-1).name);
        end
        cd ..
    end
end

end