% [image_data, meta_data]=imEx1(fileName,startFrame,endFrame)

function [image_data, meta_data]=imEx1(fileName,startFrame,endFrame)
if exist('startFrame')==0
    data=bfopen(fileName);
    grayScale=data{1,1};
    [frame_Xsize, frame_Ysize]=size(grayScale{1,1});
    frame_Tsize=size(grayScale,1);
    image_data=zeros(frame_Xsize,frame_Ysize,frame_Tsize,'uint16');
    for i=1:frame_Tsize
        image_data(:,:,i)=grayScale{i,1};
    end
else
    InfoImage=imfinfo(fileName);
    frame_Ysize=InfoImage(1).Width;
    frame_Xsize=InfoImage(1).Height;
    numFrames=length(InfoImage);
    image_data=zeros(frame_Xsize,frame_Ysize,endFrame-startFrame+1,'uint16');
    for i=1:endFrame-startFrame+1
        image_data(:,:,i)=imread(fileName,i+startFrame-1);
    end;
end
MetaDataFile=strcat(fileName(1:end-7),'.ome.xml');
try
    meta_data=xmlread(MetaDataFile)
catch
    try
        meta_data=data{1,4};
    catch
         meta_data=0;
    end
end

end