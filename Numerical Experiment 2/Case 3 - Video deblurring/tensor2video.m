function tensor2video(A,label)
    % This function takes a three-dimensional tensor and converts it into a video
    video = VideoWriter(label, 'MPEG-4'); 
    open(video);
    numFrames=size(A,3);
    
    for i = 1:numFrames
        ith_frame = im2uint8(A(:,:,i)); 
        writeVideo(video, ith_frame); 
    end
    close(video);    
end